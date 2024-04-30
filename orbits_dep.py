import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches
from math import sqrt
from operator import add

GRAV_CONST = 6.67430e-11

COLORS = ["b", "g", "r", "c", "m", "y"]


class CelestialBody:
    def __init__(
        self,
        mass=1,
        r_i: npt.NDArray = None,
        v_i: npt.NDArray = None,
        a_i: npt.NDArray = None,
        dimension=2,
        name=None,
    ):
        """Initiate a celestial body,

        Args:
            mass (_type_): mass
            r_i (npt.NDArray, optional): initial position. Defaults to None.
            v_i (npt.NDArray, optional): initial velocity. Defaults to None.
            a_i (npt.NDArray, optional): initial acceleration. Defaults to None.
            dimension (int, optional): dimension of the coordinates. Defaults to 2.
        """
        self.dimension = dimension
        if r_i is None:
            r_i = np.zeros(self.dimension)
        else:
            if r_i.size != self.dimension:
                raise ValueError(
                    f"Dimension parameter: {self.dimension}, Position's dimension {r_i.shape}, they should match"
                )

        if v_i is None:
            v_i = np.zeros(self.dimension)
        else:
            if v_i.size != self.dimension:
                raise ValueError(
                    f"Dimension parameter: {self.dimension}, Velocity's dimension {v_i.shape}, they should match"
                )

        if a_i is None:
            a_i = np.zeros(self.dimension)
        else:
            if a_i.size != self.dimension:
                raise ValueError(
                    f"Dimension parameter: {self.dimension}, Acceleration's dimension {a_i.shape}, they should match"
                )
        self.name = name
        self.r = r_i
        self.v = v_i
        self.a = a_i
        self.mass = mass
        self.color = np.random.choice(COLORS)

    @staticmethod
    def init_moon(other_body, moon, orbital_radius):
        m = moon
        if np.linalg.norm(other_body.r) == 0:
            pos = np.zeros(other_body.dimension)
            pos[0] = orbital_radius
            m.r = pos
        else:
            other_to_zero = (
                orbital_radius * (-other_body.r) / np.linalg.norm(other_to_zero)
            )
            m.r = other_to_zero
        return m

    def grav_force_from(self, other):
        me_to_them = other.r - self.r
        if np.linalg.norm(me_to_them) == 0:
            return np.zeros(self.dimension)

        force = (
            GRAV_CONST * self.mass * other.mass / (np.linalg.norm(me_to_them) ** 2)
        ) * (me_to_them / np.linalg.norm(me_to_them))
        return force

    # def potential_energy(self, other):


class CelestialSim:
    def __init__(
        self, bodies, timestep=1, integration_method="euler-cromer", strict=True
    ):
        """_summary_

        Args:
            bodies (_type_): List of celestial bodies
            timestep (int, optional): In seconds. Defaults to 1.
            integration_method (str, optional): Integration method. Defaults to 'euler-cromer'.
        """
        self.time = 0
        self.bodies = np.array(bodies)
        self.timestep = timestep
        self.strict = strict
        self.max_x = 1
        self.max_y = 1
        self.min_x = -1
        self.min_y = -1
        self.max_mass = 1
        if integration_method == "euler-cromer":
            self.integrator = EulerCromer(timestep=timestep)

        if len(bodies) != 0:
            self.dim = bodies[0].dimension
            for index, b in enumerate(bodies):
                try:
                    bodies[index] = self.fix_body(b, strict)
                except:
                    raise TypeError(
                        f"Body: {b.name}, index: {index}, has incompatible dimension {b.dimension}"
                    )
                if bodies[index].mass > self.max_mass:
                    self.max_mass = bodies[index].mass
                if bodies[index].r[0] > self.max_x:
                    self.max_x = bodies[index].r[0]
                if bodies[index].r[1] > self.max_y:
                    self.max_y = bodies[index].r[1]
                if bodies[index].r[0] < self.min_x:
                    self.min_x = bodies[index].r[0]
                if bodies[index].r[1] < self.max_y:
                    self.min_y = bodies[index].r[1]

    @staticmethod
    def init_from_file(filepath):
        """file should be structured as follows:
        # = comment (will be ignored)
        ; is separator
        1 celestial body per line with data as follows:
        name;mass;v_vec;a_vec;r_vec;moon_of;moon_radius

        if a value is empty, it will treat it like None

        vec should be (x1,x2,...,xn) up to n=dimension

        Args:
            filepath (_type_): _description_
        """
        read_bodies = []
        with open(filepath, "r") as body_data:
            for index, line in enumerate(body_data):
                if line.startswith("#"):
                    continue
                l_data = line.split(";")
                if len(l_data != 7):
                    raise ValueError(
                        f"Line {index} contains {len(l_data)} arguments, which is invalid, should be 7. Also make sure to use ; as separator"
                    )
                name = l_data[0]
                mass = float(l_data[1])
                r = None
                v = None
                a = None
                moon_of = None
                if l_data[2]:
                    v = np.array([float(num) for num in l_data[2][1:-1].split(",")])
                if l_data[3]:
                    a = np.array([float(num) for num in l_data[3][1:-1].split(",")])
                if l_data[4]:
                    r = np.array([float(num) for num in l_data[4][1:-1].split(",")])
                if l_data[5]:  # Is a moon of something
                    moon_of = l_data[5]
                    if l_data[6]:
                        moon_rad = l_data[6]
                    else:
                        raise ValueError(f"Line {index} is a moon but no radius given")
                new_body = CelestialBody(
                    name=name,
                    mass=mass,
                    r_i=r,
                    v_i=v,
                    a_i=a,
                )

        pass

    def fix_body(self, body, strict):
        b = body
        if b.dimension != self.dim:
            if strict:
                # raise TypeError(
                #     f"Body: {b.name}, index: {index}, has incompatible dimension {b.dimension}"
                # )
                raise TypeError(
                    f"Body: {b.name}, has incompatible dimension {b.dimension}, should be {self.dim}"
                )
            else:
                if b.dimension > self.dim:
                    b.r = b.r[: self.dim]
                    b.v = b.v[: self.dim]
                    b.a = b.a[: self.dim]
                else:
                    b.r = np.pad(b.r, (0, self.dim - b.dimension), constant_values=0)
        return b

    def add_body(self, body, strict=None):
        s = None
        if strict is None:
            s = self.strict
        self.bodies.append(self.fix_body(body, strict=s))

    def grav_force_on(self, body):
        force_vector = np.zeros(self.dim)
        for b in self.bodies:
            force_vector += body.grav_force_from(b)
        return force_vector

    def step_hidden(self):
        for b in self.bodies:
            grav_force = self.grav_force_on(b)
            acceleration = grav_force / b.mass
            b.a = acceleration
            # print(
            #     f"Name: {b.name}, position_before: {b.r}, position_after: {b.r}, force: {grav_force}"
            # )
            b.r = self.integrator.r_next(b.r, b.v, b.a)
            b.v = self.integrator.v_next(b.v, b.a)
        self.time += self.timestep
        pass

    def run_sim_infinite(self, steps, show_energy=False):
        fig, ax = plt.subplots(1, 1)
        ax.axis("equal")
        ax.set(
            xlim=(5 * self.min_x, 5 * self.max_x),
            ylim=(5 * self.min_y, 5 * self.max_y),
        )
        pe_data = []
        ke_data = []
        te_data = []

        def update_infinite(frame):
            # for p in ax.patches:
            #     p.remove()
            container = self.step_graphical()
            for p in container:
                ax.add_patch(p)

            pe_data.append(self.total_potential_energy())
            ke_data.append(self.total_kinetic_energy())
            te_data.append(self.total_potential_energy() + self.total_kinetic_energy())

            # if (frame + 1) % 10 == 0:
            #     print(f"Total Kinetic Energy: {self.total_kinetic_energy()}")
            #     print(f"Total Poential Energy: {self.total_potential_energy()}")
            #     print(
            #         f"Total Energy: {self.total_kinetic_energy() + self.total_potential_energy()}"
            #     )
            #     pass
            return container

        ani = animation.FuncAnimation(
            fig=fig, func=update_infinite, frames=steps, interval=50, blit=True
        )
        plt.show()
        return (ke_data, pe_data, te_data)

    def run_sim_finite(self, steps, show):
        planet_data = []
        ke_data = []
        pe_data = []
        te_data = []
        for i in range(steps):
            if show:
                planet_data.append((self.planet_patches()))
            ke_data.append(self.total_kinetic_energy())
            pe_data.append(self.total_potential_energy())
            te_data.append(self.total_kinetic_energy() + self.total_potential_energy())
            self.step_hidden()
        if show:
            return (planet_data, ke_data, pe_data, te_data)
        return (None, ke_data, pe_data, te_data)

    def run_sim(self, steps=100, finite=True, show=True):
        ke_data, pe_data, te_data = None, None, None
        if finite:
            # planet_data = []
            # ke_data = []
            # pe_data = []
            # te_data = []
            # for i in range(steps):
            #     if show:
            #         planet_data.append((self.planet_patches()))
            #     ke_data.append(self.total_kinetic_energy())
            #     pe_data.append(self.total_potential_energy())
            #     te_data.append(
            #         self.total_kinetic_energy() + self.total_potential_energy()
            #     )
            #     self.step_hidden()
            planet_data, ke_data, pe_data, te_data = self.run_sim_finite(steps, show)
            if show:
                self.plot_finite(planet_data)

            # fig, ax = plt.subplots()
            # ax.axis("equal")
            # ax.set(
            #     xlim=(5 * (self.min_x + 1), 5 * (self.max_x + 1)),
            #     ylim=(5 * (self.min_y + 1), 5 * (self.max_y + 1)),
            # )

            # def update(frame):
            #     container = PatchCollection(planet_data[frame])
            #     if (frame + 1) % 10 == 0:
            #         print(f"Total Kinetic Energy: {self.total_kinetic_energy()}")
            #         print(f"Total Poential Energy: {self.total_potential_energy()}")
            #         print(
            #             f"Total Energy: {self.total_kinetic_energy() + self.total_potential_energy()}"
            #         )
            #         pass
            #     for p in ax.collections:
            #         if isinstance(p, PatchCollection):
            #             p.remove()
            #     return ax.add_collection(container)

            # ani = animation.FuncAnimation(
            #     fig=fig, func=update, frames=steps, interval=20
            # )
            # ani.save(filename="planetmovie.mp4")
            # plt.show()

        else:
            ke_data, pe_data, te_data = self.run_sim_infinite(steps=steps)

        self.plot_energy(ke_data, pe_data, te_data)

    def plot_finite(self, planet_data):
        fig, ax = plt.subplots()
        ax.axis("equal")
        ax.set(
            xlim=(5 * (self.min_x + 1), 5 * (self.max_x + 1)),
            ylim=(5 * (self.min_y + 1), 5 * (self.max_y + 1)),
        )

        def update(frame):
            container = PatchCollection(planet_data[frame])
            for p in ax.collections:
                if isinstance(p, PatchCollection):
                    p.remove()
            return ax.add_collection(container)

        ani = animation.FuncAnimation(
            fig=fig, func=update, frames=len(planet_data), interval=20
        )
        plt.show()

        ani.save(filename="planetmovie.mp4")

    def plot_energy(self, ke_data, pe_data, te_data):
        fig, ax = plt.subplots()
        ax.plot(ke_data, color="r")
        ax.plot(pe_data, color="g")
        ax.plot(te_data, color="b")
        plt.show()

    # def plot_finite(self, planet_history):
    #     fig, ax = plt.subplots()
    #     ax.set(
    #         xlim=(2 * self.min_x, 2 * self.max_x),
    #         ylim=(2 * self.min_y, 2 * self.max_y),
    #     )

    #     def update(frame):
    #         for p in ax.collections:
    #             if isinstance(p, PatchCollection):
    #                 p.remove()
    #         print(self.total_kinetic_energy())
    #         return ax.add_collection(planet_history[frame])

    #     ani = animation.FuncAnimation(
    #         fig=fig, func=update, frames=len(planet_history), interval=10
    #     )
    #     ani.save(filename="planetmovie.mp4")
    #     plt.show()

    def total_kinetic_energy(self):
        k_e = 0
        for b in self.bodies:
            k_e += (1 / 2) * b.mass * (np.linalg.norm(b.v) ** 2)
        return k_e

    def total_potential_energy(self):
        p_e = 0
        for b in self.bodies:
            for other in self.bodies:
                dist = np.linalg.norm(other.r - b.r)
                if dist == 0:
                    p_e += 0
                else:
                    p_e += (((-1) * GRAV_CONST * b.mass * other.mass) / dist) / len(
                        self.bodies
                    )
        return p_e

    def step_graphical(self, show=False):
        planets = self.planet_patches()
        self.step_hidden()
        return planets

    def planet_patches(self):
        planets = []

        for b in self.bodies:
            planets.append(
                mpatches.CirclePolygon(
                    (b.r[0], b.r[1]),
                    radius=self.mass_to_radius(b.mass),
                    facecolor=b.color,
                    edgecolor="k",
                    alpha=0.6,
                )
            )
        return planets

    def mass_to_radius(self, mass):
        scale = ((1 + self.max_x - self.min_x) + (1 + self.max_y - self.min_y)) / (
            4 * 5
        )
        return max((mass / self.max_mass) * scale, scale * 0.2)


class EulerCromer:
    def __init__(self, timestep) -> None:
        self.timestep = timestep
        pass

    def v_next(self, v_now, a_now):
        return v_now + a_now * self.timestep

    def r_next(self, r_now, v_now, a_now):
        return r_now + self.v_next(v_now, a_now) * self.timestep
