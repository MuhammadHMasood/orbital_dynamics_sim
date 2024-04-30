# Muhammad Masood
# s2410617


import numpy as np
import numpy.typing as npt
import matplotlib.text as mtext
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches
import matplotlib.transforms as mtrans
from math import sqrt, pi
from math import pow
from operator import add

GRAV_CONST = 6.67430e-11

COLORS = ["b", "g", "r", "c", "m", "y"]

DENSITY = 5515


class EulerCromer:
    """Euler-Cromer integrator"""

    def __init__(self, timestep=1) -> None:
        self.timestep = timestep
        pass

    def v_next(self, v_now, a_now):
        return v_now + a_now * self.timestep

    def r_next(self, r_now, v_now, a_now):
        return r_now + self.v_next(v_now, a_now) * self.timestep


# class SecondOrder:
#     def __init__(self, timestep) -> None:
#         self.timestep = timestep
#         pass

#     def v_next(self, v_now, a_now):
#         return v_now + a_now * self.timestep

#     def r_next(self, r_now, v_now, a_now):
#         return r_now + self.v_next(v_now, a_now) * self.timestep


# r(t+dt) = r(t) + v(t)dt + 1/2 a(t)(dt^2)
# derive a(t+dt) but don't set it
# a(t+dt) =
#
#


class CelestialBody:
    def __init__(
        self,
        mass=1,
        r_i: npt.NDArray = None,
        v_i: npt.NDArray = None,
        a_i: npt.NDArray = None,
        dimension=2,
        name=None,
        radius=None,
        density=DENSITY,
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
        if radius == None:
            volume = mass / DENSITY
            self.radius = pow((3 / (4 * pi)) * volume, 1 / 3)
        else:
            self.radius = radius
            # v = 4/3 pi r^3
            # (3 / (4 * pi) * v)^(1/3)

    @staticmethod
    def init_moon(other_body, moon, orbital_radius):
        """Generates a moon to another body given an orbital radious

        Args:
            other_body (_type_): The body to be the moon of
            moon (_type_): some specifications of the moon
            orbital_radius (_type_): the orbital radius of the moon

        Returns:
            _type_: _description_
        """
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
        if np.linalg.norm(m.v) == 0:  # if it has
            vel = np.zeros(other_body.dimension)
            vel[-1] = sqrt(GRAV_CONST * other_body.mass / orbital_radius)
            m.v = vel

        return m

    def grav_force_from(self, other, no_explosion=False, fun=False):
        """Returns the gravitational force from another body onto self.

        Args:
            other (_type_): _description_
            no_explosion (bool, optional): No explosion tries to fix the issue of explosions when radius << 1. Defaults to False.
            fun (bool, optional): Fun chooses whether or not to have accurate radii. Defaultss to false

        Returns:
            _type_: _description_
        """
        me_to_them = other.r - self.r
        if np.linalg.norm(me_to_them) == 0:
            return np.zeros(self.dimension)
        if fun:
            tol = 0.5
        else:
            tol = self.radius
        if no_explosion and np.linalg.norm(me_to_them) < tol:
            pass
            me_to_them = self.radius
        # print(f"Apparent Distance: {me_to_them}, Calculated Radius: {self.radius}")

        force = (
            GRAV_CONST * self.mass * other.mass / (np.linalg.norm(me_to_them) ** 2)
        ) * (me_to_them / np.linalg.norm(me_to_them))
        return force

    # def potential_energy(self, other):


class CelestialSim:
    def __init__(
        self,
        bodies,
        timestep=1,
        integrator=EulerCromer,
        strict=True,
        no_explosion=False,
        fun=True,
        steps=100,
    ):
        """_summary_

        Args:
            bodies (_type_): List of bodies
            timestep (int, optional): length of a timestep. Defaults to 1.
            integrator (_type_, optional): The integrator. Defaults to EulerCromer.
            strict (bool, optional): Whether to reject bodies with incorrect dimensions. Defaults to True.
            no_explosion (bool, optional): Whether to try and fix explosions when bodies are close. Defaults to False.
            fun (bool, optional): Accurate radius or not. Defaults to True.
            steps (int, optional): Number of steps if finite simulation. Defaults to 100.

        Raises:
            TypeError: _description_
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
        self.no_explosion = no_explosion
        self.fun = fun
        self.integrator = integrator(timestep=timestep)
        self.steps = steps

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

        for steps and timestep, do !steps;timestep

        if a value is empty, it will treat it like None

        vec should be (x1,x2,...,xn) up to n=dimension

        Args:
            filepath (_type_): _description_
        """
        read_bodies = []
        steps = 0
        step_size = 1
        with open(filepath, "r") as body_data:
            for index, line in enumerate(body_data):
                if line.startswith("!"):
                    steps = int(line.strip()[1:].split(";")[0])

                    step_size = float(line.strip()[1:].split(";")[1])

                    continue
                if line.startswith("#"):
                    continue
                l_data = line.strip().split(";")
                if len(l_data) != 7:
                    raise ValueError(
                        f"Line {index} contains {len(l_data)} arguments, which is invalid, should be 7. Also make sure to use ; as separator"
                    )
                name = l_data[0]
                mass = float(l_data[1])
                r = None
                v = None
                a = None
                moon_of = None
                moon_rad = None
                if l_data[5]:  # Check if moon, equivalent to if l_data[5] != ""
                    moon_of = l_data[5]
                    if l_data[6]:
                        moon_rad = float(l_data[6])
                    else:
                        raise ValueError(f"Line {index} is a moon but no radius given")
                elif l_data[4]:  # If not a moon, get position
                    r = np.array([float(num) for num in l_data[4][1:-1].split(",")])
                else:
                    raise ValueError(
                        f"Line {index} is not a moon and not given a position"
                    )

                if l_data[2]:
                    v = np.array([float(num) for num in l_data[2][1:-1].split(",")])
                if l_data[3]:
                    a = np.array([float(num) for num in l_data[3][1:-1].split(",")])
                if moon_of:
                    if len(read_bodies) == 0:
                        raise ValueError(
                            "You can't have a moon as your first object because moons are defined in relation to other celestial bodies"
                        )
                    body_of = next(body for body in read_bodies if body.name == moon_of)
                    moon = CelestialBody(
                        name=name,
                        mass=mass,
                        r_i=r,
                        v_i=v,
                        a_i=a,
                        dimension=body_of.dimension,
                    )
                    new_body = CelestialBody.init_moon(
                        other_body=body_of, moon=moon, orbital_radius=moon_rad
                    )
                    pass

                else:
                    new_body = CelestialBody(
                        name=name, mass=mass, r_i=r, v_i=v, a_i=a, dimension=len(r)
                    )

                read_bodies.append(new_body)
        pass

        return CelestialSim(bodies=read_bodies, timestep=step_size, steps=steps)

    def fix_body(self, body, strict):
        """fix the dimensions of a body with incorrect dimensions

        Args:
            body (_type_): _description_
            strict (_type_): _description_

        Raises:
            TypeError: _description_

        Returns:
            _type_: _description_
        """
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
        """Add a body to the sim

        Args:
            body (_type_): _description_
            strict (_type_, optional): _description_. Defaults to None.
        """
        s = None
        if strict is None:
            s = self.strict
        self.bodies.append(self.fix_body(body, strict=s))

    def grav_force_on(self, body):
        force_vector = np.zeros(self.dim)
        for b in self.bodies:
            force_vector += body.grav_force_from(
                b, no_explosion=self.no_explosion, fun=self.fun
            )
        return force_vector

    def step_hidden(self):
        for b in self.bodies:
            grav_force = self.grav_force_on(b)
            acceleration = grav_force / b.mass
            b.a = acceleration
            b.r = self.integrator.r_next(b.r, b.v, b.a)
            b.v = self.integrator.v_next(b.v, b.a)
        self.time += self.timestep

    def run_sim(self, steps=100, finite=True, show=True):
        """run the simulation

        Args:
            steps (int, optional): _description_. Defaults to 100.
            finite (bool, optional): _description_. Defaults to True.
            show (bool, optional): _description_. Defaults to True.
        """
        ke_data, pe_data, te_data = None, None, None
        if finite:
            planet_data, ke_data, pe_data, te_data = self.run_sim_finite(steps, show)
            if show:
                self.plot_finite(planet_data)
        else:
            ke_data, pe_data, te_data = self.run_sim_infinite(steps)
        self.plot_energy(ke_data, pe_data, te_data)

    def run_sim_infinite(self, steps):
        fig, ax = plt.subplots(1, 1)
        ax.axis("equal")
        ax.set(
            xlim=(
                5 * min(self.min_y, self.min_x, -self.max_x, -self.max_y),
                5 * max(self.max_y, self.max_x, -self.min_y, -self.min_x),
            ),
            ylim=(
                5 * min(self.min_y, self.min_x, -self.max_x, -self.max_y),
                5 * max(self.max_y, self.max_x, -self.min_y, -self.min_x),
            ),
        )
        print(self.max_y, self.max_x, -self.min_y, -self.min_x)
        pe_data = []
        ke_data = []
        te_data = []  #

        ke_patch = mpatches.Patch(color="r", label=f"E_kinetic: {0}")
        pe_patch = mpatches.Patch(color="g", label=f"E_potential: {0}")
        te_patch = mpatches.Patch(color="b", label=f"E_total {0}")

        def update_infinite(frame):
            for p in ax.patches:
                p.remove()

            container = self.step_graphical()
            for p in container:
                ax.add_patch(p)

            pe_data.append(self.total_potential_energy())
            ke_data.append(self.total_kinetic_energy())
            te_data.append(self.total_energy())

            if (frame + 1) % 10 == 0:
                ke_patch.set(label=f"E_kinetic: {self.total_kinetic_energy():.3e}")
                pe_patch.set(label=f"E_potential: {self.total_potential_energy():.3e}")
                te_patch.set(label=f"E_total: {(self.total_energy()):0.3e}")
            L = plt.legend(
                handles=[ke_patch, pe_patch, te_patch],
                prop={"size": 6},
                framealpha=0.5,
            )
            container.append(L)
            return container

        ani = animation.FuncAnimation(
            fig=fig, func=update_infinite, frames=21, interval=5, blit=True
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
            te_data.append(self.total_energy())
            self.step_hidden()
        if show:
            return (planet_data, ke_data, pe_data, te_data)
        return (None, ke_data, pe_data, te_data)

    def plot_finite(self, planet_data):
        fig, ax = plt.subplots()
        ax.axis("equal")
        ax.set(
            xlim=(
                5 * min(self.min_y, self.min_x, -self.max_x, -self.max_y),
                5 * max(self.max_y, self.max_x, -self.min_y, -self.min_x),
            ),
            ylim=(
                5 * min(self.min_y, self.min_x, -self.max_x, -self.max_y),
                5 * max(self.max_y, self.max_x, -self.min_y, -self.min_x),
            ),
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

        # ani.save(filename="planetmovie.mp4")

    def plot_energy(self, ke_data, pe_data, te_data):
        fig, ax = plt.subplots()
        tr = mtrans.offset_copy(ax.transData, fig=fig, x=0.0, y=-1.5, units="points")
        ax.plot(ke_data, color="r", alpha=0.5)
        ax.plot(pe_data, color="g", alpha=0.5)
        ax.plot(te_data, color="b", alpha=0.5)
        red_patch = mpatches.Patch(color="r", label="Kinetic Energy")
        green_patch = mpatches.Patch(color="g", label="Potential Energy")
        blue_patch = mpatches.Patch(color="b", label="Total Energy")
        ax.legend(handles=[red_patch, green_patch, blue_patch])
        ax.set_xlabel(f"Time, in intervals of {self.timestep} seconds")
        ax.set_ylabel(f"Energy (Joules)")
        ax.title.set_text("Energy Over Time")

        plt.show()

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

    def total_energy(self):
        return self.total_kinetic_energy() + self.total_potential_energy()

    def step_graphical(self):
        planets = self.planet_patches()
        self.step_hidden()
        return planets

    def planet_patches(self):
        planets = []

        for b in self.bodies:
            # print(f"({b.r[0]}, {b.r[1]})")
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
        if self.fun:
            scale = ((1 + self.max_x - self.min_x) + (1 + self.max_y - self.min_y)) / (
                4 * 5
            )

            return max((mass / self.max_mass) * scale, scale * 0.5)
        else:
            return pow((3 / (4 * pi)) * (mass / DENSITY), 1 / 3)

        scale = ((1 + self.max_x - self.min_x) + (1 + self.max_y - self.min_y)) / (
            4 * 5
        )
        return max((mass / self.max_mass) * scale, scale * 0.2)
