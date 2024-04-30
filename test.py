from orbits import *

body_1 = CelestialBody(
    mass=3e10,
    r_i=np.array((0, 0)),
    v_i=np.array((0, 0)),
    dimension=2,
    name="Planet_1",
)

body_2 = CelestialBody(
    mass=1e8,
    r_i=np.array((-5, -7)),
    v_i=np.array((0.3, 0)),
    a_i=np.array((0, 0)),
    dimension=2,
    name="Planet_2",
)

body_3 = CelestialBody(
    mass=1e8,
    r_i=np.array((5, 7)),
    v_i=np.array((-0.3, 0.2)),
    dimension=2,
    name="Planet_3",
)


sim = CelestialSim(
    [body_1, body_2, body_3], timestep=0.05, no_explosion=False, fun=True
)

# sim = CelestialSim.init_from_file("./test.txt")
# sim.fun = True

sim.run_sim(steps=sim.steps, finite=False)
