# Celestial Body Simulation Project

This project is a Python-based simulation of celestial bodies in a 2D space. It uses principles of physics, specifically Newton's law of universal gravitation, to model the motion of celestial bodies over time.

## Key Concepts

### Gravitational Constant

The gravitational constant (GRAV_CONST) is a physical constant involved in the calculations of gravitational force between two bodies. It is approximately equal to 6.67430e-11 m³ kg⁻¹ s⁻².

### Euler-Cromer Method

The Euler-Cromer method, also known as the semi-implicit Euler method, is a numerical method used to solve ordinary differential equations. It's used in this project to calculate the next position and velocity of a celestial body.

### Celestial Body

A celestial body in this project is represented by the CelestialBody class. It has properties like mass, position (r), velocity (v), acceleration (a), and dimension. It also has methods to calculate gravitational force from another body.

### Celestial Simulation

The CelestialSim class represents a simulation of multiple celestial bodies. It uses an integrator (like Euler-Cromer) to calculate the motion of the bodies over time. It also has methods to add a body to the simulation, calculate gravitational force on a body, and run the simulation.

## Academic Relevance

This project is a practical application of physics and numerical methods in a computer simulation. It can be used as a teaching tool for understanding the principles of celestial mechanics and the effects of gravity on the motion of celestial bodies. It also demonstrates the use of numerical methods, like the Euler-Cromer method, in solving differential equations that represent physical phenomena.

The project also touches on concepts from computer science and software engineering, such as object-oriented programming, error handling, and the use of libraries for mathematical and graphical functions.

## Usage

To use this project, you can create instances of CelestialBody and CelestialSim, and use the methods provided by these classes to set up and run a simulation. You can also customize the simulation by changing the properties of the celestial bodies and the parameters of the simulation.

## Future Work

Potential improvements to this project could include adding support for 3D simulations, improving the accuracy of the numerical methods, and adding more features to the simulation, such as collisions between bodies or the effects of non-gravitational forces.
