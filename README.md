# Assignment A2

Please carefully read the [rules](rules/README.md). After reading this document, copy the [answer sheet](answer-sheet.md) to file `report.md` and provide your answers in the latter. Do not edit this document nor the [answer sheet](answer-sheet.md).

## Objectives

In this assignment, you are tasked to implement numerical integration. You must not use any ODE library, numerical integration routine, or any other implementation of an integrator that is not authored by you.

## Tasks

### 1 Update your orbit

1. Go to [celestrak](https://celestrak.org/satcat/search.php) and update TLE data of the Earth-orbiting objects you chose in assignment A0.

2. Use your SGP4 propagator to generate an orbit with at least 5 revolutions. We'll call this the *sgp4 orbit*.

### 2 Implement an Euler integrator

1. Define a function that computes the 3D gravitational acceleration resulting from a point mass as function of position.
2. Build an Euler integrator that receives:
    - initial state vector: 3D position and velocity,
    - time step (scalar),
    - stop time (scalar), and
    - the function from step 1.
3. Reproduce the orbit from [Task 1](1-Update-your-orbit), in the form of the henceforth called *integrated orbit*:
    - for initial state vector, use the first data point of the sgp4 orbit (given in Cartesian coordinates);
    - use the same start/stop/step time as in Task 1.
4. Measure the time your integrator takes to produce the orbit, as well as the average computational time per integration step.

**Opportunity for Excellence**: Implement additional integrators and the ability to easily switch between them.

### 3 Assess the errors of your integrator

1. Show in the same plot the local error as function of time for two different step sizes. You need to re-generate the sgp4 orbit (from the TLE data) with this new time domain.
2. Plot the global error (the RMS of the local error) as function of different step sizes used in your Euler integrator. Consider at least 10 different values of step size. 

**Opportunity for Excellence**: Repeat this task for one other integrator that you implemented in [Task 2](2-Implement-an-Euler-integrator). Derive observations, interpretations and conclusions on the differences between the two integrators.

### 4 Assess the stability of your integrator

For each one of the initial state vector components (both position and velocity) of the spg4 orbit derived in [Task 1](1-Update-your-orbit):

1. Change it by a very small amount and re-compute the resulting orbit, which we call *tweaked orbit*.
2. Plot the 3D different of the **position** of the tweaked orbit with respect to the spg4 orbit as function of time. 
3. Repeat steps 1-2 above for 5 others components of initial state vector, including velocity. Keep these values very close to the original one, and experiment to see how far can you go before the result diverges significantly while still showing differences that are visible in the plot.
4. Show all 6 lines in the same plot. Don't forget to include a legend that properly identifies each line. You should have a single plot with 6 lines. 
5. Show the global 3D difference (the RMS of the 3D different as function of time) as function of the value of the tweaked initial state vector component, in terms of percentage relative to the original value. This produces 6 plots, with 1 line each.

**Opportunity for Excellence**: Repeat this task for one other integrator that you implemented in [Task 2](2-Implement-an-Euler-integrator). Derive observations, interpretations and conclusions on the differences between the two integrators.

## Final remarks

Don't forget to:

- report Assignment and Code Excellence
- report the time it took you to solve the assignment
- provide feedback on the assignment, namely how interesting you found it and how challenging it was to answer it.
