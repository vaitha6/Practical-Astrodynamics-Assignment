# Assignment A4

Please carefully read the [rules](rules/README.md). After reading this document, copy the [answer sheet](answer-sheet.md) to file `report.md` and provide your answers in the latter. Do not edit this document nor the [answer sheet](answer-sheet.md).

## Objectives

In this assignment, you are tasked to improve your orbit prediction by optimizing the initial state vector, so that you get a better orbit fit to the SGP4 orbit. 

## Tasks

### 1 Update your orbit

1. Go to [celestrak](https://celestrak.org/satcat/search.php) and update TLE data of the Earth-orbiting objects you chose in assignment A0.

2. Generate the sgp4 orbit, similarly to Step 2 of Task 1 of Assignment A2 (see below details for the start, step and stop times). You may include any (or none) of the force model components you developed in Assignment A3.

3. Reproduce the orbit from step 2, in the form of the *integrated orbit*, similarly to Step 3 of Task 2 of Assignment A2, using your Euler integrator:

  - for initial state vector $`y_0`$, use the first data point of the sgp4 orbit (given in Cartesian coordinates);
  - use the same start/stop/step time as in Step 2.

Choose the parameters that define the time domain, i.e., start, step and stop times, such that each complete orbit integration (using your own integrator):

  - does not take longer than half a second of [wall-clock time](https://pythonspeed.com/articles/blocking-cpu-or-io/), with smaller execution times preferable,
  - spans over (at least) 10000 integration steps, and
  - produces an integrated orbit with a RMS of the epoch-by-epoch 3D position difference relative to the sgp4 orbit (a scalar quantity) that is at the level of a few km.

If you are not able to make this benchmark, there may be unnecessary code in the integration kernel, i.e., the part of the code that computes $`y_{i+1}`$. Unnecessary code are `if`, `print` and other non-mathematical statements. Also, try different techniques to save the state vector at each integration step: you may append a onto growing matrix (using `np.append` or `np.column_stack`, see [this](https://stackoverflow.com/a/65470570) or [this](https://www.geeksforgeeks.org/python-ways-to-add-row-columns-in-numpy-array/), for example) or directly assign the columns of a pre-sized matrix (possibly there are others techniques). Spend some time streamlining your integration function, it will pay back later.

As suggested in the requirements above, you are free to make the step/stop time as small as possible. Nevertheless, choose these parameters so that you have an orbit with a few km 3D RMS difference. This is needed in order to avoid that the search space of the optimization is too flat, i.e., you have something substantial to optimize. If you choose stop time of 1 second (assuming start time of zero) and a time step of $`10^{-5}`$, your orbit will move roughly 7 km, and your 3D position difference will be on the order of centimetres.

Note that if you do a grid search optimization, with 3 values (which is not a lot) for each one of the 6 unknowns (3 positions and 3 velocities), your search space has a size of $`3^6=729`$ and each optimization run will take roughly 6 minutes (assuming half a second per orbit integration). If you build your grid with 5 values for each unknown, your optimization may take over two hours to finish. It is not a bad idea to adapt the values of the grid manually after each optimization run, or implement a more advanced optimization algorithm.

**Opportunity for Excellence**: Assess the effect of the integration time step (at least 10 different values) to the resulting 3D RMS difference and integration wall-clock time, by plotting these two quantities in the same plot (using two y-axis).

### 2 Optimize the initial state vector of your orbit

Based on (one of) the integrators you implemented in assignment A2:

1. Build an objective function $`f(\vec x)`$ that computes the global position error of an orbit integrated from a set of initial state vector differences $`\vec x = \delta \vec y_0`$ that change the original one $`y_0`$ by the small amount $`\delta \vec y_0`$. In other words, the input for this objective function is $`\delta \vec y_0`$ and the output is the 3D RMS position difference, represented by $`\Upsilon^{({\rm pos})}`$:

$`\Upsilon^{({\rm pos})}=f(\delta \vec y_0)`$

Recall that $`\Upsilon^{({\rm pos})}`$ is computed as the RMS of the magnitude of the vector difference between a) the position of the integrated orbit resulting from the initial state vector $`\vec y_0+\delta \vec y_0`$ and b) the position of the sgp4 orbit.

2. Minimize $`f`$ using any of the optimization algorithms discussed during the lectures, thus producing the *optimized orbit*. 

3. Plot the x,y,z position and velocity coordinates of the sgp4 and optimized orbits, along with their difference in a second y-axis. You should have 6 plots, one for each position and velocity coordinate, with 3 lines each: 

  - one line for the sgp4 orbit, 
  - one line for the optimized orbit and 
  - one line for their difference (connected to a second y-axis).

**Opportunity for Excellence**: Implement a second optimization algorithm and repeat steps 2 and 3 of this task.

**Opportunity for Excellence**: Use a different orbit integrator with (any of) the optimization algorithm(s) you implemented and repeat steps 2 and 3 of this task.

## Final remarks

Don't forget to:

- report Assignment and Code Excellence
- report the time it took you to solve the assignment
- provide feedback on the assignment, namely how interesting you found it and how challenging it was to answer it.
