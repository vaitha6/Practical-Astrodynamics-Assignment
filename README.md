# Assignment A3

Please carefully read the [rules](rules/README.md). After reading this document, copy the [answer sheet](answer-sheet.md) to file `report.md` and provide your answers in the latter. Do not edit this document nor the [answer sheet](answer-sheet.md).

## Objectives

In this assignment, you are tasked to implement disturbing forces in your orbit integrator. As before, you must not use any library that computes forces, orbital disturbances or celestial body locations. In case of doubt, please ask in the Brightspace forum.

In step 2 of Tasks 2 to 5, you are asked to plot radial, longitudinal and latitudinal quantities. This is to be done **in the same plot**, with two different y-axis: one y-axis for the length quantity (radial difference) and another y-axis for the angular quantities (longitudinal and latitudinal differences). There are numerous online tutorials on how to do this with python, e.g.:

- [matplotlib.axes.Axes.twinx](https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.twinx.html)
- [How to Make a Plot with Two Different Y-axis in Python with Matplotlib?](https://cmdlinetips.com/2019/10/how-to-make-a-plot-with-two-different-y-axis-in-python-with-matplotlib/)

Additionally, you must choose a y-axis scale for distance (for radial) and angular (for longitudinal and latitudinal) quantities that allow all lines to be clearly visible in the same plot.

## Tasks

### 1 Update your orbit

1. Go to [celestrak](https://celestrak.org/satcat/search.php) and update TLE data of the Earth-orbiting objects you chose in assignment A0.

2. Generate the sgp4 orbit, similarly to Step 2 of Task 1 of Assignment A2, with at least 3 orbit revolutions.

3. Generate the integrated orbit, similarly to Step 3 of Task 3 of Assignment A2, with at the same time domain as the sgp4 orbit.

### 2 Implement the gravitational disturbance resulting from Earth's oblateness

1. Under the assumption that the only gravitational disturbance affecting the orbit of your satellite is due to Earth's oblateness, as quantified by the $`J_2`$ geopotential coefficient equal to $`1.083\times 10^{-3}`$ (with the gravitational constant equal to $`G=398600.4415 km^3/s^2`$ and $`R = 6378.1363 km`$), compute the resulting orbit considering as initial conditions the state given by the first epoch of the integrated orbit computed in Task 1.

2. Plot the radial, longitudinal and latitudinal difference between the disturbed orbit and the integrated orbit computed in Task 1, as function of time. 

3. Plot the magnitude of gravitational acceleration resulting from Earth's oblateness, as function of time.

**Opportunity for Excellence**: In addition to the effect of $`J_2`$, also include the effect of $`J_{2,2}=1.8155628 10^{-6}`, with $`\Lambda_{2,2}=-14.9287`$ degrees, when addressing steps 2 and 3 of this Task. For the assignment of excellence, plot the radial, longitudinal and latitudinal difference between the orbit perturbed with only the effect of $`J_2`$ and the orbit perturbed by the effect of both $`J_2`$ and $`J_{2,2}`$. Also report the code of that allows you to compute the gravitational attraction caused by $`J_{2,2}`$.


### 3 Implement the effect of drag

1. Considering the following assumptions:

- Your satellite is a perfect sphere,
- Its ballistic coefficient is $`0.01 m^2/kg`$, and 
- Earth's atmosphere is spherical, non-rotating and its density follows the exponential density distribution $`\rho(z)=\rho_0 e^{-z/H}`$, with $`\rho_0=1.2kg/m^3`$ and $`H=8km`$,

compute the orbit under the influence of atmospheric drag, considering as initial conditions the state given by the first epoch of the integrated orbit computed in Task 1.

2. Plot the radial, longitudinal and latitudinal difference between the disturbed orbit and the integrated orbit computed in Task 1, as function of time. 

3. Plot the magnitude of the drag acceleration, as function of time.

**Opportunity for Excellence**: Assume that Earth's atmosphere is rotating with the same angular speed as the Earth, when addressing steps 2 and 3 of this Task. For the assignment of excellence, plot the radial, longitudinal and latitudinal difference between the orbit perturbed with the non-rotating Earth atmosphere and the orbit perturbed by the rotating atmosphere. 

### 4 Implement the effect of radiation pressure force

1. Considering the following assumptions:

- Your satellite is a perfect sphere,
- Its radiation pressure ballistic coefficient is $`0.02 m^2/kg`$, 
- The solar radiation flux is 1361 $`W/m^2`$, and
- The solar illumination of the satellite follows the shadow function proposed by Doornbos (2012)[Thermospheric Density and Wind Determination from Satellite Dynamics](http://repository.tudelft.nl/view/ir/uuid%3A33002be1-1498-4bec-a440-4c90ec149aea/), pp. 58 to 61, assuming total eclipse when the satellite is in both the penumbra and umbra regions.

2. Plot the radial, longitudinal and latitudinal difference between the disturbed orbit and the integrated orbit computed in Task 1, as function of time. 

3. Plot the magnitude of the radiation pressure acceleration, as function of time.

**Opportunity for Excellence**: Implement the geometric shadowing factor $`f_g`$ given by Eq. 3.32 of Doornbos (2012), or Eq. 26 in slide 35 of Lecture 4, thus contemplating the effect of the penumbra, assuming no absorption of solar flux by Earth's atmosphere (the factor $`f_a`$ in Eq. 3.31 or Eq. 27 in the same slide, is 1), when addressing steps 2 and 3 of this Task. For the assignment of excellence, plot the radial, longitudinal and latitudinal difference between the orbit perturbed with the effect of the umbra and the orbit perturbed by both the effect of penumbra and umbra atmosphere. 

### 5 Implement the effect of 3rd body perturbations

1. Considering the following assumptions, defined in an Earth-Centred (pseudo) Inertial reference frame:

- The sun is located at position [-2.141664114391220E+07 -1.336901089766156E+08 -5.795020557937574E+07] km
- The moon is located at position [-3.316594688659666E+05 1.896187710272793E+05 1.235349432352076E+05] km

compute the orbit under the influence of the gravitational attraction of the Sun and the Moon, considering as initial conditions the state given by the first epoch of the integrated orbit computed in Task 1.

2. Plot the radial, longitudinal and latitudinal difference between the disturbed orbit and the integrated orbit computed in step 2 of Task 1, as function of time. 

3. Plot the magnitude of the 3rd body accelerations caused by the Sun and Moon separately, as function of time.

**Opportunity for Excellence**: Retrieve the actual positions of the Sun and Moon provided by JPL's [Horizons System](https://ssd.jpl.nasa.gov/horizons/) in an automated way and replace the fixed Sun and Moon coordinates given above with the ones provided by JPL, when addressing steps 2 and 3 of this Task. For the assignment of excellence, report the code of that allows you to retrieve the Sun and Moons ephemerides.

## Final remarks

Don't forget to:

- report Assignment and Code Excellence
- report the time it took you to solve the assignment
- provide feedback on the assignment, namely how interesting you found it and how challenging it was to answer it.
