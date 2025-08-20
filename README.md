# Assignment A1

Please carefully read the [rules](rules/README.md). After reading this document, copy the [answer sheet](answer-sheet.md) to file `report.md` and provide your answers in the latter. Do not edit this document nor the [answer sheet](answer-sheet.md).


## Objectives

This assignment intends to show you how to represent orbital data under different time and coordinate systems.

## Tasks

### 1 Update your orbit

1. Go to [celestrak](https://celestrak.org/satcat/search.php) and update TLE data of the Earth-orbiting objects you chose in assignment A0.

2. Use your SGP4 propagator to generate an orbit with at least 10 revolutions.

### 2 Implement time system conversions

Convert the time domain of your orbit to the following time standards:

1. UT1: find the time correction tables [here](https://maia.usno.navy.mil/products/daily.htm) and note the following:
    - you only need to retrieve the values of DUT1=UT1-UTC from these tables;
    - you may choose any of the available tables:
        - either those using 1976 Precession/1980 Nutation Theory or IAU2000A   Nutation/Precession Theory, 
        - the tables with predictions for 90 days, 15 days or 1 year.
2. GPS: find the Leap Seconds table [here](https://maia.usno.navy.mil/products/leap-second). 

Represent the different time standards as the number of seconds since 2000/01/01 12:00:00 (also know as *J2000 epoch*). Use the UTC time standard represented as J2000 for your answers to Tasks 3 to 5.

**Opportunity for Excellence**: have your code automatically download the leap-second and time correction tables (unless the local files with this data are already available).

**Opportunity for Excellence**: share your routines that read (and possibly download) the leap-second and time correction tables in a public GitLab repository:

- _only_ share the code that download and reads the data; how this data is used, i.e. the implementation of the actual time standard conversion, _must not_ be shared
- include a `README.md` file that clearly describes:
    - how to use your code
    - a working example showing your code in use
    - the data type returned by your routines
- publish these routines in a public GitLab repository created by you
- announce the availability of the code in the Brightspace forums
- keep track of colleagues that use your code and report it in the answer sheet.

All students are encouraged to use each other's routines described above (and only the code abiding to those limits). There is no penality for doing so.

### 3 Represent your orbit in terms of Keplerian elements 

1. Implement a program that:
    1. reads your orbit data file represented in Cartesian coordinates (henceforth *Cartesian orbit*), produced in Task 1
    2. converts it to Keplerian elements (henceforth *Keplerian orbit*)
    3. write the Keplerian orbit to another file.

2. Represent the Keplerian orbit you produce in 6 plots, with the x-axis showing time and the y-axis showing the values of a different Keplerian element.

**Opportunity for Excellence**: use `stdin` to read the orbit data and `stdout` to write the converted orbit data in your conversion program; use `stdin` to read the converted orbit data in your plotting program; [pipe](https://www.delftstack.com/howto/linux/pipes-in-bash/) the two programs together.

### 4 Assess the errors of the Cartesian to Keplerian and Keplerian to Cartesian coordinates transformation

1. Implement a program that:
    1. reads your Keplerian orbit, produced in Task 3
    2. converts it back to Cartesian coordinates (henceforth *back-converted Cartesian orbit*)
    3. write the back-converted Cartesian orbit to another file.

2. Represent the difference between the original (from Task 1) and back-converted Cartesian orbits in 2 plots, as function of time:
    1. One plot showing the difference in position, along the x, y and z axis
    2. One plot showing the difference in velocity, along the x, y and z axis.

### 5 Represent your orbit in terms of (inertial) Spherical coordinates

1. Implement a program that:
    1. reads your Cartesian orbit, produced in Task 1
    2. converts the position (ignore the velocity) to (inertial) Spherical coordinates (henceforth *Spherical orbit*)
    3. write the Spherical orbit to another file.

2. Represent the Spherical orbit in 3 plots:
    1. One plot showing the two angular coordinates (inertial latitude and intertial longitude) as function of time
    2. One plot showing the altitude (not radius) as function of time assuming:
        1. a spherical Earth
        2. an oblate Earth modelled as an ellipsoid of revolution with the flattening defined by the [WGS84 ellipsoid reference](https://earth-info.nga.mil/index.php?dir=wgs84&action=wgs84)
    3. One plot showing the position of the satellite in the inertial longitude - inertial latitude domain

### 6 Assess the errors of the Cartesian to Spherical and Spherical to Cartesian coordinates transformation

1. Implement a program that:
    1. reads your orbit data file represented in Spherical coordinates (positions only), produced in Task 5
    2. converts it back to Cartesian positions
    3. write the back-converted Cartesian positions to another file.

2. Represent the difference between the original (from Task 1) and back-converted (from this Task) Cartesian positions, as function of time, along the x, y and z axis.

### 7 Represent your orbit in terms of co-rotating Spherical coordinates 

1. Implement a program that:
    1. reads your orbit data file represented in (inertial) Spherical coordinates (positions only), produced in Task 5
    2. converts it to Earth-Centred, Earth-Fixed positions, i.e. by including the mean constant rotation of the Earth in the longitude coordinate (henceforth *Co-rotating orbit*)
    3. write the Co-rotating orbit to another file.

2. Represent the Co-rotating orbit in one plot showing the position of the satellite in the longitude-latitude domain.

**Opportunity for Excellence**: also plot the geographic outline of the continents in the plot of Task 7. You are encouraged to use libraries or data retrieved elsewhere to accomplish this.

3. Determine the time of closest approach of your satellite to the Aerospace Faculty (latitude 51.98988681616558, longitude 4.375752354503242), in the UT1 time standard using the Modified Julian Date representation.

## Final remarks

Don't forget to:

- report Assignment and Code Excellence
- report the time it took you to solve the assignment
- provide feedback on the assignment, namely how interesting you found it and how challenging it was to answer it.
