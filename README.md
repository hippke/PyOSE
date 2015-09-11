## Stacked exomoons with the Orbital Sampling Effect 
*We are currently preparing a paper to introduce this new code. Until then, please refer to Heller 2014, ApJ 787*

How to get started:

1. Define parameters for the star, planet and moon [all as floats]:
  * Star: Stellar radius, limb-darkening parameters
  * Planet: Radius, semimajor axis, impact parameter, period
  * Moon: Radius, semimajor axis, eccentricity, ascending node, longitude at periastron, inclination

2. Define simulation parameters:
  *  Account for mutual planet/moon eclipses [bool]
  *  Show the planet in the plots [bool]
  *  Add noise to the simulation [ppm/minute, can be zero]
  *  How many transits are to be sampled [integer]
  *  How many (randomly chosen) transits are observed [integer, 0=all]
  *  Moon phase to highlight in plots [0..1]
  *  Quality setting [integer, defines pixel radius of star; scales other bodies]
   
#### Generate a 3D-view of these parameters for visual verification:
 
![ScreenShot](http://www.jaekle.info/osescreenshots/git1.png)

#### Generate a riverplot covering phase 0..1:

![ScreenShot](http://www.jaekle.info/osescreenshots/git2.png)

#### Generate the stacked lightcurve:

![ScreenShot](http://www.jaekle.info/osescreenshots/git3.png)

#### Calculate the total occulted stellar flux: 
```
-7916.17381004 ppm hrs
```

You can save all figures as PDF, PNG, EPS, etc., and the time/flux data as CSV, Excel etc.
