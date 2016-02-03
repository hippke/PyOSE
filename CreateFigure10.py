"""20.07.2015 PyOSE: Stacked exomoons with the Orbital Sampling Effect."""

import PyOSE
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
#import xlwt
#from tempfile import TemporaryFile
from numpy import pi


# Set stellar parameters
StellarRadius = 0.7 * 696342.  # km
limb1 = 0.5971
limb2 = 0.1172

# Set planet parameters
#PlanetRadius = 63700.  # [km]
PlanetRadius = 6371 * 4.8
PlanetAxis = 0.1246 * 149597870.700  # [km]
PlanetImpact = 0.0# 0.25  # [0..1.x]; central transit is 0.
PlanetPeriod = 16.96862  # [days]

# Set moon parameters
MoonRadius = 6371 * 0.7  # [km]
MoonAxis = 238912.5  # [km] <<<<------
MoonEccentricity = 0.0  # 0..1
MoonAscendingNode = 0.0  # degrees
MoonLongitudePeriastron  = 50.0  # degrees
MoonInclination = 83.0  # 0..90 in degrees. 0 is the reference plain (no incl).

# Set other parameters
ShowPlanetMoonEclipses = True  # True: the reality; False would be no mutual
#  eclipses. Of course unphysical, but useful for tests and comparisons)
ShowPlanet = False  # True: Planet+Moon; False: Moon only
Noise = 0  # [ppm per minute]; 0 = no noise is added
NumberOfTransits = 0  # How many (randomly chosen) transits are observed;
#  if 0 then all available are sampled (and their number is 10*Quality).
PhaseToHighlight = 0.2  # If no highlighting is desired, choose value < 0
Quality = 250  # Radius of star in pixels --> size of numerical sampling grid
NumberOfSamples = 250  # How many transits are to be sampled


# Curve
PlanetImpact = 0.0# 0.25  # [0..1.x]; central transit is 0.
MoonInclination = 90.0  # 0..90 in degrees. 0 is the reference plain (no incl).
MyNewCurve = PyOSE.curve(StellarRadius, limb1, limb2, PlanetRadius, PlanetAxis, 
    PlanetImpact, PlanetPeriod, MoonRadius, MoonAxis, MoonEccentricity, 
    MoonAscendingNode, MoonLongitudePeriastron, MoonInclination, 
    ShowPlanetMoonEclipses, ShowPlanet, Quality, NumberOfSamples, Noise, 
    NumberOfTransits)
Time = MyNewCurve[0][1:]
Flux = MyNewCurve[1][1:]
ax = plt.axes()
plt.plot(Time, Flux, color = 'k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.tick_params(axis='both', which='major', labelsize=16)
plt.xlabel('time around planetary mid-transit [days]',fontsize=16)
plt.ylabel('normalized stellar brightness [ppm]',fontsize=16)
plt.axis([-0.2, +0.2, -100, 1], set_aspect='equal', fontsize=16)
ax.tick_params(direction='out')

"""
MoonInclination = 0  # 0..90 in degrees. 0 is the reference plain (no incl).
MyNewCurve = PyOSE.curve(StellarRadius, limb1, limb2, PlanetRadius, PlanetAxis, 
    PlanetImpact, PlanetPeriod, MoonRadius, MoonAxis, MoonEccentricity, 
    MoonAscendingNode, MoonLongitudePeriastron, MoonInclination, 
    ShowPlanetMoonEclipses, ShowPlanet, Quality, NumberOfSamples, Noise, 
    NumberOfTransits)
Time = MyNewCurve[0][1:]
Flux = MyNewCurve[1][1:]
plt.plot(Time, Flux, color = 'k')
"""

"""
PlanetImpact = 0.95# 0.25  # [0..1.x]; central transit is 0.
MoonInclination = 90.0  # 0..90 in degrees. 0 is the reference plain (no incl).
MyNewCurve = PyOSE.curve(StellarRadius, limb1, limb2, PlanetRadius, PlanetAxis, 
    PlanetImpact, PlanetPeriod, MoonRadius, MoonAxis, MoonEccentricity, 
    MoonAscendingNode, MoonLongitudePeriastron, MoonInclination, 
    ShowPlanetMoonEclipses, ShowPlanet, Quality, NumberOfSamples, Noise, 
    NumberOfTransits)
Time = MyNewCurve[0][1:]
Flux = MyNewCurve[1][1:]
plt.plot(Time, Flux, color = 'k', linestyle='dashed')
"""
"""
MoonInclination = 0  # 0..90 in degrees. 0 is the reference plain (no incl).
MyNewCurve = PyOSE.curve(StellarRadius, limb1, limb2, PlanetRadius, PlanetAxis, 
    PlanetImpact, PlanetPeriod, MoonRadius, MoonAxis, MoonEccentricity, 
    MoonAscendingNode, MoonLongitudePeriastron, MoonInclination, 
    ShowPlanetMoonEclipses, ShowPlanet, Quality, NumberOfSamples, Noise, 
    NumberOfTransits)
Time = MyNewCurve[0][1:]
Flux = MyNewCurve[1][1:]
plt.plot(Time, Flux, color = 'k', linestyle='dashed')
"""
plt.annotate(r"$b=0, i_s=90^{\circ}$ ", xy=(0.04, -90), size=16)
plt.annotate(r"$b=0, i_s=0^{\circ}$ ", xy=(-0.04, -60), size=16)
plt.annotate(r"$b=0.95, i_s=0^{\circ}$ ", xy=(-0.045, -44), size=16)
plt.annotate(r"$b=0.95, i_s=90^{\circ}$ ", xy=(-0.05, -5), size=16)


#plt.annotate(r"$b_P$ ", xy=(0.04, 0.02), size=16)
#plt.annotate(r"$b_P$ ", xy=(0.04, 0.02), size=16)
#plt.annotate(r"$b_P$ ", xy=(0.04, 0.02), size=16)


plt.savefig("fig_10.pdf", bbox_inches='tight')
plt.show()
