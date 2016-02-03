"""Create Figure 7a (top panel) of the paper"""

import PyOSE
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
from numpy import pi


# Set stellar parameters
StellarRadius = 0.7 * 696342.  # km
limb1 = 0.5971
limb2 = 0.1172

# Set planet parameters
PlanetRadius = 6371 * 4.8
PlanetAxis = 0.1246 * 149597870.700  # [km]
PlanetImpact = 0.25  # [0..1.x]; central transit is 0.
PlanetPeriod = 16.96862  # [days]

# Set moon parameters
MoonRadius = 6371 * 0.7  # [km]
MoonAxis = 238912.5  # [km] <<<<------
MoonEccentricity = 0.0  # 0..1
MoonAscendingNode = -30.0  # degrees
MoonLongitudePeriastron  = 50.0  # degrees
MoonInclination = 83.0  # 0..90 in degrees. 0 is the reference plain (no incl).

# Set other parameters
ShowPlanetMoonEclipses = True  # True: the reality; False would be no mutual
#  eclipses. Of course unphysical, but useful for tests and comparisons
ShowPlanet = False  # True: Planet+Moon; False: Moon only
Noise = 0  # [ppm per minute]; 0 = no noise is added
NumberOfTransits = 0  # How many (randomly chosen) transits are observed;
#  if 0 then all available are sampled (and their number is 10*Quality).
PhaseToHighlight = 0.2  # If no highlighting is desired, choose value < 0
Quality = 5000  # Radius of star in pixels --> size of numerical sampling grid
NumberOfSamples = 250  # How many transits are to be sampled



# River
MyRiverKepler = PyOSE.river(
    StellarRadius, limb1, limb2,
    PlanetRadius, PlanetAxis, PlanetImpact, PlanetPeriod,
    MoonRadius, MoonAxis, MoonEccentricity, MoonAscendingNode, 
    MoonLongitudePeriastron, MoonInclination,
    ShowPlanetMoonEclipses, Quality, NumberOfSamples, Noise)

# River function returns pixel map. To plot time axis, call function timeaxis
MyTime = PyOSE.timeaxis(
    PlanetPeriod, PlanetAxis, MoonRadius, StellarRadius, Quality)
plt.imshow(MyRiverKepler, cmap=cm.gray, interpolation='none', aspect='auto',
    extent=[MyTime[0], -MyTime[0], 1, 0])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
ax = plt.axes()

T_dur_P = 0.1477 / 2

plt.plot([T_dur_P, T_dur_P], [0, 1], 'k--', linewidth = 1.5)
plt.plot([-T_dur_P, -T_dur_P], [0, 1], 'k--', linewidth = 1.5)
ax.arrow(-0.2, PhaseToHighlight, 0.001, 0, head_width=0.05, head_length=0.01, 
    fc='k', ec='k')
plt.tick_params(axis='both', which='major', labelsize=16)
plt.xlabel('time around planetary mid-transit [days]',fontsize=16)
plt.ylabel('phase',fontsize=16)
plt.axis([-0.2, 0.2, 1, 0], fontsize=16)
ax.tick_params(direction='out')
plt.savefig("fig_7b.pdf", bbox_inches='tight')
#plt.savefig("fig_7b.eps", bbox_inches='tight')
plt.show()