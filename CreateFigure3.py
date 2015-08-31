"""31.0/.2015 PyOSE: Stacked exomoons with the Orbital Sampling Effect."""

import PyOSE
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
from numpy import pi


# Set stellar parameters
StellarRadius = 696342.  # km
limb1 = 0.3643
limb2 = 0.2807

# Set planet parameters
#PlanetRadius = 63700.  # [km]
PlanetRadius = 63750.
PlanetAxis = 149597870.700  # [km]
PlanetImpact = 0.4  # [0..1.x]; central transit is 0.
PlanetPeriod = 365.25  # [days]

# Set moon parameters
MoonRadius = 17500.  # [km]
MoonAxis = 384000.  # [km]
MoonEccentricity = 0.7  # 0..1
MoonAscendingNode = 30.0  # degrees
MoonLongitudePeriastron  = 50.0  # degrees
MoonInclination = 83.0  # 0..90 in degrees. 0 is the reference plain (no incl).

# Set other parameters
ShowPlanetMoonEclipses = True  # True: the reality; False would be no mutual
#  eclipses. Of course unphysical, but useful for tests and comparisons)
ShowPlanet = False  # True: Planet+Moon; False: Moon only
Noise = 0  # [ppm per minute]; 0 = no noise is added
NumberOfTransits = 0  # How many (randomly chosen) transits are observed;
#  if 0 then all available are sampled (and their number is 10*Quality).
PhaseToHighlight = 0.09  # If no highlighting is desired, choose value < 0
Quality = 250  # Radius of star in pixels --> size of numerical sampling grid
NumberOfSamples = 250  # How many transits are to be sampled


# 3D model
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
MyModelview = PyOSE.modelview(
	StellarRadius, limb1, limb2, PlanetRadius, PlanetImpact,
    MoonRadius, MoonAxis, MoonEccentricity, MoonAscendingNode, 
    MoonLongitudePeriastron, MoonInclination, PhaseToHighlight, Quality)
ax = plt.axes()
ax.arrow(0, 0.4, 0.364, 0, head_width=0.0, head_length=0.0, fc='k', ec='k', 
    zorder = 10)  # direction
#ax.arrow(0.364, 0.4, 0, -0.1, head_width=0.0, head_length=0.0, fc='k', ec='k')
ax.arrow(0.356, 0.0, 0.0, 0.205, head_width=0.0, head_length=0.0, fc='k', ec='k',
    zorder = 10)  # b_S
ax.tick_params(direction='out')

plt.tick_params(axis='both', which='major', labelsize=16)


plt.annotate(r"direction", xy=(0.17, 0.42), size=16)
plt.annotate(r"$b_S$", xy=(0.4, 0.05), size=16)
#plt.annotate(r"$b_P$ ", xy=(-0.12, 0.1), size=16)

plt.xlabel('distance [stellar radii]',fontsize=16)
plt.ylabel('distance [stellar radii]',fontsize=16)
ax.set_aspect('equal')
plt.savefig("figure3-3dview.pdf", bbox_inches='tight')


MyModelview.show()



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
plt.imshow(MyRiverKepler, cmap=cm.gray, interpolation='none', 
	extent=[MyTime[0], -MyTime[0], 1, 0])
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
ax = plt.axes()
ax.arrow(MyTime[0], PhaseToHighlight, 0.001, 0, head_width=0.1, head_length=0.1, 
	fc='k', ec='k')
ax.arrow(-0.5, PhaseToHighlight, 0.001, 0, head_width=0.05, head_length=0.05, 
    fc='k', ec='k')
plt.tick_params(axis='both', which='major', labelsize=16)
plt.xlabel('time around planetary mid-transit [days]',fontsize=16)
plt.ylabel('phase [0..1]',fontsize=16)
plt.axis([-0.5, 0.5, 1, 0], fontsize=16)
ax.tick_params(direction='out')

plt.savefig("figure3-river.pdf", bbox_inches='tight')
plt.show()


# Curve
MyNewCurve = PyOSE.curve(StellarRadius, limb1, limb2, PlanetRadius, PlanetAxis, 
    PlanetImpact, PlanetPeriod, MoonRadius, MoonAxis, MoonEccentricity, 
    MoonAscendingNode, MoonLongitudePeriastron, MoonInclination, 
    ShowPlanetMoonEclipses, ShowPlanet, Quality, NumberOfSamples, Noise, 
    NumberOfTransits)
Time = MyNewCurve[0][1:]
Flux = MyNewCurve[1][1:]
plt.plot(Time, Flux, color = 'k')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.tick_params(axis='both', which='major', labelsize=16)
plt.xlabel('time around planetary mid-transit [days]',fontsize=16)
plt.ylabel('normalized stellar brightness [$10^{-6}$]',fontsize=16)
plt.axis([-0.5, +0.5, -700, 1], set_aspect='equal', fontsize=16)
ax.tick_params(direction='out')

plt.savefig("figure3-curve.pdf", bbox_inches='tight')
plt.show()

# Integral
MyIntegral = PyOSE.integral(StellarRadius, limb1, limb2, PlanetRadius, 
    PlanetAxis, PlanetImpact, PlanetPeriod, MoonRadius, MoonAxis, 
    MoonEccentricity, MoonAscendingNode, MoonLongitudePeriastron, 
    MoonInclination, ShowPlanetMoonEclipses, ShowPlanet, 
    Quality, NumberOfSamples, Noise, NumberOfTransits)
print 'Total occulted stellar flux:', MyIntegral, 'ppm hrs'


#Save curve to Excel file
#book = xlwt.Workbook()
#sheet1 = book.add_sheet('sheet1')
#for i in range(len(Time)):
#    sheet1.write(i,0,Time[i])
#    sheet1.write(i,1,Flux[i])
#book.save("curve.xls")
#book.save(TemporaryFile())


# Debug OSE duration. Use only with noiseless data and eccentricity = 0
v = 2 * pi * PlanetAxis / PlanetPeriod
D = (2 * MoonAxis + 2 * MoonRadius + 2 * StellarRadius) / v
print 'Moon transit duration (analytical):', D, '[days]'

#Get first data point that contains negative flux value
for i in range(len(Time)):
    if Flux[i] < 0:
        FirstDataPoint = Time[i]
        break

#Get last data point that contains negative flux value
for i in range(int(0.5 * len(Time)), len(Time)):
    if Flux[i] == 0:
        LastDataPoint = Time[i-1]
        break

MoonTransitDurationFromData = LastDataPoint - FirstDataPoint
print 'Moon transit duration (from data): ', \
    MoonTransitDurationFromData, '[days]'

ErrorInDuration = ((MoonTransitDurationFromData / D) - 1) * 100
print 'Error:', ErrorInDuration, '[%]. Note: Assumes eccentricity = impact = 0.'
