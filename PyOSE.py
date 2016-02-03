"""20.07.2015 PyOSE: Stacked exomoons with the Orbital Sampling Effect."""

# Stacked exomoons with the Orbital Sampling Effect (OSE, Heller 2014, ApJ 787)
#
#             _.--"~~ __"-.     Copyright 2015 Michael Hippke and contributors
#          ,-"     .-~  ~"-\    OSE-Sampler is free software
#        .^       /       ( )   made available under the MIT License
#       {_.---._ /         ~        http://opensource.org/licenses/MIT
#       /    .  Y
#      /      \_j
#     Y     ( --l__
#     |            "-.
#     |      (___     \         If you make use of OSE-Sampler in your work,
#     |        .)~-.__/         please cite our paper:
#     l        _)                   Hippke & Heller 2015, arXiv, ADS, BibTeX
#      \      "l
#       \       \
#        \       ^.
#         ^.       "-.
#           "-._      ~-.___,
#               "--.._____.^
#
# Contributors: Add your name here if you like
# Michael Hippke - initiator
# Rene Heller - requirements and OSE Theory (Heller 2014, ApJ 787)
# Stefan Czesla - conversion of some functions from Pascal to Python


import matplotlib.pylab as plt
import numpy as np
from numpy import pi, sqrt, sin, cos, arcsin, arccos, radians, power, degrees
from PyAstronomy import pyasl
from PyAstronomy.modelSuite import forTrans
from matplotlib.pyplot import cm

def modelview(
    StellarRadius, limb1, limb2, PlanetRadius, PlanetImpact, MoonRadius, 
    MoonAxis, MoonEccentricity, MoonAscendingNode,  MoonLongitudePeriastron, 
    MoonInclination, PhaseToHighlight, Quality):
    """Calculate the 3D model view. This is the vector version (v2)"""

    # Convert values from km to internal measure (stellar radius = 1)
    PlanetRadius = PlanetRadius / StellarRadius
    MoonRadius = MoonRadius / StellarRadius
    MoonAxis = MoonAxis / StellarRadius

    # Make background unicolor black or white
    #allwhite = plt.Circle((0, 0), 100, color=(1, 1, 1))
    #allblack = plt.Circle((0, 0), 100, color=(0, 0, 0))    
    #plt.gcf().gca().add_artist(allwhite)

    # Alternatively plot a gradient map as background
    X = [[-1, 0], [0, 1]]
    plt.imshow(X, interpolation='bicubic', cmap=cm.gray,
          extent=(-1.1, 1.1, -1.1, 1.1), alpha=1)


    # Star
    StarQuality = Quality
    StarQuality = 100

    for i in range(StarQuality):
        Impact = i / float(StarQuality)
        LimbDarkening = QuadraticLimbDarkening(Impact, limb1, limb2)
        Sun = plt.Circle((0, 0), 1 - i / float(StarQuality), 
            color=(LimbDarkening, LimbDarkening, LimbDarkening))
        # for a yellow shaded star, replace the last LimbDarkening with 0
        plt.gcf().gca().add_artist(Sun)

    # Moon's orbit: Kepler ellipse normalized to 1x1 stellar radii
    Ellipse = pyasl.KeplerEllipse(
        MoonAxis, MoonAxis, 
        e = MoonEccentricity, 
        Omega = MoonAscendingNode, 
        w = MoonLongitudePeriastron, 
        i = MoonInclination)
    NumerOfSamples = 1 / float(Quality)
    time = np.linspace(0., MoonAxis, Quality)
    coordinates = np.zeros((len(time), 3), dtype=np.float) 
    for i in xrange(len(time)):
        coordinates[i,::] = Ellipse.xyzPos(time[i]) 

    OccultedDataPoints = []  # To clip the moon ellipse "behind" the planet
    for i in range(Quality / 2):
        CurrentMoonImpact = coordinates[i, 1]
        CurrentHorizontalMoonPosition = coordinates[i, 0]
        CurrentDistanceMoonPlanet = sqrt(CurrentMoonImpact ** 2
            + CurrentHorizontalMoonPosition ** 2) * StellarRadius
        if CurrentDistanceMoonPlanet < PlanetRadius * StellarRadius:
            OccultedDataPoints.append(i)
    if len(OccultedDataPoints) > 0:
        FirstOcculted = OccultedDataPoints[0]
        LastOcculted  = OccultedDataPoints[-1] + 1
        plt.plot(-coordinates[:FirstOcculted,0], 
                  coordinates[:FirstOcculted,1] + PlanetImpact, 'k', zorder = 5)
        plt.plot(-coordinates[LastOcculted:,0], 
                  coordinates[LastOcculted:,1]  + PlanetImpact, 'k', zorder = 5)
    else:
        plt.plot(-coordinates[::,0], 
                  coordinates[::,1] + PlanetImpact, 'k', zorder = 5)

    # Planet
    PlanetCircle = plt.Circle((0, PlanetImpact), PlanetRadius, color = 'k', 
        zorder = 4)
    plt.gcf().gca().add_artist(PlanetCircle)

    # Moon
    PosPhase = PhaseToHighlight / float(NumerOfSamples)
    coordinates[1,::] = Ellipse.xyzPos(time[PosPhase]) 
    if PhaseToHighlight < 0.5:
        CurrentOrder = 2  # behind the planet
    else:
        CurrentOrder = 4  # in front of the planet
    MoonCircle = plt.Circle((-coordinates[1,0], 
        coordinates[1,1] + PlanetImpact), MoonRadius, color = 'k', 
        zorder = CurrentOrder)
    plt.gcf().gca().add_artist(MoonCircle)
    # Square not functional?
    plt.axis([-1.1, +1.1, -1.1, +1.1], set_aspect='equal', fontsize=16)  
    return plt


def riverwithoutnoise(
    StellarRadius, limb1, limb2, PlanetRadius, PlanetAxis, PlanetImpact, 
    PlanetPeriod, MoonRadius, MoonAxis, MoonEccentricity, MoonAscendingNode, 
    MoonLongitudePeriastron, MoonInclination, ShowPlanetMoonEclipses, Quality, 
    NumberOfSamples):
    """Core function projecting Kepler Moon Ellipse onto PixelGrid"""

    # Moon's orbit: Kepler ellipse normalized to 1x1 stellar radii
    NormalizedMoonAxis = MoonAxis / StellarRadius
    NormalizedMoonRadius = MoonRadius / StellarRadius
    CounterFullEclipses = 0
    CounterPartialEclipses = 0
    CurrentEclipsedRatio = 0
    Ellipse = pyasl.KeplerEllipse(
        NormalizedMoonAxis, NormalizedMoonAxis, 
        e = MoonEccentricity, 
        Omega = MoonAscendingNode, 
        w = MoonLongitudePeriastron, 
        i = MoonInclination)
    time = np.linspace(0., NormalizedMoonAxis, NumberOfSamples)
    coordinates = np.zeros((len(time), 3), dtype=np.float) 
    for i in xrange(len(time)):
        coordinates[i,::] = Ellipse.xyzPos(time[i]) 
    StretchSpace = 5  # 5 Grid size [multiples of planet transit duration]
    UnStretchedTransitDuration = PlanetPeriod / pi * arcsin(sqrt(
        (MoonRadius + StellarRadius) ** 2) / PlanetAxis)
    cache = np.zeros((NumberOfSamples, Quality), dtype=np.float)
    output = np.zeros((NumberOfSamples, StretchSpace * Quality), dtype=np.float)
   
    ma = forTrans.MandelAgolLC()  # Prepare light curves
    ma["T0"] = 0
    ma["b"] = 0. 
    ma["linLimb"] = limb1
    ma["quadLimb"] = limb2
    ma["per"] = PlanetPeriod
    ma["a"] = PlanetAxis / StellarRadius
    ma["p"] = MoonRadius / StellarRadius
    time = np.linspace(ma["per"] - (0.5 * UnStretchedTransitDuration), 
                       ma["per"] + (0.5 * UnStretchedTransitDuration), Quality)
    
    for i in range(NumberOfSamples):  # Fetch curves: The core of this function
        CurrentMoonImpact = coordinates[i, 1]  # Position-dependent moon impact
        CurrentHorizontalMoonPosition = coordinates[i, 0]
        ma["i"] = ImpactToInclination(
            CurrentMoonImpact + PlanetImpact, StellarRadius, PlanetAxis)
        cache[i,::] = ma.evaluate(time)  # Fetch each curve
        
        # Mutual eclipses: calculate distance moon --> planet [km]
        if ShowPlanetMoonEclipses:
            CurrentDistanceMoonPlanet = sqrt(CurrentMoonImpact ** 2
                + CurrentHorizontalMoonPosition ** 2) * StellarRadius
            CurrentEclipsedRatio = EclipsedRatio(
                CurrentDistanceMoonPlanet, PlanetRadius, MoonRadius)
        
        for k in range(Quality): # Transform flux from e.g. 0.995 to -5ppm
            cache[i, k] = -(1 - cache[i, k]) * 10 ** 6
            # And reduce the flux according to the eclipsed area
            if ShowPlanetMoonEclipses and cache[i, k] < 0:
                cache[i, k] = -(1 - cache[i, k]) * (1 - CurrentEclipsedRatio)

    for i in range(NumberOfSamples):  # Apply time shift due to moon orbit
        for k in range(Quality):
            HorizontalPixelPosition = (coordinates[i, 0] * Quality) / 2
            MidShift = (0.5 * StretchSpace - 0.5) * Quality
            output[i, k + HorizontalPixelPosition + MidShift] = cache[i, k]

    return output


def river(
    StellarRadius, limb1, limb2, PlanetRadius, PlanetAxis, PlanetImpact, 
    PlanetPeriod, MoonRadius, MoonAxis, MoonEccentricity, MoonAscendingNode, 
    MoonLongitudePeriastron, MoonInclination, ShowPlanetMoonEclipses, 
    Quality, NumberOfSamples, Noise):
    """Draw the river plot. Get map from RiverWithoutNoise."""

    StretchSpace = 5  # 5 times longer than a transit duration   
    map = riverwithoutnoise(  # Get noiseless map
        StellarRadius, limb1, limb2,
        PlanetRadius, PlanetAxis, PlanetImpact, PlanetPeriod,
        MoonRadius, MoonAxis, MoonEccentricity, MoonAscendingNode, 
        MoonLongitudePeriastron, MoonInclination,
        ShowPlanetMoonEclipses, Quality, NumberOfSamples)

    # Calculate helper variables for noise level per pixel
    PlanetTransitDuration = PlanetPeriod / pi * arcsin(sqrt(
        (2 * StellarRadius + 2 * PlanetRadius) ** 2) / PlanetAxis)
    Pixeltimestep = PlanetTransitDuration / float(Quality)
    timeResolution = Pixeltimestep * 24 * 60

    # Prepare the amount of noise per pixel (=integration)
    # NoisePerIntegration [ppm/integration] = Noise [ppm/minute]
    #     / sqrt(Integrationtime / 1 [Minute])
    NoisePerIntegration = Noise / sqrt(timeResolution)

    for i in range(NumberOfSamples): # Now spread the noise
        for k in range(StretchSpace * Quality):  # Check if noise is desired.           
            if Noise <= 0:  # Skip the (computationally expensive) noise calc
                RandomNumber = 0
            else:
                # Gaussian noise of mean=0; stdev=NoisePerIntegration
                # Can be swapped for any other noise law, or even real data
                RandomNumber = np.random.normal(0, NoisePerIntegration, 1)            
            map[i, k] = map[i, k] + RandomNumber  # Add the noise to the map

    return map


def timeaxis(PlanetPeriod, PlanetAxis, MoonRadius, StellarRadius, Quality):
    Duration = PlanetPeriod / pi * arcsin(sqrt(
        (StellarRadius + MoonRadius) ** 2) / PlanetAxis)
    StretchSpace = 5 * Quality
    time = np.linspace(-2.5 * Duration, 2.5 * Duration, StretchSpace)
    return time


def curve(
    StellarRadius, limb1, limb2, PlanetRadius, PlanetAxis, PlanetImpact, 
    PlanetPeriod, MoonRadius, MoonAxis, MoonEccentricity, MoonAscendingNode, 
    MoonLongitudePeriastron, MoonInclination, PlanetMoonEclipses, ShowPlanet, 
    Quality, NumberOfSamples, Noise, NumberOfTransits):
    """Calculate the curve. Use map as source."""
    HorizontalResolution = 5 * Quality
    curve = np.zeros((2, HorizontalResolution), dtype=np.float)
    curve[0,::] = timeaxis(
        PlanetPeriod, PlanetAxis, MoonRadius, StellarRadius, Quality)

    # Get the MapMoon with noise; if noise=0 then none will be added
    MapMoon = river(
        StellarRadius, limb1, limb2,
        PlanetRadius, PlanetAxis, PlanetImpact, PlanetPeriod,
        MoonRadius, MoonAxis, MoonEccentricity, MoonAscendingNode, 
        MoonLongitudePeriastron, MoonInclination, PlanetMoonEclipses, Quality, 
        NumberOfSamples, Noise)

    if ShowPlanet:  # If wanted: Add noiseless planet curve
        MapPlanet = river(StellarRadius, limb1, limb2,
                    PlanetRadius, PlanetAxis, PlanetImpact, PlanetPeriod,
                    PlanetRadius, 1., 0., 0., 0., 90., False, Quality, 1, 0)

    if NumberOfTransits == 0:  # Full OSE curve is requested
        for i in range(HorizontalResolution):
            Collector = 0
            for k in range(NumberOfSamples):
                Collector = Collector + MapMoon[k,i]
            Chartflux = Collector / (NumberOfSamples)
            curve[1, i] = Chartflux
    else:  # Get a subset
        MapMoonObserved = np.zeros((NumberOfTransits, HorizontalResolution), 
            dtype=np.float)

        # Fetch #NumberOfTransits transits and push them to #MapMoonObserved
        for k in range(NumberOfTransits):
            # Instead of the random selection, one might introduce a rule to
            # select orbits according to phase time, e.g. observing the transit
            # at phase=0.17 at first visit, then 0.34, then 0.51...
            MyRandomTransit = np.random.uniform(0, NumberOfSamples)
            # <== apply other rule here
            for i in range(HorizontalResolution - 1):
                MapMoonObserved[k, i] = MapMoonObserved[k, i] + \
                                        MapMoon[MyRandomTransit, i]

    # Collector routine now using subset:MapMoonObserved; compensate for flux
    for i in range(HorizontalResolution):
        Collector = 0
        for k in range(NumberOfTransits):
            Collector = Collector + MapMoonObserved[k, i]
            Chartflux = Collector / float(NumberOfSamples)
            curve[1, i] = Chartflux * (NumberOfSamples / float(NumberOfTransits))

    if ShowPlanet:
        curve[1, ::] = curve[1, ::] + MapPlanet[0, ::]

    return curve


def integral(
    StellarRadius, limb1, limb2,
    PlanetRadius, PlanetAxis, PlanetImpact, PlanetPeriod,
    MoonRadius, MoonAxis, MoonEccentricity, MoonAscendingNode, 
    MoonLongitudePeriastron, MoonInclination, PlanetMoonEclipses, ShowPlanet, 
    Quality, NumberOfSamples, Noise, NumberOfTransits):
    """Calculate the integral. Use curve as source."""
    MyCurve = curve(
        StellarRadius, limb1, limb2, PlanetRadius, PlanetAxis, PlanetImpact, 
        PlanetPeriod, MoonRadius, MoonAxis, MoonEccentricity, MoonAscendingNode,
        MoonLongitudePeriastron, MoonInclination, PlanetMoonEclipses, 
        ShowPlanet, Quality, NumberOfSamples, Noise, NumberOfTransits)

    time = MyCurve[0][::]
    flux = MyCurve[1][::]
    TimeLength = time[0]
    TimeStep = time[1] - time[0]
    NumberOfBins = TimeLength / TimeStep
    TotalFlux = np.sum(flux)

    return TotalFlux * TimeLength / NumberOfBins * 24  # [ppm hrs]


def CircleCircleIntersect(radius1, radius2, distance):
    """Calculates area of asymmetric "lens" in which two circles intersect
    Source: http://mathworld.wolfram.com/Circle-CircleIntersection.html"""
    return radius1 ** 2 * (arccos(((distance ** 2) + (radius1 ** 2) -
        (radius2 ** 2)) / (2 * distance * radius1))) + ((radius2 ** 2) *
        (arccos((((distance ** 2) + (radius2 ** 2) - (radius1 ** 2)) /
        (2 * distance * radius2))))) - (0.5 * sqrt((-distance + radius1 + 
        radius2) * (distance + radius1 - radius2) * (distance - radius1 + 
        radius2) * (distance + radius1 + radius2)))


def EclipsedRatio(CurrentDistanceMoonPlanet, PlanetRadius, MoonRadius):
    """Returns eclipsed ratio [0..1] using CircleCircleIntersect"""
    HasFullEclipse = False
    HasMutualEclipse = False
    CurrentEclipsedRatio = 0
    if abs(CurrentDistanceMoonPlanet) < (PlanetRadius + MoonRadius):
        HasMutualEclipse = True
        if ((PlanetRadius - MoonRadius) > abs(CurrentDistanceMoonPlanet)):
            HasFullEclipse = True
    # For partial eclipses, get the fraction of moon eclipse using transit
    if HasMutualEclipse and not HasFullEclipse:
        CurrentEclipsedRatio = CircleCircleIntersect(
            PlanetRadius, MoonRadius, CurrentDistanceMoonPlanet)
        # ...and transform this value into how much AREA is really eclipsed
        CurrentEclipsedRatio = CurrentEclipsedRatio / (pi * MoonRadius ** 2)
    if HasFullEclipse:
        CurrentEclipsedRatio = 1.0
    return CurrentEclipsedRatio


def ImpactToInclination(PlanetImpact, StellarRadius, PlanetAxis):
    """Converts planet impact parameter b = [0..1.x] to inclination [deg]"""
    return degrees(arccos(PlanetImpact * StellarRadius / PlanetAxis))


def QuadraticLimbDarkening(Impact, limb1, limb2): 
    """Quadratic limb darkening. Kopal 1950, Harvard Col. Obs. Circ., 454, 1"""
    return 1 - limb1 * (1 - Impact) - limb2 * (1 - Impact) ** 2
