
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Levellogger\_Correction

Levellogging dataloggers based on pressure transducers have been shown
to have temperature dependent artifacts, even when manufacturers attemp
to account for temperature effects. Within a logger there are systematic
measurement errors induced by these temperature effects. Systematic
measurement errors exist in both the water and air loggers and must be
addressed separately. Additionally if the loggers are not deployed in
similar thermal regimes then the temperature difference between the two
loggers also needs to be accounted for. Depending upon the strength of
impact of temperature difference on compensation, it may be necessary to
correct even small temperature differences

Within a given logger the temperature compensation models vary for air
and water compenstation.

## Solinst Info

Solinst records pressure at sea level (altitude adjusted) so

WLrecorded = WLlocal \* 0.121\*Elevation Where WLlocal and elevation are
in common units and WLrecorded is pressure at sea level adjusted for
elevation

## External Barometric Data

Use altimeter pressure :

STATION PRESSURE: This is the pressure that is observed at a specific
elevation and is the true barometric pressure of a location. It is the
pressure exerted by the atmosphere at a point as a result of gravity
acting upon the “column” of air that lies directly above the point.
Consequently, higher elevations above sea level experience lower
pressure since there is less atmosphere on which gravity can act. Put
another way, the weight of the atmosphere decreases as one increases in
elevation. Consequently then, in general, for every thousand feet of
elevation gain, the pressure drops about 1 inch of mercury. For example,
locations near 5000 feet (about 1500 meters) above mean sea level
normally have pressures on the order of 24 inches of mercury.

ALTIMETER SETTING: This is the pressure reading most commonly heard in
radio and television broadcasts. It is not the true barometric pressure
at a station. Instead it is the pressure “reduced” to mean sea level
using the temperature profile of the “standard” atmosphere, which is
representative of average conditions over the United States at 40
degrees north latitude. The altimeter setting is the pressure value to
which an aircraft altimeter scale is set so that it will indicate the
altitude (above mean sea level) of the aircraft on the ground at the
location for which the pressure value was determined. The altimeter
setting is an attempt to remove elevation effects from pressure readings
using “standard” conditions.

MEAN SEA LEVEL PRESSURE: This is the pressure reading most commonly used
by meteorologists to track weather systems at the surface. Like
altimeter setting, it is a “reduced” pressure which uses observed
conditions rather than “standard” conditions to remove the effects of
elevation from pressure readings. This reduction estimates the pressure
that would exist at sea level at a point directly below the station
using a temperature profile based on temperatures that actually exist at
the station. In practice the temperature used in the reduction is a mean
temperature for the preceding twelve hours. Mean sea level pressure
should be used with caution at high elevations as temperatures can have
a very profound effect on the reduced pressures, sometimes giving rise
to fictitious pressure patterns and anomalous mean sea level pressure
values.

  - from <https://www.weather.gov/bou/pressure_definitions>

# Compensation Models

Are they unique to each logger or to logger model? Do they vary between
in water and in air? Are temperature difference models specific to
logger pairs, or can they be generalized to a single model?

# Assessing Uncertainty

Actual water level uncertainty consists of sensor uncertainty as well as
any processing uncertainty. 1. Measurement uncertainty of water pressure
2. Measurement uncertainty of barometric pressure 3. Parameter
uncertainty of temperature compensation models (for both models if
unvented loggers are used) 3a. Measurement uncertainty of reference
barometric pressure 3b. Measurement uncertainty of actual water height
4. Paremeter uncertainty of temperature difference models 4a.
Measurement uncertainty of air temperature 4b. Measurement uncertainty
of air temperature

Using vented levelloggers may reduce total uncertainty by removing
multiple sources of pressure and temperature measurement and
compensation uncertainty. I need to explore literature to see if anyone
has done this work previously. I would also need to get some
experimental data from a vented levellogger. Of particular importance
would be knowing if air temperature variation affects the vented
levellogger.

# Water Temperature Impacts

Water level calculated with water\_pressure - air\_pressure
dens\_water\_pressure - dens\_air\_pressure dens\_water\_pressure -
pred\_air\_pressure

They all have linear relationships, but the last one shows the best
alignment between var-sim and stat-sim experiments (fall on the same
linear trends). And not accounting for density is not an option, you can
see the curve from density impacts at higher water and air temperatures.
Once this curve is removed the relationship between air\_pressure and
ex\_air\_pressure requires a linear correction.
