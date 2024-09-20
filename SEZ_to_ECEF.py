## Script Name: SEZ_to_ECEF.py

## Usage: python3 SEZ_to_ECEF.py o_lat_deg o_long_deg o_hae_km s_km e_km z_km

## Parameters:
# o_lat_deg: Original SEZ Latitude, in degrees
# o_long_deg: Original SEZ Longitude, in degrees
# o_hae_km: Original SEZ Height Above Ellipsoid, in kilometers
# s_km: S-magnitude, in kilometers
# e_km: E-magnitude, in kilometers
# z_km: Z-magnitude, in kilometers

## Output: Converts given SEZ-vector into the equivalent ECEF-vector

## Written by Carl Hayden

## Importing Libraries
import math # Importing Mathematics Library
import sys # Importing Argument-Reading Library
import numpy # Importing numpy Library (matrix math!)

## Defining Constants
R_Earth = 6378.1363 # Radius of Earth in km
e_Earth = 0.081819221456 # Eccentricity of Earth

## Defining Other Dependent Functions
def calc_denom(e_Earth, o_lat_rad):
    return math.sqrt(1.0 - e_Earth ** 2.0 * math.sin(o_lat_rad) ** 2.0)

## Initialize Script Arguments
o_lat_deg = float('nan') # Latitude in degrees 
o_long_deg = float('nan') # Longitude in degrees
o_hae_km = float('nan') # Height above ellipsoid in kilometers
s_km = float('nan') # S-magnitude, in kilometers
e_km = float('nan') # E-magnitude, in kilometers
z_km = float('nan') # Z-magnitude, in kilometers

## Parse Script Arguments
if len(sys.argv)==7:
    o_lat_deg = float(sys.argv[1])
    o_long_deg = float(sys.argv[2])
    o_hae_km = float(sys.argv[3])
    s_km = float(sys.argv[4])
    e_km = float(sys.argv[5])
    z_km = float(sys.argv[6])

else:
    print(\
        'Usage: '\
        'python3 SEZ_to_ECEF.py o_lat_deg o_long_deg o_hae_km s_km e_km z_km'\
        )
    exit()

## Main Script
o_lat_rad = o_lat_deg * math.pi / 180
o_long_rad = o_long_deg * math.pi / 180
denom = calc_denom(e_Earth, o_lat_rad)
C_E = R_Earth/denom
S_E = (R_Earth * (1-e_Earth**2)) / denom
r_x_km = (C_E + o_hae_km) * math.cos(o_lat_rad) * math.cos(o_long_rad)
r_y_km = (C_E + o_hae_km) * math.cos(o_lat_rad) * math.sin(o_long_rad)
r_z_km = (S_E + o_hae_km) * math.sin(o_lat_rad)

SEZ_Vect = numpy.array([[s_km], [e_km], [z_km]])
Rotation1 = numpy.array([[math.sin(o_lat_rad), 0, math.cos(o_lat_rad)], 
             [0, 1, 0], 
             [-math.cos(o_lat_rad), 0, math.sin(o_lat_rad)]])
Rotation2 = numpy.array([[math.cos(o_long_rad), -math.sin(o_long_rad), 0],
             [math.sin(o_long_rad), math.cos(o_long_rad), 0],
             [0, 0, 1]])

Calc1 = numpy.dot(Rotation1, SEZ_Vect)
Calc2 = numpy.dot(Rotation2, Calc1)

x_ECEF_km = Calc2[0] + r_x_km
y_ECEF_km = Calc2[1] + r_y_km
z_ECEF_km = Calc2[2] + r_z_km

print('\nSEZ Vector: [' + str(Calc2[0])+ ', ' + str(Calc2[1]) + ', ' + str(Calc2[2]) + ']')
print('X-Val ECEF: ' + str(x_ECEF_km) + ' km')
print('Y-Val ECEF: ' + str(y_ECEF_km) + ' km')
print('Z-Val ECEF: ' + str(z_ECEF_km) + ' km\n')