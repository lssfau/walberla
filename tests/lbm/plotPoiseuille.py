#!/usr/bin/python

import numpy as np
import matplotlib
matplotlib.use('QT4Agg')
import matplotlib.pyplot as plt
import sys

import py_waLBerla.Config



matplotlib.rcParams.update({'font.size': 12})




def analyticalSolution( config ):
    c = config
    
    # Parameters from Config file
    diameter         = float( c['diameter'] )
    aspectRatio      = float( c['aspectRatio'] ) 
    cellsPerDiameter = float( c['cellsPerDiameter'] )
    viscosity        = float( c['viscosity'] )
    omega            = float( c['omega'] )
    acceleration_l   = float( c['acceleration_l'] )
    density          = float( c['density'] )
    
    # Derived quantities
    height = diameter * aspectRatio
    dx = diameter / cellsPerDiameter
    dt = ( 2.0 - omega ) * dx * dx / ( 6 * omega * viscosity );
    dm = density * dx * dx * dx
    acceleration = acceleration_l * dx / ( dt * dt )
    radius = diameter / 2.0 
    viscosity_l = viscosity * dt / ( dx * dx )

    cellsPerRadius = cellsPerDiameter / 2
    x  = np.arange( - cellsPerDiameter/2 , cellsPerDiameter/2  ) 
    x += 0.5
    
    if c['geometry'] == "pipe":
        geometryFactor = 1.0
    else:
        geometryFactor = 2.0
    
    print "geomFactor " + str( geometryFactor )
    
    vel = geometryFactor * acceleration_l / ( 4 * viscosity_l ) * ( cellsPerRadius * cellsPerRadius - x * x )
    return ( x + cellsPerRadius - 0.5, vel )
    
    

nr = sys.argv[1]

configFile = open('poiseuille.prm').read()
config = py_waLBerla.Config.parse( configFile )["Poiseuille"]   

zVelX, zVelY    = np.loadtxt( 'velocity_' + str(nr) + ".dat", unpack=True )
rhoX, rhoY      = np.loadtxt( 'density_'  + str(nr) + ".dat", unpack=True )

fig = plt.gcf()

ax = fig.add_subplot( 411 )
ax.plot( rhoX, rhoY )
plt.title("Density")

ax = fig.add_subplot( 412 )
ax.plot( zVelX, zVelY  )
plt.title("Z Velocity")


r, analytical = analyticalSolution(config)
ax = fig.add_subplot( 413 )
ax.plot( r, analytical )
plt.title("Analytical")


ax = fig.add_subplot( 414 )
print (  zVelY  )
print ( analytical )
diff = zVelY - analytical
print ( diff  )
ax.plot( zVelX, diff )

#sim, = ax.plot( zVelX, zVelY  )
#ana, = ax.plot( r, analytical )
#plt.legend( [sim, ana], ['Sim', 'Ana'], loc=8 )

plt.title("Analytical vs Simulation")


diff = analytical - zVelY
error = np.sum (diff * diff) / len( diff )
print error

plt.show()




