from numpy import *
from mpl_toolkits.mplot3d import axes3d
from pylab import *
import pdb

def x_coord(RA,VRay):
  return multiply(-sin(radians(RA)),VRay)

def y_coord(RA,VRay):
  return multiply(cos(radians(RA)),VRay)

def z_coord(DEC,VRay):
  return multiply(sin(radians(DEC)),VRay)

sorttype = raw_input('Data set (f,*,+): ')

finallist = genfromtxt('finallist.dat',dtype=[('ASTK',str_,16),('RA','f8'),('DEC','f8'),('VRAY','f8')],usecols=(4,9,10,11))

finallist = finallist[where(finallist['ASTK'] == sorttype)]
x = x_coord(finallist['RA'],finallist['VRAY'])
y = y_coord(finallist['RA'],finallist['VRAY'])
z = z_coord(finallist['DEC'],finallist['VRAY'])
for i in range(len(x)-1):
    a = asarray((x[i],y[i],z[i]))
    b,dist = [],[]
    for j in range(len(x)-1):
        if i!=j:
            b =(x[j],y[j],z[j])
            dist.append(norm(a-b))

    dist = sort(asarray(dist))
