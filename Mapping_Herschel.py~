#!/astro/apps/pkg/python64/bin/python

####################################
#Ramon Sharma                      #
#University of Washington          #
#Herschel 3-D Galactic Map Project #
#11/25/2013                        #
####################################

from pylab import *
from numpy import *
import pdb
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import axes3d
import os, sys, inspect

home = '/astro/users/lambdacd/Herschel/'

#### Herschel's famous slice uses DEC = 0, STD = 18
#### To use J1690, replace all relevant manipulations of NEWRA/NEWDEC with OLDRA/OLDDEC

### The sizes of the sets are as follows:
# Published Data (SET1): 551
# Asterisk Data (SET2): 130
# Unpublished Data (SET3): 514

FOV = 15 #Arcminutes
Count = 0
X,Y,Z,DECArray,RAArray,XG,YG,ZG = [],[],[],[],[],[],[],[]
XS,YS,ZS = [],[],[]
VRAY1AR,FIELDS,STARS,GLAT,GLONG = [],[],[],[],[]
PROX,PROY,PROZ = [],[],[]
Xp,Yp,Zp,DECArrayp,RAArrayp,XGp,YGp,ZGp = [],[],[],[],[],[],[],[]
XSp, YSp, ZSp = [],[],[]
VRAY1ARp, FIELDSp, STARSp, GLATp, GLONGp = [],[],[],[],[]
VRAY1ARp1, VRAY1ARp2 = [], []
GLATp1, GLATp2, GLONGp1, GLONGp2 = [], [], [], []
star_mag, gage_mag = [], []
span = .4
slice = 50
SPHRAD = 10.

#### Galactic Coordinate Axes in J2000 RA&DEC
REFRA=[266.417,192.859,318.004] #RA of SAG A*, l=0 b=90, l=90 b=0
REFDEC=[-29.008,27.128,48.330] #DEC of SAG A*, l=0 b=90, l=90 b=0

#### Galactic Coordinate Axes in J1690 RA&DEC
#REFRA=[263.004367,190.231959,316.153294] #RA of SAG A*, l=0 b=90, l=90 b=0
#REFDEC=[-28.897618,28.301223,47.453070] #DEC of SAG A*, l=0 b=90, l=90 b=0

REFDIST=[500,500,500] #Arbitrary (large) visual rays

SIRDIST=2.54 #Distance to Sirius in parsecs
#try:
#  filename = sys.argv[1]
#except:
#  print 'Invalid input file'

def ToRad(Deg):
  Rad = []
  for i in range(len(Deg)):
    if Deg[i]>0:
      Rad.append(radians(Deg[i]) % 2*pi)
    else:
      Rad.append(radians(360+Deg[i]) % 2*pi)
      #Rad.append(radians(Deg[i]))
  return Rad
  
def XCoord(RA,VRay):
  return -sin(radians(RA))*VRay

def YCoord(RA,VRay):
  return cos(radians(RA))*VRay

def ZCoord(DEC,VRay):
  return sin(radians(DEC))*VRay

def XSCoord(RA, DEC, VRay):
  return VRay*sin(DEC)*cos(RA)

def YSCoord(RA, DEC, VRay):
  return VRay*sin(DEC)*sin(RA)

def ZSCoord(RA, DEC, VRay):
  return VRay*cos(DEC)



def Project(GLAT,GLONG):
  theta = pi / 2.
  x,y = [],[]
  GLAT = radians(GLAT)
  GLONG = radians(GLONG)
  for i in range(len(GLAT)):
    if GLAT[i] > pi/2.:
      GLAT[i] = pi - GLAT[i]
    if GLAT[i] < -pi/2.:
      GLAT[i] = -pi - GLAT[i]
    if GLONG[i] < -pi:
      GLONG[i] = GLONG[i] + 2.*pi
    elif GLONG[i] > pi:
      GLONG[i] = GLONG[i] - 2.*pi
      
    if GLAT[i] == pi / 2.:
      theta = pi/2.
    elif GLAT[i] == -pi / 2.:
      theta = -pi/2.
    else:
      t1 = GLAT[i]
      while 1:
        t2 = t1 - (2.*t1+sin(2.*t1)-pi*sin(GLAT[i]))/(2.+2.*cos(2.*t1))
        if fabs(t2-t1) < float(1.0e-10):
          break
        else:
          t1 = t2
      theta = t2
    x.append((-2.*sqrt(2.)/pi)*GLONG[i]*cos(theta))
    y.append(sqrt(2.)*sin(theta))
  return x, y

def Hipparcos(whra,whdec,vray1,hipra,hipdec,hipdist):
  ra,dec,dist = [],[],[]
  whra = float64(asarray(whra))
  whdec = float64(asarray(whdec))
  vray1 = asarray(vray1)

  for i in range(len(hipra)):
    for j in range(len(whra)):
      if (fabs(hipra[i] - whra[j]) <= span ) and (fabs(hipdec[i] - whdec[j]) <= span):
        ra.append(hipra[i])
        dec.append(hipdec[i])
        dist.append(hipdist[i])
  return asarray(ra),asarray(dec),asarray(dist)/SIRDIST


def ang_sep(ra1, dec1, ra2, dec2): #Outputs degree separation of a plate object from the boresight

  DEG_PER_HR = 360. / 24.             # degrees per hour
  DEG_PER_MIN = DEG_PER_HR / 60.      # degrees per min
  DEG_PER_S = DEG_PER_MIN / 60.       # degrees per sec
  DEG_PER_AMIN = 1./60.               # degrees per arcmin
  DEG_PER_ASEC = DEG_PER_AMIN / 60.   # degrees per arcsec
  RAD_PER_DEG = pi / 180.
  
  ra1 = asarray(ra1);  ra2 = asarray(ra2)
  dec1 = asarray(dec1);  dec2 = asarray(dec2)

  ra1 = ra1 * RAD_PER_DEG           # convert to radians
  ra2 = ra2 * RAD_PER_DEG
  dec1 = dec1 * RAD_PER_DEG
  dec2 = dec2 * RAD_PER_DEG

  sra1 = sin(ra1);  sra2 = sin(ra2)
  cra1 = cos(ra1);  cra2 = cos(ra2)
  sdec1 = sin(dec1);  sdec2 = sin(dec2)
  cdec1 = cos(dec1);  cdec2 = cos(dec2)

  csep = cdec1*cdec2*(cra1*cra2 + sra1*sra2) + sdec1*sdec2

  # An ugly work-around for floating point issues.
  #if np.any(csep > 1):  print csep
  csep = where(csep > 1., 1., csep)

  degsep = arccos(csep) / RAD_PER_DEG
  # only works for separations > 0.1 of an arcsec or  >~2.7e-5 dec
  degsep = where(degsep < 1e-5, 0, degsep)
  return degsep




finallist = open('finallist.dat', 'r')

DataSet = raw_input('Data set (f,*,+): ')
DecSlice = raw_input('Input a slicing declination (DEG): ')
DecRange = raw_input('Input a slicing range (DEG): ')
HipCall = raw_input('Hipparcos data? ')
if HipCall != '':
     SuperImp = raw_input('Superimpose data when plotting? ')
else:
     SuperImp = ''


for line in finallist:
  element = line.split()
  try:
    OLDRA = float(element[0])
    OLDDEC = float(element[1])
    ASTK = str(element[4])
    BRA = float(element[5])
    BDEC = float(element[6])
    GLATE = float(element[8])
    GLONGE = float(element[7])
    NEWRA = float(element[9])
    NEWDEC = float(element[10])
    VRAY1 = float(element[11])
    if ('*' in DataSet) and ASTK=='*':
      if (NEWDEC <= (float(DecSlice)+float(DecRange))) and (NEWDEC >= (float(DecSlice)-float(DecRange))):
	
	if abs(ZCoord(GLATE,VRAY1)) < 1500.:
          XG.append(XCoord(GLONGE,VRAY1))
          YG.append(YCoord(GLONGE,VRAY1))
          ZG.append(ZCoord(GLATE,VRAY1))
	#else:
          X.append(XCoord(NEWRA,VRAY1))
          Y.append(YCoord(NEWRA,VRAY1))
          Z.append(ZCoord(NEWDEC,VRAY1))

          XS.append(XSCoord(NEWRA, NEWDEC, SPHRAD))
          YS.append(YSCoord(NEWRA, NEWDEC, SPHRAD))
          ZS.append(ZSCoord(NEWRA, NEWDEC, SPHRAD))
	  Xp.append(XCoord(GLONGE,VRAY1))
          Yp.append(YCoord(GLONGE,VRAY1))
          Zp.append(ZCoord(GLATE,VRAY1))
          
          DECArray.append(NEWDEC)
          RAArray.append(NEWRA)
          STARS.append(float(element[2]))
          FIELDS.append(float(element[3]))
          VRAY1AR.append(VRAY1)
          GLAT.append(GLATE)
          GLONG.append(GLONGE)
	  VRAY1ARp.append(VRAY1)
	  GLATp.append(GLATE)
	  GLONGp.append(GLONGE)
          Count+=1
    if ('+' in DataSet) and ASTK=='+':
      if (NEWDEC <= (float(DecSlice)+float(DecRange))) and (NEWDEC >= (float(DecSlice)-float(DecRange))):
	
        if abs(ZCoord(GLATE,VRAY1)) < 1500.:
	  XG.append(XCoord(GLONGE,VRAY1))
          YG.append(YCoord(GLONGE,VRAY1))
          ZG.append(ZCoord(GLATE,VRAY1))
          
	#else:
          X.append(XCoord(NEWRA,VRAY1))
          Y.append(YCoord(NEWRA,VRAY1))
          Z.append(ZCoord(NEWDEC,VRAY1))       
          XS.append(XSCoord(NEWRA,NEWDEC, SPHRAD))
          YS.append(YSCoord(NEWRA,NEWDEC, SPHRAD))
          ZS.append(ZSCoord(NEWRA,NEWDEC, SPHRAD))
          
          
          DECArray.append(NEWDEC)
          RAArray.append(NEWRA)
          STARS.append(float(element[2]))
          FIELDS.append(float(element[3]))
          VRAY1AR.append(VRAY1)
          GLAT.append(GLATE)
          GLONG.append(GLONGE)
          Xp.append(XCoord(GLONGE,VRAY1))
          Yp.append(YCoord(GLONGE,VRAY1))
          Zp.append(ZCoord(GLATE,VRAY1))
          XSp.append(XSCoord(NEWRA,NEWDEC, SPHRAD))
          YSp.append(YSCoord(NEWRA,NEWDEC, SPHRAD))
          ZSp.append(ZSCoord(NEWRA,NEWDEC, SPHRAD))
          
          XGp.append(XCoord(GLONGE,VRAY1))
          YGp.append(YCoord(GLONGE,VRAY1))
          ZGp.append(ZCoord(GLATE,VRAY1))
          DECArrayp.append(NEWDEC)
          RAArrayp.append(NEWRA)
          STARSp.append(float(element[2]))
          FIELDSp.append(float(element[3]))
          VRAY1ARp1.append(VRAY1)
          GLATp1.append(GLATE)
          GLONGp1.append(GLONGE)
          
        Count+=1
    if ('f' in DataSet) and ASTK=='f':
      if (NEWDEC <= (float(DecSlice)+float(DecRange))) and (NEWDEC >= (float(DecSlice)-float(DecRange))):
	
	if abs(ZCoord(GLATE,VRAY1)) < 1500.:
          XG.append(XCoord(GLONGE,VRAY1))
          YG.append(YCoord(GLONGE,VRAY1))
          ZG.append(ZCoord(GLATE,VRAY1))
        #else:
          X.append(XCoord(NEWRA,VRAY1))
          Y.append(YCoord(NEWRA,VRAY1))
          Z.append(ZCoord(NEWDEC,VRAY1))
          
          XS.append(XSCoord(NEWRA, NEWDEC, SPHRAD))
          YS.append(YSCoord(NEWRA, NEWDEC, SPHRAD))
          ZS.append(ZSCoord(NEWRA, NEWDEC, SPHRAD))
   
	  Xp.append(XCoord(GLONGE,VRAY1))
          Yp.append(YCoord(GLONGE,VRAY1))
          Zp.append(ZCoord(GLATE,VRAY1))
          
          DECArray.append(NEWDEC)
          RAArray.append(NEWRA)
          STARS.append(float(element[2]))
          FIELDS.append(float(element[3]))
	  GLAT.append(GLATE)
          GLONG.append(GLONGE)
          VRAY1AR.append(float(element[11]))
	  VRAY1ARp2.append(VRAY1)
          GLATp2.append(GLATE)
          GLONGp2.append(GLONGE)
          Count+=1
    else:
      pass
  except:
    print 'Invalid'
    pdb.set_trace()
    pass
'''
print Count
if HipCall != '':
    HX,HY,HZ = [],[],[]
    hip = genfromtxt(home+'hiptrunk.dat',delimiter=' ',usecols=(0,1,2,3,4),dtype=('i8','i8','f8','f8','f8'),names=('HIP','VMAG','RA','DEC','PARA'))
    for i in range(len(hip['PARA'])):
      if hip['PARA'][i] < 0.002:
        delete(hip,i,0)
      else:
        hip['PARA'][i] = 1000./(hip['PARA'][i])
      if hip['PARA'][i] > 500:
        delete(hip,i,0)

    HipRA, HipDEC, HipDIST = Hipparcos(RAArray,DECArray,VRAY1AR,hip['RA'],hip['DEC'],hip['PARA'])
    HX = XCoord(HipRA,HipDIST)
    HY = YCoord(HipRA,HipDIST)
    HZ = ZCoord(HipDEC,HipDIST)
    HX = HX[where(ang_sep(HipRA, HipDEC, RAArray, DECArray) < FOV )]
    HY = HY[where(ang_sep(HipRA, HipDEC, RAArray, DECArray) < FOV )]
    HZ = HZ[where(ang_sep(HipRA, HipDEC, RAArray, DECArray) < FOV )]
    for i in range(len(VRAY1AR)):
      for j in range(len(HipDIST)):
        match_count = 0
        if ang_sep(HipRA[j], HipDEC[j], RAArray[i],DECArray[i]) < (FOV):
          star_mag.append(hip['VMAG'][j])
          match_count+=1
      print match_count
      gage_mag.append(average(array(star_mag)))
    print 'Hipparcos Objects:',len(HipRA)

##### Galactic Reference Points
REFX = XCoord(REFRA,REFDIST)
REFY = YCoord(REFRA,REFDIST)
REFZ = ZCoord(REFDEC,REFDIST)
#####


ion()

figure(0)
plot(X, Y, 'bo')
plot(X, Y, 'r')
title('SET2 Equatorial Position')
ylabel('Position (VRU)')
xlabel('Position (VRU)')

figure(1)
hist(DECArray, 10)
xlabel('Declination')
title('Declination Distribution (SET2)')

figure(1)
plot(GLAT, VRAY1AR, 'bo',label='SET1+SET2')
plot(GLATp, VRAY1ARp, 'ro',label='SET3')
title('SET1 + SET2 + SET3')
ylabel('Visual Ray Distance (VRU)')
xlabel('Galactic Latitude (deg)')
legend()

figure(2)
plot(XG, YG, 'ro')
plot([0], [0], 'bo')
title('SET1 + SET2 + SET3 (<5' + r'$\degree$' + ' Galactic Lat)')
ylabel('Position ')
xlabel('Position ')


##### 3D Plot Equatorial Coordinates
fig1 = figure()
ax = fig1.add_subplot(111,projection='3d')
ax.set_xlabel('X (VRU)')
ax.set_ylabel('Y (VRU)')
ax.set_zlabel('Z (VRU)')
ax.set_title('SET1+SET2+SET3 Equatorial Plot')
ax.scatter(Xp,Yp,Zp,s=10,c='b')
ax.plot([400], [0], [0],'g-')
ax.plot([0], [400], [0],'g-')
ax.plot([0], [0], [400],'g-')
ax.scatter([0],[0],[0],s=100,c='g',marker='x')
ax.scatter(X,Y,Z,s=10,c='r')
#####

##### 3D Plot Galactic Coordinates
fig2 = figure()
bx = fig2.add_subplot(111,projection='3d')
bx.set_xlabel('X (VRU)')
bx.set_ylabel('Y (VRU)')
bx.set_zlabel('Z (VRU)')
bx.set_title('SET1+SET2+SET3 Galactic Plot')
bx.scatter(XG,YG,ZG,s=10,c='r')
bx.scatter(Xp,Yp,Zp,s=10,c='b')
bx.plot([400], [0], [0],'g-')
bx.plot([0], [400], [0],'g-')
bx.plot([0], [0], [400],'g-')
bx.scatter([0],[0],[0],s=100,c='g',marker='x')
#bx.scatter(X,Y,Z,s=10,c='r')
#####

##### Hipparcos Plotting
if HipCall != '' and SuperImp != '':
     ax.scatter(HX,HY,HZ,s=10,c='y',marker='s') #Superimpose data onto 3D Equatorial Map of Herschel Gages
elif HipCall != '' and (SuperImp == '' or 'n' in SuperImp):
     fig2 = figure()
     bx = fig2.add_subplot(111,projection='3d')
     bx.set_xlabel('X (Roughly Sirius Dist')
     bx.set_ylabel('Y (Roughly Sirius Dist')
     bx.set_zlabel('Z (Roughly Sirius Dist')
     bx.set_title('3D Equatorial Plot Hipparcos (J2000)')
     bx.scatter(HX,HY,HZ,s=10)
     bx.scatter([0],[0],[0],s=40,c='g',marker='x')
     bx.scatter(REFX,REFY,REFZ,s=20,c='k')
     bx.plot([0,REFX[0]],[0,REFY[0]],[0,REFZ[0]],'g--') #SAG A* AXIS RED LINE
     bx.plot([0,REFX[1]],[0,REFY[1]],[0,REFZ[1]],'r--') #GLAT = 90 AXIS GREEN LINE
     bx.plot([0,REFX[2]],[0,REFY[2]],[0,REFZ[2]],'b--') #GLAT = 0 AXIS BLUE LINE

#####

##### 3D Plot Galactic Coordinates
fig3 = figure()
cx = fig3.add_subplot(111,projection='3d')
cx.set_xlabel('X (VRU)')
cx.set_ylabel('Y (VRU)')
cx.set_zlabel('Z (VRU)')
cx.set_title('3D Galactic Plot (B1950)')
cx.plot(XG,YG,ZG,'r')
cx.scatter(XG,YG,ZG,c='b')
cx.plot([0],[400],[0],'g--')
cx.scatter([0],[0],[0],s=40,c='g',marker='x') #HOME REFERENCE POINT
#####

##### Mollweide Projection Galactic
try:
  fig4 = figure()
  dx = fig4.add_subplot(111,projection='mollweide')
  dx.set_title('Galactic Mollweide Projection')
  dx.set_xlabel('GALACTIC LONGITUDE')
  dx.set_ylabel('GALACTIC LATITUDE')
  dx.grid()
  dx.set_xticklabels(r_[210:360:30, 0:180:30][::-1])
  PR = dx.scatter(Project(GLAT,GLONG)[0],Project(GLAT,GLONG)[1],s=20,c=VRAY1AR,marker='o',cmap=cm.gist_rainbow)
  fig4.colorbar(PR)
  fig4.suptitle('FULL DATA')
except:
  pass

#####


##### Mollweide Projection Equatorial
try:
  print np.sort(RAArray)
  fig5 = figure()
  ex = fig5.add_subplot(111,projection='mollweide')
  ex.set_title('Equatorial Mollweide Projection')
  ex.set_xlabel('RIGHT ASCENSION')
  ex.set_ylabel('DECLINATION')
  ex.grid()
  ex.set_xticklabels(r_[240:360:30, 0:180:30][::-1])
  cbar1 = ex.scatter(Project(DECArray,RAArray)[0],Project(DECArray,RAArray)[1],s=20,c=VRAY1AR,marker='o',cmap=cm.gist_rainbow)
  fig5.colorbar(cbar1)
  fig5.suptitle('FULL DATA')
except:
  pass

#####


##### 2D Plot
try:
  pdb.set_trace()
  fig6 = figure()
  fx = fig6.add_subplot(111)
  fx.set_title('VISUAL RAY LENGTH VS GALACTIC LATITUDE')
  fx.set_xlabel('GALACTIC LAT (DEG)')
  fx.set_ylabel('VISUAL RAY (VRU)')
  #fx.tick_params(axis='x', pad=-10)
  #fx.tick_params(axis='y', pad=-20)
  #fx.grid()
  fx.plot(GLAT, VRAY1AR, 'bo')
  fx.plot(GLATp, VRAY1ARp, 'ro')
  #cbar2 = fx.scatter(GLONG,VRAY1AR, s=20, c=VRAY1AR, marker='o',cmap=cm.gist_rainbow)
  #fig6.colorbar(cbar2)
  fig6.suptitle('FULL DATA')
  #fx.set_adjustable('box')
  #fx.set_aspect(1)
  #fx.set_xlim([0,360])
  #fx.set_ylim([-90,90])
except:
  pass

#####

##### Galactic Cartesian Projection
try:
  fig7 = figure()
  gx = fig7.add_subplot(111)
  gx.set_title('Galactic Cartesian Projection')
  gx.set_xlabel('LONGITUDE')
  gx.set_ylabel('LATITUDE')
  gx.tick_params(axis='x', pad=-10)
  gx.tick_params(axis='y', pad=-20)
  gx.grid()
  cbar3 = gx.scatter(GLONG,GLAT, s=20, c=VRAY1AR, marker='o',cmap=cm.gist_rainbow)
  fig7.colorbar(cbar3)
  fig7.suptitle('FULL DATA')
  gx.set_adjustable('box')
  gx.set_aspect(1)
  gx.set_xlim([0,360])
  gx.set_ylim([-90,90])
except:
  pass
#####

##### Histograms
fig10 = figure()

zx = fig10.add_subplot(221)
zx.hist(VRAY1AR,bins=50,color='g')
zx.set_title('Visual Ray Count')

zx = fig10.add_subplot(222)
zx.hist(FIELDS,bins=[0.,0.25,0.5,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.],color='r',rwidth=.5)
zx.set_title('Field Count')

zx = fig10.add_subplot(223)
zx.hist(STARS,bins=100,color='m')
zx.set_title('Star Count')

zx = fig10.add_subplot(224)
zx.hist(DECArray,bins=amax(DECArray))
zx.set_title('Declination Dispersion')
#zx.hist(GLAT,bins=90)
#zx.set_title('Galactic Lat Dispersion')

fig10.suptitle('SET 1,2')
#####


##### Write Equatorial 3d coordinates to a text file, capped at 'slice' number of gages per file


G = array(GLONG)
A = array(zip(X,Y,Z))
A = A[A[:,0].argsort()]
GLAT = radians(np.array(GLAT))
figure(10)
title('SET1+SET2+SET3')
ylabel('Visual Ray Length (VRU)')
xlabel('Galactic Longitude (degrees)')

for i in range(len(G)):
  if G[i] <= (90./cos(GLAT[i])):
    plot([G[i]],[VRAY1AR[i]],'bo')
    #cloud = open(home+'LONGSLICE/'+'XYZ1.txt','a+')
    #cloud.write("%s " % str(A[i][0])+'	'+str(A[i][1])+'	'+str(A[i][2])+"\n")
  if G[i] > (90./cos(GLAT[i])) and G[i] <= (180./cos(GLAT[i])):
    plot([G[i]],[VRAY1AR[i]],'ro')
    #cloud = open(home+'LONGSLICE/'+'XYZ2.txt','a+')
    #cloud.write("%s " % str(A[i][0])+'	'+str(A[i][1])+'	'+str(A[i][2])+"\n")
  if G[i] > (180./cos(GLAT[i])) and G[i] <= (270./cos(GLAT[i])):
    plot([G[i]],[VRAY1AR[i]], 'go')
    #cloud = open(home+'LONGSLICE/'+'XYZ3.txt','a+')
    #cloud.write("%s " % str(A[i][0])+'	'+str(A[i][1])+'	'+str(A[i][2])+"\n")
  if G[i] > (270./cos(GLAT[i])) and G[i] <= (360./cos(GLAT[i])):
    plot([G[i]],[VRAY1AR[i]],'ko')
    #cloud = open(home+'LONGSLICE/'+'XYZ4.txt','a+')
    #cloud.write("%s " % str(A[i][0])+'	'+str(A[i][1])+'	'+str(A[i][2])+"\n")
raw_input()
count = 0

'''
ion()
figure(0)
plot(GLATp2, VRAY1ARp2, 'ro',label='SET1')
plot(GLATp, VRAY1ARp, 'bo',label='SET2')
plot(GLATp1, VRAY1ARp1, 'go',label='SET3')
axis('tight')
ylabel('Visual Ray Length (VRU)')
xlabel('Galactic Latitude (deg)')
title('Visual Ray vs Galactic Lat (SET1+SET2+SET3)')

legend()

figure(1)
hist(GLAT, 50, color='Blue',label='SET1',histtype='step')
axis('tight')
xlabel('Galactic Latitude (deg)')
title('Galactic Latitude Distribution (SET1+SET2+SET3)')


figure(2)
plot(GLONG, GLAT, 'bo')
axis('tight')
ylabel('Galactic Latitude (deg)')
xlabel('Galactic Longitude (deg)')
title('Galactic Lat vs Galactic Long (SET1+SET2+SET3)')


figure(3)
plot(RAArray, DECArray, 'bo')
axis('tight')
xlabel('Right Ascension (deg)')
ylabel('Declination (deg)')
title('Declination vs Right Ascension (SET1+SET2+SET3)')

fig4 = figure()
ax = fig4.add_subplot(111,projection='mollweide')
ax.set_title('Galactic Mollweide Projection')
ax.set_xlabel('GALACTIC LONGITUDE')
ax.set_ylabel('GALACTIC LATITUDE')
ax.grid(True)
#ax.set_xticklabels(r_[210:360:30, 0:180:30][::-1])
#PR = ax.scatter(Project(GLAT,GLONG)[0],Project(GLAT,GLONG)[1],s=20,c=VRAY1AR,marker='o',cmap=cm.gist_rainbow)
ax.scatter(Project([0,0,0,0,0,0,0],[0,30,60,90,120,150,180])[0],Project([0,0,0,0,0,0,0],[0,30,60,90,120,150,180])[1])
#fig4.colorbar(PR)
ax.axis('tight')
fig4.suptitle('FULL DATA')

raw_input()

'''
count = 0
A = array(zip(XG,YG,ZG))
cloud = open(home+'XYZ4/'+'XYZ.txt','w+')
ID = 0
A = A[A[:,0].argsort()]
for i in range(len(A)):
  cloud.write("%s " % str(A[i][0])+'	'+str(A[i][1])+'	'+str(A[i][2])+"\n")
  

#for i in range(len(A)):
#  count+=1
#  if(count == slice):
#    count = 0
#    ID+=1
#    cloud = open(home+'XYZ4/'+'XYZ'+str(ID)+'.txt','w+')
#  cloud.write("%s " % str(A[i][0])+'	'+str(A[i][1])+'	'+str(A[i][2])+"\n")

#####

##### Statistics
AVG=average(DECArray)
SIG=std(DECArray)

print 'Objects loaded in: ' +str(Count)
print 'Slicing declination, range: ' +str(DecSlice)+', '+str(DecRange)
print 'Average data declination, standard deviation:',AVG,SIG

#####
'''
