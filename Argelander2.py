from pylab import *
import pdb

data = [15.1,30.3,29.2,20.5,27.4,10.7,10.5,33.1,14.9,9.6,19.9,21.8,12.9,22.2,22.8,37.9,35.0,54.7,62.8,56.8,49.4,31.3,39.6,27.9,23.4,22.9,15.9,25.2,18.4,15.1,10.5,14.4,16.5,15.3,15.7,17.1,18.9,21.5,20.4,21.1,17.3,20.6,27.6,15.4,27.5,33.5,40.4,31.0]

ra = arange(0,270,5)

hist(data,bins=len(data))
title('Herschel vs Argelander Star count per 1 deg')
xlabel('Star Count / 1 degree')
ylabel('Counts')
pdb.set_trace()
show()
