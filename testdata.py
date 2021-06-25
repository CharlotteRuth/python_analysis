import pynbody
import pylab
s = pynbody.load('/home/christenc/Code/testdata/g15784.lr.01024.gz')
#s = pynbody.load('/home/christensen/Code/testdata/g15784.lr.01024.gz')
h = s.halos()
h1 = h[1]
pynbody.analysis.halo.center(h1,mode='hyb')
h1.physical_units()
s.physical_units()
pynbody.plot.image(h1.g, width=20, cmap='Blues', clear = True) #,threaded = False)
pynbody.plot.image(h1.g, width='20 kpc', cmap='Blues', clear = True) #,threaded = False)
pynbody.plot.image(h1.g,width = 20,x1 = -5,y1 = -5, y2 = 10, cmap='Blues', clear = True) #,threaded = False)

