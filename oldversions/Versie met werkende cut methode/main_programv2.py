from meshingv2 import *
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab

import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.tri as tri
#matplotlib.rcParams['text.usetex'] = True
#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

hstep = 1/16
threshold = 0.05 * hstep

print('-----------------------------------------------------')

def cos_times_sin(p):               # the level-set function
    (x,y) = p
    (a,b) = (1.5,0.5)
    if isinstance(x,list) or isinstance(x,np.ndarray):
        f = []
        for i,xval in enumerate(x):
            f.append( np.cos(a*np.pi*xval)*np.sin(b*np.pi*y[i]) )
        return f
    else:
        return np.cos(a*np.pi*x)*np.sin(b*np.pi*y)

def linear(p):
    (x,y)   = p
    (a,b,c) = (3,2,1)
    if isinstance(x,list) or isinstance(x,np.ndarray):
        f = []
        for i,xval in enumerate(x):
            f.append( a*xval+b*y[i]+c )
        return f
    else:
        return a*x+b*y+c

def circle(p):
    (x,y) = p
    (a,b) = (1.1,1.3)
    if isinstance(x,list) or isinstance(x,np.ndarray):
        f = []
        for i,xval in enumerate(x):
            f.append( (xval/a)**2+(y[i]/b)**2-1 )
        return f
    else:
        return (x/a)**2+(y/b)**2-1

def star(p):
    (x,y) = p
    if isinstance(x,list) or isinstance(x,np.ndarray):
        f = []
        for i in range(len(x)):
            if x[i] == 0:
                if y[i] > 0:
                    angle = np.pi/2
                else:
                    angle = -np.pi/2
            else:
                angle = np.arctan(y[i]/x[i])
            f.append( (3./5) + (1./4)*np.cos(6*angle)**2 - (x[i]**2+y[i]**2) )
        return f
    else:
        if x == 0:
            if y > 0:
                angle = np.pi/2
            else:
                angle = -np.pi/2
        else:
            angle = np.arctan(y/x)
        return ( (3./5) + (1./4)*np.cos(6*angle)**2 - (x**2+y**2) )

def dumbbell(p):
    (x,y)      = p
    (ax,ay)    = (-0.5,0)
    (bx,by)    = (0.5,0)
    (c1,c2,c3) = (0.3,0.25,0.1)
    
    if isinstance(x,list) or isinstance(x,np.ndarray):
        f = []
        for i in range(len(x)):
            dxa = np.sqrt( (ax-x[i])**2+(ay-y[i])**2 )
            dxb = np.sqrt( (bx-x[i])**2+(by-y[i])**2 )
            if ax<x[i] and x[i]<bx:
                f.append(max( c1-dxa, c2-dxb,c3-abs(y[i]) ))
            else:
                f.append(max( c1-dxa, c2-dxb ))
        return f
    else:
        dxa = np.sqrt( (ax-x)**2+(ay-y)**2 )
        dxb = np.sqrt( (bx-x)**2+(by-y)**2 )
        if ax<x and x<bx:
            return max( c1-dxa, c2-dxb, c3-abs(y) )
        else:
            return max( c1-dxa, c2-dxb )

def qualityplot(data,names):
    ind = np.arange(len(data))
    width = len(data)*(0.35/5)
    
    for j in range(len(names)):
        fig = pyplot.figure(figsize=(7,7))
        for i in range(len(data)):
            p1 = pyplot.bar(i-width/2, data[i][1][j], width, color='xkcd:salmon')
            p2 = pyplot.bar(i+width/2, data[i][2][j], width, color='xkcd:light blue')
        pyplot.title(names[j])
        xticks = []
        for i in range(len(data)):
            xticks.append(data[i][0].replace("_"," "))
        pyplot.xticks(ind,xticks)
        pyplot.ylabel(names[j])
        pyplot.legend((p1[0],p2[0]),('Before','After'))
        string = 'Saved_figures/' + str(names[j]).replace(" ", "_") + '.png'
        fig.savefig(string, dpi=fig.dpi, bbox_inches = "tight")
    return

def simplesol(z,thresh, newtriang = False):
    m = mesh(hstep,thresh)
    m.initzeroedgesinfo(z)
    m.updatezeroelementinfo(z)
    m.findzero(z)
    m.plotzeroz('Saved_figures/' + z.__name__ + '_mesh.png')
    m.plotz(z)
    pyplot.close()

    [avgskewness, maxskewness, stdskewness] = m.skewness()
    [minsize, maxsize, stdsize] = m.sizes()
    before = [avgskewness, maxskewness, stdskewness, minsize, maxsize, stdsize]

    m.simplegridtozeropoint(z,thresh)

    if newtriang:
        m.newtriang()
        m.updatezeroelementinfo(z)

    m.findzero(z)

    [avgskewness, maxskewness, stdskewness] = m.skewness()
    [minsize, maxsize, stdsize] = m.sizes()
    after = [avgskewness, maxskewness, stdskewness, minsize, maxsize, stdsize]

    m.plotzeroz('Saved_figures/' + z.__name__ + '_adjmesh.png')
    return [before, after]

def simplemeshfit(farr,thresh, newtriang = True):
    print('hstep = '+str(hstep))
    if isinstance(farr,list):
        soldata = []
        for i in range(len(farr)):
            soldata.append([farr[i].__name__])
        for i in range(len(farr)):
            print('Calculating mesh for function '+farr[i].__name__)
            soldata[i].extend(simplesol(farr[i],thresh,newtriang = True))
        return soldata
    else:
        print('Calculating mesh for function '+farr.__name__)
        return simplesol(farr,thresh)

thresh = 0.3 * hstep
data = simplemeshfit([cos_times_sin,linear,circle,star,dumbbell],thresh, newtriang = True)

#pyplot.show()
pyplot.close('all')

names = ['Average skewness', 'Maximum skewness', 'Standard deviation of skewness',
         'Minimum triangle size', 'Maximum triangle size', 'Standard deviation of triangle size']
qualityplot(data,names)

#pyplot.show()
pyplot.close('all')


