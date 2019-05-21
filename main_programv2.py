"""

author: Timo Wortelboer

written for Python 3.6.5

to be used in accordance with the Creative Commons Share-Alike licence, so all derivative works
should be published with the same kind of licence and should be properly accredited

"""

from meshingv3 import *
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pylab

import os

import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.tri as tri

import math

try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle


def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


##########################################################################
##### Test level-set functions ###########################################
##########################################################################

def cos_times_sin(p):
    (x, y) = p
    (a, b) = (1.5, 0.5)
    if isinstance(x, list) or isinstance(x, np.ndarray):
        f = []
        for i, xval in enumerate(x):
            f.append(np.cos(a * np.pi * xval) * np.sin(b * np.pi * y[i]))
        return f
    else:
        return np.cos(a * np.pi * x) * np.sin(b * np.pi * y)


def linear(p):
    (x, y) = p
    (a, b, c) = (3, 2, 1)
    if isinstance(x, list) or isinstance(x, np.ndarray):
        f = []
        for i, xval in enumerate(x):
            f.append(a * xval + b * y[i] + c)
        return f
    else:
        return a * x + b * y + c


def circle(p):
    (x, y) = p
    (a, b) = (1.1, 1.3)
    if isinstance(x, list) or isinstance(x, np.ndarray):
        f = []
        for i, xval in enumerate(x):
            f.append((xval / a) ** 2 + (y[i] / b) ** 2 - 1)
        return f
    else:
        return (x / a) ** 2 + (y / b) ** 2 - 1


def star(p):
    (x, y) = p
    if isinstance(x, list) or isinstance(x, np.ndarray):
        f = []
        for i in range(len(x)):
            if x[i] == 0:
                if y[i] > 0:
                    angle = np.pi / 2
                else:
                    angle = -np.pi / 2
            else:
                angle = np.arctan(y[i] / x[i])
            f.append((3. / 5) + (1. / 4) * np.cos(6 * angle) ** 2 - (x[i] ** 2 + y[i] ** 2))
        return f
    else:
        if x == 0:
            if y > 0:
                angle = np.pi / 2
            else:
                angle = -np.pi / 2
        else:
            angle = np.arctan(y / x)
        return ((3. / 5) + (1. / 4) * np.cos(6 * angle) ** 2 - (x ** 2 + y ** 2))


def dumbbell(p):
    (x, y) = p
    (ax, ay) = (-0.5, 0)
    (bx, by) = (0.5, 0)
    (c1, c2, c3) = (0.3, 0.25, 0.1)

    if isinstance(x, list) or isinstance(x, np.ndarray):
        f = []
        for i in range(len(x)):
            dxa = np.sqrt((ax - x[i]) ** 2 + (ay - y[i]) ** 2)
            dxb = np.sqrt((bx - x[i]) ** 2 + (by - y[i]) ** 2)
            if ax < x[i] and x[i] < bx:
                f.append(max(c1 - dxa, c2 - dxb, c3 - abs(y[i])))
            else:
                f.append(max(c1 - dxa, c2 - dxb))
        return f
    else:
        dxa = np.sqrt((ax - x) ** 2 + (ay - y) ** 2)
        dxb = np.sqrt((bx - x) ** 2 + (by - y) ** 2)
        if ax < x and x < bx:
            return max(c1 - dxa, c2 - dxb, c3 - abs(y))
        else:
            return max(c1 - dxa, c2 - dxb)


##########################################################################
##### Quality measurement ################################################
##########################################################################


class qualitycontainer:
    # This object keeps the quality measurement data and can be saved to a file

    def __init__(self, methodname, hstep):
        self.methodname = methodname  # name of the applied method
        self.hstep = hstep  # typical grid edge length
        self.before = False  # holds the quality data from before adjsuting the mesh
        self.fnames = []  # array of names of the used test-functions
        self.after = []  # array of the data after applying the method on the different functions

    def save(self, hstep):
        # This function saves the qualitycontainer to a file based on hstep and methodname
        directory = 'Saved_qualitydata/' + str(round(hstep, 8)) + '/'
        try:
            os.makedirs(directory)
        except OSError:
            None
        savename = directory + self.methodname.replace(" ", "_") + '.pkl'
        save_object(self, savename)
        print('Saved the qdata object with name \n  ' + savename)


def qualityplot(allmethodnames, hstep, fnames=False):
    # This function fetches the qualitycontainers from the saved files on the disk and plots bar charts of the quality
    # The fetching process is based on the given allmethodnames and prints a warning if a file for some methodname does not exist
    qualityobjects = []
    for q, methodname in enumerate(allmethodnames):
        loadname = 'Saved_qualitydata/' + str(round(hstep, 8)) + '/' + methodname.replace(" ", "_") + '.pkl'
        try:
            with open(loadname, 'rb') as input:
                qdata = pickle.load(input)
                qualityobjects.append(qdata)
        except:
            print('Quality dataset of method <' + methodname + '> is not detected.')
    if fnames == False and len(qualityobjects) > 0:
        fnames = qualityobjects[0].fnames
    else:
        print('Error, no quality object files found')

    nfuns = len(fnames)
    nbars = len(qualityobjects) + 1
    width = 0.7 / nbars
    colors = ['xkcd:salmon', 'xkcd:light blue', 'xkcd:green', 'xkcd:yellow', 'xkcd:lilac', 'xkcd:royal blue']
    plegend = ['Before']

    for m, beforedataset in enumerate(qualityobjects[0].before):  # we plot a figure for each measuretype with index m
        measurename = beforedataset[0]  # every datapackage in qualityobjects is a list of the measurename and the data
        beforedata = beforedataset[1]
        fig = pyplot.figure(figsize=(8, 8))
        for f, fun in enumerate(fnames):
            if nbars % 2 == 0:
                plist = [pyplot.bar(f - nbars / 2 * width, beforedata, width, color=colors[0], align='edge')]
                for q, qdata in enumerate(qualityobjects):
                    plegend.append(qdata.methodname)
                    if fun in qdata.fnames:
                        findex = qdata.fnames.index(fun)
                        plist.append(pyplot.bar(f + (q + 1 - nbars / 2) * width, qdata.after[findex][m][1], width,
                                                color=colors[q + 1], align='edge'))
            else:
                plist = [pyplot.bar(f - nbars // 2 * width, beforedata, width, color=colors[0], align='center')]
                for q, qdata in enumerate(qualityobjects):
                    plegend.append(qdata.methodname)
                    if fun in qdata.fnames:
                        findex = qdata.fnames.index(fun)
                        plist.append(pyplot.bar(f + (q + 1 - nbars // 2) * width, qdata.after[findex][m][1], width,
                                                color=colors[q + 1], align='center'))

        pyplot.title(measurename)
        pyplot.xticks(range(len(qualityobjects[0].fnames)), qualityobjects[0].fnames)
        pyplot.ylabel(measurename)
        xfmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
        xfmt.set_powerlimits((-3, 4))
        pyplot.gca().yaxis.set_major_formatter(xfmt)
        pyplot.legend(plist, plegend)

        savename = 'Saved_figures/' + str(round(hstep, 8)) + '/Quality_' + measurename.replace(" ", "_") + '.png'
        fig.savefig(savename, dpi=fig.dpi, bbox_inches="tight")
    return


def simplesol(z, thresh, qdata, relaxthresh, newtriang=False, flip=False, distplot=False, simplepointrelax=False,
              simplegridtozeropoint=False, redistribute=False):
    # This function creates a mesh and applies the given methods that are given as True in the input with the test function z
    # It adds the resulting quality data in the qualitycontainer given in the input

    m = mesh(hstep, thresh)
    m.initzeroinfo(z)

    mname = qdata.methodname.replace(" ", "_")
    directory = 'Saved_figures/' + str(round(qdata.hstep, 8)) + '/' + mname + '/'
    directory2 = 'Saved_figures/' + str(round(qdata.hstep, 8)) + '/'
    try:
        os.makedirs(directory)
    except OSError:
        None
    savestring = directory2 + 'Before_' + z.__name__ + '_mesh.png'
    if not os.path.isfile(savestring):  # plots and saves only if the file does not yet exist
        m.plotzeroz(savestring)
    checkstring = directory2 + z.__name__ + '_contour.png'
    if not os.path.isfile(checkstring):
        m.plotz(z, directory2)
    pyplot.close()

    if qdata.before == False:  # calculates the quality before adjusting the mesh only if it does not yet exist in the qualitycontainer
        qdata.before = m.skewness()
        qdata.before.extend(m.sizes())

    if simplegridtozeropoint:  # facilitates the movement of meshpoints when closer than 0.3h to a zeropoint
        m.simplegridtozeropoint(thresh, z=z)
        m.initzeroedgesinfo(z)
        m.updatezeroelementinfo(z)

    if newtriang and flip:  # facilitates the cut method (newtriang) and flip method
        m.newtriangwithflip()
        m.updatezeroelementinfo(z)
    elif newtriang:
        m.newtriang()
        m.updatezeroelementinfo(z)

    if simplepointrelax:  # facilitates a pointwise mesh-relaxation method
        m.pointrelax(relaxthresh)

    if shortestpath:  # facitates the shortest path method with some additional plots
        (boundarytype, boundaryplist) = m.setupboundaryplist()
        # savename = qdata.methodname
        # savename.replace(" ", "_")
        # print('---Plotting 1')
        # m.plotzeroz(directory + 'After_setuppoints_' + mname + '_' + z.__name__ + '_adjmesh.png', dist=distplot,
        #           plotchosenboundary=True)
        zeropathlist = m.findshortestpath(boundaryplist, boundarytype)
        # print('---Plotting 2')
        # m.plotzeroz(directory + 'After_shortestpath_' + mname + '_' + z.__name__ + '_adjmesh.png', dist=distplot)

        m.movetolevelset(zeropathlist)
        if redistribute:
            m.redistribute_projectedpoints(zeropathlist, z)
        m.finallevelsetinfoupdate()
        # print('---Plotting 3')
        # m.plotzeroz(directory + 'After_movetolevelset_' + mname + '_' + z.__name__ + '_adjmesh.png', dist=False, mode=2)

    if fixedpointrelaxation:  # facilitates a fixed-point iteration scheme for mesh relaxation (divergent)
        m.fixedpoint()
        print('---Plotting 4')
        m.plotzeroz(directory + 'After_fixedpoint_' + mname + '_' + z.__name__ + '_adjmesh.png', dist=False, mode=2)

    if eulerrelaxation:  # facilitates a iterative mesh-relaxation scheme that seems to always converge
        m.eulerforward()
        print('---Plotting 4')
        m.plotzeroz(directory + 'After_euler_' + mname + '_' + z.__name__ + '_adjmesh.png', dist=False, mode=2)

    quality = m.skewness()
    quality.extend(m.sizes())
    qdata.after.append(quality)

    savename = qdata.methodname
    savename.replace(" ", "_")
    print('---Plotting 4')
    m.plotzeroz(directory + 'After_' + mname + '_' + z.__name__ + '_adjmesh.png', dist=distplot)
    return m


def simplemeshfit(thresh, methodname, hstep, relaxthresh, farr=[linear, circle, star, dumbbell],
                  newtriang=False, flip=False, distplot=False, simplepointrelax=False, simplegridtozeropoint=False,
                  redistribute=False):
    # This function creates the qualitycontainer for a specific hstep size and method and applies the method through the simplesol function for each test-function in farr
    print('hstep = ' + str(hstep))
    qdata = qualitycontainer(methodname, hstep)

    qdata.fnames.append(farr[0].__name__)
    print('Applying method <' + methodname + '> for function ' + farr[0].__name__)
    mf = simplesol(farr[0], thresh, qdata, relaxthresh, newtriang=newtriang, flip=flip, distplot=distplot,
                   simplepointrelax=simplepointrelax, simplegridtozeropoint=simplegridtozeropoint,
                   redistribute=redistribute)
    qdata.save(hstep)
    return mf


##########################################################################
##### Program execution ##################################################
##########################################################################

print('-----------------------------------------------------')

allmethodnames = ['Cut method', 'Cut and flip method', 'Shortest path fit', 'Shortest path fit with redistribution',
                  'Shortest path fit with Euler forward relaxation']

# Choose and uncomment the right set of method from the allmethodnames list below (and comment out the rest)


# methodstoapply = []
# methodstoapply = ['Cut method']
# methodstoapply = ['Cut and flip method']
# methodstoapply = ['Shortest path fit']
methodstoapply = ['Shortest path fit with redistribution']
# methodstoapply = ['Shortest path fit with fixed-point relaxation']
# methodstoapply = ['Shortest path fit with Euler forward relaxation']
# methodstoapply = ['Cut method','Cut and flip method']
# methodstoapply = ['Shortest path fit', 'Shortest path fit with redistribution',
                  #'Shortest path fit with Euler forward relaxation']
# methodstoapply = allmethodnames


# Choose the right set of sizes to apply from the list below

# sizestoapply = [1 / 8, 1 / 16]
sizestoapply = [1/8]
# sizestoapply = [1/16]


# Choose the right set of test functions to apply below

# farr = [cos_times_sin, linear, circle, star, dumbbell]
# farr = [cos_times_sin,linear,circle]
# farr = [star,dumbbell]
# farr = [cos_times_sin]
# farr = [circle]
farr = [linear]
# farr = [linear, circle, star, dumbbell]


for hstep in sizestoapply:  # sets the right parameters for each method and applies the method through simplemeshfit
    for methodname in methodstoapply:
        newtriang = False
        flip = False
        simplepointrelax = False
        distplot = False
        simplegridtozeropoint = False
        shortestpath = False
        fixedpointrelaxation = False
        eulerrelaxation = False
        redistribute = False
        if methodname == 'Cut method':
            simplegridtozeropoint = True
            newtriang = True
        elif methodname == 'Cut and flip method':
            simplegridtozeropoint = True
            newtriang = True
            flip = True
        elif methodname == 'Cut, flip, point relaxation':
            simplegridtozeropoint = True
            newtriang = True
            flip = True
            simplepointrelax = True
        elif methodname == 'Shortest path fit':
            shortestpath = True
            distplot = True
        elif methodname == 'Shortest path fit with redistribution':
            shortestpath = True
            distplot = True
            redistribute = True
        elif methodname == 'Shortest path fit with fixed-point relaxation':
            shortestpath = True
            distplot = True
            fixedpointrelaxation = True
        elif methodname == 'Shortest path fit with Euler forward relaxation':
            shortestpath = True
            distplot = True
            eulerrelaxation = True

        relaxthresh = hstep / 100  # minimum displacement threshold for relaxing gridpoints in point relaxation method
        thresh = 0.3 * hstep  # maximum distance threshold for moving gridpoints in simplegridtozeropoint function
        simplemeshfit(thresh, methodname, hstep, relaxthresh, farr=farr,
                      newtriang=newtriang, flip=flip, distplot=distplot, simplepointrelax=simplepointrelax,
                      simplegridtozeropoint=simplegridtozeropoint, redistribute=redistribute)

        pyplot.close('all')  # The plots can be found in the Saved_figures directory

    print('Calculating and plotting the quality measurement data')
    qualityplot(allmethodnames, hstep, fnames=False)

    pyplot.close('all')  # The plots can be found in the Saved_figures directory
