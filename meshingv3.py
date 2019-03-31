"""

author: Timo Wortelboer

written for Python 3.6.5

to be used in accordance with the Creative Commons Share-Alike licence, so all derivative works
should be published with the same kind of licence and should be properly accredited

"""

import numpy as np

import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.tri as tri

import sys

import math

import matlab
import matlab.engine as me
from mpl_toolkits.mplot3d import Axes3D

import colorsys

from copy import copy


##########################################################################
##### Some useful functions ##############################################
##########################################################################


def create_mesh(x,y,h):
    # This function checks if the mesh already exists otherwise creates one using Matlab
    matlab_engine= me.start_matlab()
    x = matlab_engine.eval(np.array2string(x))
    y = matlab_engine.eval(np.array2string(y))
    h = matlab_engine.eval(np.array2string(h))
    matlab_engine.workspace['x'] = x
    matlab_engine.workspace['y'] = y
    matlab_engine.workspace['h'] = h
    matlab_engine.mesh_gen(nargout=0)
    elements = matlab_engine.workspace['elements']
    elements = np.array(elements).astype(int)
    edges = matlab_engine.workspace['edges']
    edges = np.array(edges).astype(int)
    points = matlab_engine.workspace['points']
    points = np.array(points)        
    matlab_engine.quit()
    return elements, edges, points

def pythagoras(p1,p2):
    # Calculates the distance between two points and/or sets of coordinates
    if isinstance(p1,point):
        p1coords = np.array( p1.coordinates )
    else:
        p1coords = np.array( p1 )
    if isinstance(p2,point):
        p2coords = np.array( p2.coordinates )
    else:
        p2coords = np.array( p2 )
    return np.linalg.norm( p1coords - p2coords )

def zerodistmeasure( angle, distance, h ):
    # Sets up the distance measure of the shortest path method
    return ( 1 + angle ) * ( 1 + distance/h ) - 1


def zzerointerp(p1,p2,zfun):
    # Finds the zeropoints along zeroedges by linear interpolation
    [x1,y1] = p1
    [x2,y2] = p2
    zp1 = zfun([x1,y1])
    zp2 = zfun([x2,y2])
    if x1 == x2:
        x0 = x1
    else:
        x0 = x1+(x2-x1)*abs(zp1/(zp1-zp2))
    if y1 == y2:
        y0 = y1
    else:
        y0 = y1+(y2-y1)*abs(zp1/(zp1-zp2))
    return [x0,y0]

def edgescommonp(e1,e2):
    # Finds the point p that two edges have in common (if any)
    for p in e1.points:
        if p in e2.points:
            return p

def angle(a1,a2):
    # Finds the angle between two edges and/or absolute angles
    if isinstance(a1,edge):
        a1 = a1.absangle
    if isinstance(a2,edge):
        a2 = a2.absangle
    dangle = a1 - a2
    return min( abs( dangle ), abs( dangle - np.pi ), abs( dangle + np.pi ) )

def project( point, line ):
    # Calculates the orthogonal projection of a point on a line, or returns the closest boundary point of the line
    (p1,p2)     = line.points
    linevec     = np.array( np.subtract( p1.coordinates, p2.coordinates ) )
    pvec        = np.array( np.subtract( point.coordinates, p2.coordinates ) )
    lineveclen  = np.dot( linevec, linevec )
    if lineveclen > 0.001*point.hstep:
        pproj       = ( np.dot( pvec, linevec ) / np.dot( linevec, linevec ) ) * linevec
        newp        = list( np.add( p2.coordinates, pproj ) )
        minx        = min( p1.coordinates[0], p2.coordinates[0] )
        maxx        = max( p1.coordinates[0], p2.coordinates[0] )
        miny        = min( p1.coordinates[1], p2.coordinates[1] )
        maxy        = max( p1.coordinates[1], p2.coordinates[1] )
        if minx < newp[0] < maxx and miny < newp[1] < maxy:
            return newp
    if pythagoras( point, p1 ) < pythagoras( point, p2 ):
        return p1
    else:
        return p2


##########################################################################
##### Point class ########################################################
##########################################################################

class point:
    # Builds a point class consisting of coordinates and properties
    
    def __init__(self,x,y,mesh):
        
        self.coordinates        = [x,y]
        self.nbedges            = []
        self.nbelements         = []
        self.nbpoints           = []
        self.levelsetedges      = []
        self.zeropoint          = False
        self.iszeropoint        = False
        self.ismoved            = False
        self.isboundary         = False
        self.iszeropath         = False
        self.iscorner           = False
        self.ischosenboundary   = False
        self.isconsidered       = False
        self.isrelaxpoint       = False
        self.hstep              = mesh.hstep
        self.nbzerodist         = math.inf
        mesh.points.append(self)

    def __repr__(self):
        # Adds a string representation with some basic info
        
        #string  = 'Point object with \n'
        string  = 'Coordinates: x = ' + str(self.coordinates[0]) + '\n'
        string += '             y = ' + str(self.coordinates[1]) + '\n'
        #string += 'Number of incident edges = ' + str(len(self.nbedges)) + '\n'
        #string += 'Number of indicent elements = ' + str(len(self.nbedges))
        #string = 'Point'
        return string

    def __eq__(self,other):
        # Decides when two points are the same object or have the same coordinates
        if isinstance( other, point ):
            if self is other:
                return True
            else:
                return False
        else:
            if isinstance( other, zeropoint ):
                other = other.coordinates
            if self.coordinates[0] == other[0] and self.coordinates[1] == other[1]:
                return True
            else:
                return False

    def delete( self, mesh ):
        # Deletes a point and cleans up properties of nbpoints and nbedges
        for e in self.nbedges:
            e.delete(mesh)
        mesh.points = [p for p in mesh.points if p is not self]
        del self
        return

    def move( self, p, mesh, zfun = False, zpmove = False, updateinfo = True ):
        # Performs the necessary operations to move a point to new coordinates p

        if isinstance( p, point ):
            p = p.coordinates
        self.coordinates = p
        self.ismoved     = True

        if zpmove:                    # use this part if the levelsetlines info and zerotypes is important for the rest of the algorithm;
            self.iszeropoint = True   # it cleans up the level-set propertiesof the surrounding elements
            for el in self.nbelements:
                if el.zerotype in (1,2,3):
                    for levelsetline in el.levelsetlines:
                        if p in levelsetline.points:
                            op = levelsetline.otherp( p )
                            if isinstance( op, zeropoint ):
                                levelsetline.points = ( op, self )
                                p.nbedges.remove(levelsetline)
                            else:
                                levelsetline.delete( mesh )
            p.delete( mesh )
        
        if updateinfo:                # use this part to update the info and zeroinfo of surrounding elements and edges after moving
            for e in self.nbedges:    # without altering structure (no level-set lines are deleted)
                e.length = e.calclength()
                if zfun != False:
                    e.findzero( zfun, mesh )
            for el in self.nbelements:
                el.size     = el.calcsize()
                el.skewness = el.calcskew()
                if zfun != False:
                    el.findzero( zfun, mesh )
            return

    def relaxposition(self,relaxthresh):
        # Does not give good results, it relaxes the position of each mesh point by moving it to the average of surrounding points
        coords = np.array([0,0])
        for p in self.nbpoints:
            coords = np.add( coords, p.coordinates )
        coords = coords / len(self.nbpoints)
        if self.isboundary:
            if self.iscorner:
                return False
            else:
                coords = np.dot( coords - self.coordinates, self.freedomvect ) * self.freedomvect + self.coordinates

        if pythagoras( coords, self.coordinates ) > relaxthresh:
            self.move( coords )
            return True
        else:
            return False

    def relaxpositionv2(self,relaxthresh,uthresh = 0.01):
        # Tries to use fixed-point iteration to relax each mesh point individually, but diverges
        
        h    = self.hstep
        z    = np.array(self.coordinates)
        r    = []
        for p in self.nbpoints:
            r.append( p.coordinates )
        n    = len(r)
        c    = np.sum(r,axis=0)/n
        
        usum = 0
        for i in r:
            usum += pythagoras(z,i)
        d    = np.subtract( c,z )
        u    = d[0]**2 + d[1]**2 - (2*h/n)*usum
        usum = 0
        for i in r:
            usum += ( z[0]-i[0] ) / pythagoras(z,i)
        ux   = 2*( z[0]-c[0] ) - (2*h/n)*usum
        usum = 0
        for i in r:
            usum += ( z[1]-i[1] ) / pythagoras(z,i)
        uy   = 2*( z[1]-c[1] ) - (2*h/n)*usum
        
        dummy = 1
        while ux**2 + uy**2 > uthresh and dummy < 10000:
            dummy += 1

            usum = 0
            for i in r:
                usum += ( (z[1]-i[1])**2 ) / ( pythagoras(z,i)**3 )
            uxx = (-2*h/n)*usum + 2

            usum = 0
            for i in r:
                usum += ( (z[0]-i[0])**2 ) / ( pythagoras(z,i)**3 )
            uyy = (-2*h/n)*usum + 2

            usum = 0
            for i in r:
                usum +=  (z[0]-i[0])*(z[1]-i[1]) / ( pythagoras(z,i)**3 )
            uxy = (2*h/n)*usum

            J = np.matrix([[uxx,uxy],[uxy,uyy]])
            Jinv = np.linalg.inv(J)

            z = (z - np.matmul( Jinv, np.array([ux,uy]) )).T
            
            usum = 0
            for i in r:
                usum += pythagoras(z,i)
            u    = np.dot((c-z).T,c-z) - (2*h/n)*usum
            usum = 0
            for i in r:
                usum += ( z[0]-i[0] ) / pythagoras(z,i)
            ux   = 2*( z[0]-c[0] ) - (2*h/n)*usum
            usum = 0
            for i in r:
                usum += ( z[1]-i[1] ) / pythagoras(z,i)
            uy   = 2*( z[1]-c[1] ) - (2*h/n)*usum

        if n == 10000:
            print('newton mehod did not converge')

        if pythagoras( z, self.coordinates ) > relaxthresh:
            self.move( z )
            return True
        else:
            return False

class zeropoint(point):
    # Sets up a seperate class for zero points
    def __init__(self,x,y,mesh):
        self.coordinates    = [x,y]
        self.nbedges        = []
        self.nbelements     = []
        self.nbpoints       = []
        self.levelsetedges  = []
        self.onedge         = False
        self.onelement      = False
        self.isboundary     = False
        self.iszeropath     = False
        self.iscorner       = False
        self.isgrid         = False
        self.gridpoint      = False
        self.hstep          = mesh.hstep
        mesh.zeropoints.append(self)

    def delete( self, mesh ):
        # Deletes a zeropoint also cleaning up the surrounding mesh properties
        for e in self.nbedges:
            e.delete(mesh)
        mesh.zeropoints = [p for p in mesh.zeropoints if p is not self]
        del self
        return


##########################################################################
##### Edge class #########################################################
##########################################################################


class edge:
    # Builds an edge class consisting primarily of two points
    def __init__(self,p1,p2,mesh):
        self.points     = ( p1, p2 )
        for p in self.points:
            p.nbedges += (self,)
        p1.nbpoints.append(p2)
        p2.nbpoints.append(p1)
        self.nbelements     = []
        self.pprojected     = []
        self.length         = self.calclength()
        self.islevelset     = False
        self.levelsetedge   = False
        self.iszeroedge     = False
        self.isboundary     = False
        self.hstep          = mesh.hstep
        self.absangle       = self.calcangle()
        self.zmid           = self.findmid()
        mesh.edges.append(self)

    def __repr__(self):
        # Sets up a string method showing some basic properties
        string  = 'Edge object with \n'
        string += 'Number of incident elements = ' + str(len(self.nbelements)) + '\n'
        return string

    def __eq__(self,other):
        # Decides if some other line is equal based on point objects
        (s1,s2) = self.points
        if isinstance(other,edge):
            if self is other:
                return True
            else:
                return False
        elif isinstance(other,bool):
            return False
        else:
            (o1,o2) = other
            if ( (s1 == o1) and (s2 == o2) ) or ( (s1 == o2) and (s2 == o1) ):
                return True
            else:
                return False

    def delete(self,mesh,delobject = True):
        #Deletes the edge object and cleans up the surrounding mesh properties in the process
        for el in self.nbelements:
            el.delete(mesh)
        (p1,p2) = self.points
        p1.nbpoints = [p for p in p1.nbpoints if p is not p2]
        p2.nbpoints = [p for p in p2.nbpoints if p is not p1]
        for p in self.points:
            p.nbedges = [e for e in p.nbedges if e is not self]
        mesh.edges = [e for e in mesh.edges if e is not self]
        if delobject:
            del self
        return
            

    def calclength(self):
        # Calculates the length of the edge
        return pythagoras( self.points[0], self.points[1] )

    def calcangle(self):
        # Calculates the absolute angle that the edge makes with the horizontal line
        (dx,dy) = np.subtract( self.points[0].coordinates, self.points[1].coordinates )
        if dx == 0:
            return np.pi/2
        else:
            return np.arctan(dy/dx)

    def findzero( self, zfun, mesh ):
        # Finds out whether the zero-levelset line crosses the edge and sets up level-set objects and properties
        (p1, p2) = self.points
        fp1      = zfun( p1.coordinates )
        p1.zdist = fp1
        fp2      = zfun( p2.coordinates )
        p2.zdist = fp2
        if fp1 * fp2 < 0 and abs(fp1) > self.hstep/10 and abs(fp2) > self.hstep/10:
            self.iszeroedge         = True
            mesh.zeroedges.append(self)
            (x,y)                   = zzerointerp( p1.coordinates, p2.coordinates, zfun )
            zp                      = zeropoint( x, y, mesh )
            self.zeropoint          = zp
            zp.onedge               = self
            if self.isboundary:
                zp.isboundary = True
            self.islevelset = False
            return
        else:
            self.iszeroedge = False
            lset = True
            for i,fp in enumerate( (fp1,fp2) ):
                if not abs(fp) > self.hstep/10:
                    p = self.points[i]
                    if not p.iszeropoint:    # (have to be two seperate if statements)
                        p.iszeropoint   = True
                        ( px, py )      = p.coordinates
                        zp              = zeropoint( px, py, mesh )
                        zp.isgrid       = True
                        zp.gridpoint    = p
                        p.zeropoint     = zp
                        mesh.zpointongrid.append( self.points[i] )
                        if p.isboundary:
                            zp.isboundary = True
                else:
                    self.points[i].iszeropoint  = False
                    lset                        = False
            if lset:
                self.islevelset     = True
                zp1                 = p1.zeropoint
                zp2                 = p2.zeropoint
                lsetedge            = levelsetedge( zp1, zp2, mesh )
                if self.isboundary:
                    lsetedge.isboundary = True
                self.levelsetedge   = lsetedge
                p1.levelsetedges.append( lsetedge )
                p2.levelsetedges.append( lsetedge )
            else:
                self.islevelset = False
            return
        
    def otherp(self,selfp):
        # Finds the other point of an edge other than the given point
        for p in self.points:
            if p is selfp:
                None
            else:
                return p

    def findmid(self):
        # Finds the midpoint of the edge
        return np.array( np.add( self.points[0].coordinates, self.points[1].coordinates ) ) / 2

    def split(self,mesh,delobject = True):
        # Creates a new point on the zeropoint and two new lines and returns the created elements.
        p1 = self.points[0]                             
        p2 = self.points[1]
        zp = self.zeropoint.coordinates
        self.zeropoint.delete(mesh)

        np  = point( zp[0], zp[1], mesh )
        ne1 = edge( np, p1, mesh )
        ne2 = edge( np, p2, mesh )

        self.delete(mesh,delobject = delobject)
        
        return [np,ne1,ne2]

class levelsetedge(edge):
    # Creates a seperate class for the level-set edges
    
    def __init__( self, p1, p2, mesh ):
        self.pprojected = []
        self.nbedges    = []
        self.points     = ( p1, p2 )
        for p in self.points:
            p.levelsetedges.append( self )
        self.length     = self.calclength()
        self.absangle   = self.calcangle()
        self.hstep      = mesh.hstep
        mesh.levelsetedges.append(self)
        self.zmid       = self.findmid()
        self.isboundary = False

    def delete(self,mesh,delobject = True):
        # Deletes the level-set edge while also cleaning up the surrounding mesh properties
        (p1,p2) = self.points
        p1.nbpoints = [p for p in p1.nbpoints if p is not p2]
        p2.nbpoints = [p for p in p2.nbpoints if p is not p1]
        for p in self.points:
            p.nbedges = [e for e in p.nbedges if e is not self]
        mesh.levelsetedges = [e for e in mesh.levelsetedges if e is not self]
        if delobject:
            del self
        return



##########################################################################
##### Element class ######################################################
##########################################################################

class element:
    # Builds an element class primarily consisting of three points and three lines
    
    def __init__(self, p1, p2, p3, mesh, edges = None):
        # Creates the element and finds the nbedges and nbelements, and sets up the basic info
        self.points     = ( p1, p2, p3 )
        self.zerotype   = None
        self.isboundary = False
        self.isflipped  = False
        self.iscut      = False
        self.pprojected = []
        self.hstep      = mesh.hstep
        for p in self.points:
            p.nbelements.append(self)

        if not edges == None:
            self.edges = edges
        else:
            self.edges  = []
            for e in p1.nbedges:
                if e == [p1,p2] or e == [p1,p3]:
                    self.edges.append(e)
                    e.nbelements.append(self)
            for e in p2.nbedges:
                if e == [p2,p3]:
                    self.edges.append(e)
                    e.nbelements.append(self)
        self.edges = tuple(self.edges)
        self.nbelements = []
        for e in self.edges:
            for el in e.nbelements:
                if el is not self:
                    self.nbelements.append(el)
        
        if len(self.edges)<3:
            print('Error not enough edges found in init method of element')
            for p in self.points:
                print(p)
            for e in self.edges:
                print(e)

        self.size     = self.calcsize()
        self.skewness = self.calcskew()
        mesh.elements.append(self)


    def __repr__(self):
        # Adds a string representation showing some properties of the element
        (p1,p2,p3) = self.points
        string  = 'Element object with \n'
        string += 'The coordinates of point 1 are: ' + str(p1.coordinates) + '\n'
        string += '        ,,         point 2 are: ' + str(p2.coordinates) + '\n'
        string += '        ,,         point 3 are: ' + str(p3.coordinates) + '\n'
        n = 1
        if len(self.zpointonedge) > 0:
            substring =     'The coordinates of zpoint ' + str(n) + ' on an edge are: '
            for zp in self.zpointonedge:
                string += substring + str(zp) + '\n'
                n += 1
                substring = '         ,,        zpoint ' + str(n) + ' on an edge are: '
        if len(self.zpoints) > 0:
            substring =     'The coordinates of zpoint ' + str(n) + ' on an edge are: '
            for zp in self.zpoints:
                string += substring + str(zp.coordinates)
                n += 1
                substring = '         ,,        zpoint ' + str(n) + ' on an edge are: '
        string += 'Points:\n'
        for p in self.points:
            string += repr(p)
        string += '\n Edges:\n'
        for e in self.edges:
            string += repr(e)
        return string

    def __eq__(self,other):
        # Finds out if two elements are the same object
        if isinstance(other,element):
            if self is other:
                return True
            else:
                return False

    def delete(self,mesh, delobject = True):
        # Deletes an element and cleans up the surrounding mesh properties in the process
        for el in self.nbelements:
            el.nbelements = [el for el in el.nbelements if el is not self]
        for p in self.points:
            p.nbelements = [el for el in p.nbelements if el is not self]
        for e in self.edges:
            e.nbelements = [el for el in e.nbelements if el is not self]
        mesh.elements = [el for el in mesh.elements if el is not self]
        if delobject:
            del self
        return

    def findmid(self):
        # Finds the centroid of the three corner points of the element
        (p1,p2,p3) = self.points
        (x1,y1) = p1.coordinates
        (x2,y2) = p2.coordinates
        (x3,y3) = p3.coordinates
        return ( (x1+x2+x3)/3, (y1+y2+y3)/3 )

    def calcskew(self):
        # Calculates the skewness of the element, only depends on minimal angle for a two dimensional simplex (triangle)
        d1 = self.edges[0].length
        d2 = self.edges[1].length
        d3 = self.edges[2].length
        for d in (d1,d2,d3):
            if d < 0.001*self.hstep:
                return 1
        a1 = np.arccos( (d2**2+d3**2-d1**2) / (2*d2*d3) )
        a2 = np.arccos( (d1**2+d3**2-d2**2) / (2*d1*d3) )
        a3 = np.pi - a1 - a2
        self.angles = (a1,a2,a3)
        amin = min(a1,a2,a3)
        return (np.pi/3 - amin) / (np.pi/3)

    def calcsize(self):
        # Calculates the size of the element using Heron's law
        d1 = self.edges[0].length
        d2 = self.edges[1].length
        d3 = self.edges[2].length
        p  = (d1+d2+d3)/2
        return np.sqrt( abs(p*(p-d1)*(p-d2)*(p-d3)) )   # the absolute value is to suppress faulty negative values
                                                        # ...for very skewed triangles because of rounding errors
    def otherp(self,p1,p2):
        # Finds the other point of the element other than the two specified points
        for p in self.points:
            if not ( (p == p1) or (p == p2) ):
                return p
            
    def othere(self,e1,e2):
        # Finds the other edge of the element other than the two specified edges
        for e in self.edges:
            if not ( (e == e1) or (e == e2) ):
                return e

    def findzero( self, zfun, mesh ):
        # Sets up the zero properties of the element and finds and creates level-set edges and zeropoints when necessary
        self.zedges        = []                 # Edges crossing the zero-levelset line
        self.zpointonedge  = []                 # Interpolated zero-levelset point on edge
        self.levelsetlines = []
        
        nze               = 0                   # Look for crossings of zero-levelset line with edges
        for e in self.edges:
            if e.iszeroedge == True:
                self.zedges.append(e)
                nze               += 1
                self.zpointonedge.append( e.zeropoint )

        nzp          = 0                        # Look for element points on the zero-levelset line
        self.zpoints = []
        for p in self.points:
            if p.iszeropoint == True:
                self.zpoints.append( p.zeropoint )
                nzp += 1

        if nze == 0 and nzp == 0:               # Return false when no zero-levelset line incident
            self.iszeroelement = False
            self.zerotype      = 0
            return

        elif nzp == 0:
            if nze == 2:                        # ... true when line crosses two edges
                self.iszeroelement = True
                self.zerotype      = 1
                ( zp1, zp2 )       = self.zpointonedge
                lsetedge           = levelsetedge( zp1, zp2, mesh )
                self.levelsetlines.append( lsetedge )
                return

        elif nzp == 1:
            if nze == 1:                        # ... true when line crosses an edge and a point
                self.iszeroelement = True
                self.zerotype      = 2
                ( zp1, zp2 )       = ( self.zpointonedge[0], self.zpoints[0] )
                lsetedge            = levelsetedge( zp1, zp2, mesh )
                self.levelsetlines.append( lsetedge )
                return 
            elif nze == 0:                      # ... false when line incident on just one point
                self.iszeroelement = False
                self.zerotype      = 5
                return

        elif nzp == 2:
            if nze == 0:                        # ... false when line incident on just two points
                self.iszeroelement = False
                self.zerotype      = 4
                return

        elif nzp == 3:                          # ... true when line incident on all three points
            self.iszeroelement = True
            self.zerotype      = 3
            self.midpoint       = self.findmid()
            (p1,p2,p3)          = self.points
            ( mpx, mpy )        = self.midpoint
            mp                  = zeropoint( mpx, mpy, mesh )
            mp.onelement        = self
            lsedge1             = levelsetedge( p1, mp, mesh )
            lsedge2             = levelsetedge( p2, mp, mesh )
            lsedge3             = levelsetedge( p3, mp, mesh )
            
            self.levelsetlines  = [ lsedge1, lsedge2, lsedge3 ]

            lsedge1.nbedges.append( lsedge2 )
            lsedge1.nbedges.append( lsedge3 )
            lsedge2.nbedges.append( lsedge1 )
            lsedge2.nbedges.append( lsedge3 )
            lsedge3.nbedges.append( lsedge1 )
            lsedge3.nbedges.append( lsedge2 )
            return
        else:
            raise NameError('in findzero of element, not nzp=0,1 or 2')


    def cut( self, mesh, flip = False ):
        # Applied the cut method to the element based on zerotype and zero properties
        if flip == True and self.zerotype == 2:
            ze = self.zedges[0]
            if len(ze.nbelements) == 2:                    # Checks for necessary conditions for flip
                if ze.nbelements[0] == self:
                    otherel = ze.nbelements[1]
                else:
                    otherel = ze.nbelements[0]

                if otherel.zerotype == 2:                  # Finds out which points are which
                    self.delete(mesh, delobject = False)
                    otherel.delete(mesh, delobject = False)
                    ze.delete(mesh, delobject = False)
                    
                    (pz1,pz2) = ze.points
                    opself    = self.otherp( pz1, pz2 )
                    opother   = otherel.otherp( pz1, pz2 )

                    for e in self.edges:                   # Finds out which edges are which for the creation of the new elements
                        if e == ( opself, pz1 ):
                            se1 = e
                            se2 = self.othere( se1, ze )
                            break
                        elif e == ( opself, pz2 ):
                            se2 = e
                            se1 = self.othere( se2, ze )
                            break

                    for e in otherel.edges:
                        if e == ( opother, pz1 ):
                            oe1 = e
                            oe2 = otherel.othere( oe1, ze )
                            break
                        elif e == ( opother, pz2 ):
                            oe2 = e
                            oe1 = otherel.othere( oe2, ze )
                            break

                    ne   = edge( opself, opother, mesh )       # Finally adds the new objects with the appropriate properties
                    ne.islevelset = True
                    nel1 = element( opself, pz1, opother, mesh, edges = (ne,se1,oe1) )
                    nel1.zerotype = 4
                    nel1.isflipped = True
                    nel2 = element( opself, pz2, opother, mesh, edges = (ne,se2,oe2) )
                    nel2.zerotype = 4
                    nel2.isflipped = True

                    del ze
                    del self
                    return otherel
                
        if self.zerotype == 2:                  # the no flip method for zeroedge and zeropoint

            ze = self.zedges[0]
            newcomponents = ze.split(mesh,delobject = False)
            
            self.delete(mesh,delobject = False)
            
            [np,ne1,ne2] = newcomponents
            (p1,p2) = ze.points
            op = self.otherp(p1,p2)
            otheredges = copy(self.edges)
            for e in otheredges:
                if e is ze:
                    otheredges = [ edge for edge in otheredges if edge is not e ]
                    break
            if otheredges[0] == (op,p1):
                oe1 = otheredges[0]
                oe2 = otheredges[1]
            else:
                oe1 = otheredges[1]
                oe2 = otheredges[0]

            ne            = edge( op, np, mesh )
            ne.islevelset = True
            nel1          = element( op, np, p1, mesh, edges = (ne,ne1,oe1) )
            nel1.zerotype = 4
            nel1.iscut    = True
            nel2          = element( op, np, p2, mesh, edges = (ne,ne2,oe2) )
            nel2.zerotype = 4
            nel2.iscut    = True

            del ze
            del self
            return

        elif self.zerotype == 1:                # no flip method for two zeroedges

            newcomponents = []                     
            splitedges = []
            
            self.delete(mesh,delobject = False)
            
            for e in self.zedges:
                newcomponents.append( e.split(mesh,delobject = False) )
                splitedges.append(e)
            
            if splitedges[0].points[0] in splitedges[1].points:
                cp = splitedges[0].points[0]
            else:
                cp = splitedges[0].points[1]
            if pythagoras( cp, splitedges[0].zeropoint ) > pythagoras( cp, splitedges[1].zeropoint ):
                [sp,se1,se2] = newcomponents[0]
                [np,ne1,ne2] = newcomponents[1]
                op1 = self.otherp( splitedges[0].points[0], splitedges[0].points[1] )
                op2 = self.otherp( splitedges[1].points[0], splitedges[1].points[1] )
            else:
                [sp,se1,se2] = newcomponents[1]
                [np,ne1,ne2] = newcomponents[0]
                op1 = self.otherp( splitedges[1].points[0], splitedges[1].points[1] )
                op2 = self.otherp( splitedges[0].points[0], splitedges[0].points[1] )
            if cp in se2.points:
                dummy = se1
                se1 = se2
                se2 = dummy
            if cp in ne2.points:
                dummy = ne1
                ne1 = ne2
                ne2 = dummy
            for e in self.edges:
                if not cp in e.points:
                    re = e
            de1            = edge( op1, sp, mesh )
            de2            = edge( sp, np, mesh )
            de2.islevelset = True
            nel1           = element( sp, op1, op2, mesh, edges = (re,se2,de1) )
            nel1.zerotype  = 5
            nel1.iscut     = True
            nel2           = element( sp, op1, np, mesh, edges = (de1,de2,ne2) )
            nel2.zerotype  = 4
            nel2.iscut     = True
            nel3           = element( sp, np, cp, mesh, edges = (se1,ne1,de2) )
            nel3.zerotype  = 4
            nel3.iscut     = True

            for e in splitedges:
                del e
            del self
            return



##########################################################################
##### Mesh class #########################################################
##########################################################################

class mesh:
    # The mesh object holds all the points, edges and elements
    
    def __init__(self,hstep,threshold,
                     x = np.array([-1,1,1,-1]),
                     y = np.array([-1,-1,1,1])  ):
        # Sets up the mesh points, edges and elements
        h               = np.array(hstep)
        self.hstep      = hstep
        self.threshold  = threshold
        self.relaxpoints = False
            # check whether saved mesh exists
        file_str = "Grid/"+"h"+np.array2string(h)+"_x"+np.array2string(x)+"_y"+np.array2string(y)+".npz"
        try:
            npzfile = np.load(file_str)
            for key in npzfile.keys():
                exec("self.f" + key + " = npzfile[key]")
        except IOError:
            print("Mesh does not exist\nCreating the mesh")
            # create mesh, random numbers
            self.felements, self.fedges, self.fpoints = create_mesh(x,y,h)
            print("Saving the mesh")
            np.savez(file_str, elements=elements,edges=edges,points=points)
            print("Mesh saved")

        self.points = []
        for i,x in enumerate(self.fpoints[0]):
            point( x, self.fpoints[1][i], self )
            
        self.edges = []
        for e in self.fedges:
            newedge = edge( self.points[e[0]], self.points[e[1]], self )
            newedge.threshold = threshold

        self.elements = []
        for el in self.felements:
            p1 = self.points[el[0]]
            p2 = self.points[el[1]]
            p3 = self.points[el[2]]
            newelem = element( p1, p2, p3, self )
            newelem.threshold = threshold

        self.boundarypoints = []
        self.boundaryedges  = []
        for e in self.edges:
            if len(e.nbelements) == 1:
                e.isboundary           = True
                e.points[0].isboundary = True
                e.points[1].isboundary = True
                e.nbelements[0].isboundary = True
                self.boundaryedges.append(e)
                self.boundarypoints.append(e.points[0])
                self.boundarypoints.append(e.points[1])

        for p in self.boundarypoints:     # sets up freedom vectors that describe in which direction the point is still allowed to move
            dvects = []
            for e in p.nbedges:
                if e.isboundary:
                    ( p1, p2 )   = e.points
                    displacement = np.subtract( p1.coordinates, p2.coordinates )
                    dvects.append( displacement )
            if np.cross( dvects[0], dvects[1] ) == 0:
                p.freedomvect = dvects[0] / np.sqrt( dvects[0][0]**2 + dvects[0][1]**2 )
            else:
                p.iscorner = True

        self.zeropoints    = []
        self.zpointongrid  = []
        self.zeroedges     = []
        self.levelsetedges = []

    def initzeroedgesinfo( self, zfun ):
        # Calculates the zero info for the edges in the mesh
        self.zeroedges     = []
        for e in self.edges:
            e.findzero( zfun, self )
        return

    def updatezeroelementinfo( self, zfun ):
        # Calculates the zero info for the elements in the mesh
        self.zeropoints    = []
        self.levelsetedges = []
        for el in self.elements:
            el.findzero( zfun, self )
        return

    def initzeroinfo( self, zfun ):
        # Caculates the zero info for the edges and element in the mesh
        self.zeropoints     = []
        self.zeroedges      = []
        self.levelsetedges  = []
        for e in self.edges:
            e.findzero( zfun, self )
        for el in self.elements:
            el.findzero( zfun, self )
        

##########################################################################
##### Mesh adjustment methods ############################################
##########################################################################

    def simplegridtozeropoint( self, thresh, z = False ):
        # Moves points that are closer to the level-set line than some threshold thresh
        adjustments = []
        for ze in self.zeroedges:
            zp = ze.zeropoint
            for p in ze.points:
                pdist  = pythagoras( p, zp )
                if pdist < thresh:
                    if p.isboundary:
                        if p.iscorner:
                            None
                        else:
                            if zp.isboundary:
                                adjustments.append( [ p, zp, pdist ] )
                    else:
                        adjustments.append( [ p, zp, pdist ] )
        adjustments.sort(key=lambda x: id(x[0]), reverse=True)      # making use of the sorted adjustments
        if len(adjustments)>0:                                      # list to choose most favorable
            adjustments2 = [adjustments[0]]                         # adjustment based on distance
            for i in range(1,len(adjustments)):
                if not adjustments[i][0] is adjustments2[-1][0]:
                    adjustments2.append(adjustments[i])
                elif adjustments[i][2] < adjustments[i-1][2]:
                    adjustments2[-1] = adjustments[i]
        else:
            adjustments2 = []

        for adj in adjustments2:                                    # making the adjustment
            if z != False:
                adj[0].move( adj[1], self, zfun = z, zpmove = True )
                
            else:
                adj[0].move( adj[1], self, zpmove = True )
        return

    def newtriang(self):
        # Applies the cut method to each zero element in the mesh
        for el in self.elements[:]:
            if el.iszeroelement == True:
                el.cut(self)
        return

    def newtriangwithflip(self):
        # Applies the cut and flip method to each element in the mesh
        # The cut method with Flip = True decides whether flip is applied or returns None
        ellist = copy(self.elements)
        ellist = [ el for el in ellist if el.iszeroelement == True ]
        while ellist != []:
            element = ellist[0]
            otherel = element.cut(self, flip = True)
            if otherel == None:
                ellist = [ el for el in ellist if el is not element ]
            else:
                ellist = [ el for el in ellist if el is not element ]
                ellist = [ el for el in ellist if el is not otherel ]
                del otherel
        return

    def setuppdistproperties(self):
        # Calculates how many edges each point is away from a zeroedge and sorts the points in the array pdistarray
        print('---Setting up zero distance properties')
        for p in self.points:
            p.nbzerodist = math.inf
            p.isclosetozero = False
        pdistarray = [[],[],[],[]]
        for e in self.zeroedges:
            for p0 in e.points:
                if p0.nbzerodist > 0:
                    p0.nbzerodist = 0
                    p0.isclosetozero = True
                    pdistarray[0].append(p0)
                    for p1 in p0.nbpoints:
                        if p1.nbzerodist > 1:
                            p1.nbzerodist = 1
                            p1.isclosetozero = True
                            pdistarray[1].append(p1)
                            for p2 in p1.nbpoints:                                
                                if p2.nbzerodist > 2:
                                    p2.nbzerodist = 2
                                    p2.isclosetozero = True
                                    pdistarray[2].append(p2)
                                    for p3 in p2.nbpoints:
                                        if p3.nbzerodist > 3:
                                            p3.nbzerodist = 3
                                            p3.isclosetozero = True
                                            pdistarray[3].append(p3)
        for i in range(len(pdistarray)):
            pdistarray[i] = [ p for p in pdistarray[i] if p.nbzerodist == i ]  # filters out the points that ended up in other distance arrays in pdistarray
        return pdistarray

    def setupedgedistproperties(self):
        """ setupedgedistproperties(self):
        The distance properties of the edges are initialized  using a branching algorithm where some branches will be cut
        off according to how much further away from the current zeroedge they are than from the closest zeroedge. The
        edge distance measures are updated along the branches of neighbouredges of the zeroedges, at most four edges deep."""

        print('---Setting up zero distance properties for edges')
        h  = self.hstep                                     # h2 defines how much further away the edge is allowed to be from <h2 is not used>
        h2 = 2*h                                            # ...the current zeroedge in comparison with other zeroedges, 
        for e in self.edges:                                # ...otherwise the algorithm is not allowed to continue updating
            e.zerodist = math.inf                           # ...in the current branch of nbedges.
            e.zdist    = math.inf
        for e in self.edges:
            if e.islevelset or e.iszeroedge:
                if e.islevelset:                            # First finding the reference zeroedge properties
                    e.zerodist = 0                          # ...when e.islevelset
                    zeroangle  = e.absangle
                    zmid       = e.findmid()
                elif e.iszeroedge:                          # ...and when e.iszeroedge
                    if e.isboundary:                        # ...using just one levelsetedge when e.isboundary
                        zedge         = e.nbelements[0].levelsetlines[0]
                        zeroangle     = zedge.absangle
                        zmid          = zedge.zmid
                    else:                                   # ...otherwise using two levelsetedges
                        zedge1        = e.nbelements[0].levelsetlines[0]
                        zedge2        = e.nbelements[1].levelsetlines[0]
                        zmid1         = zedge1.zmid
                        zmid2         = zedge2.zmid
                        zangle1       = zedge1.absangle
                        zangle2       = zedge2.absangle
                        
                        dangle        = angle( zangle1, zangle2 )       # ...calculating the average angle of the two lines
                        if dangle == abs( zangle1 -  zangle2 ):
                            zeroangle = min( zangle1, zangle2 ) + dangle/2
                        elif dangle == abs( zangle1 - zangle2 + np.pi ) or dangle == abs( zangle1 - zangle2 - np.pi ):
                            zeroangle = abs( min( zangle1, zangle2 ) - dangle/2 )
                        zmid        = ( zmid1 + zmid2 ) / 2
                    e.zangle    = angle( e, zeroangle )
                    e.closeto   = zmid
                    e.zdist     = pythagoras( zmid, e.findmid() )
                    
                    e.zerodist = zerodistmeasure( e.zangle, e.zdist, h )
                for p0 in e.points:                            # Then applying the branching type algorithm
                    for e0 in p0.nbedges:
                        zdist  = pythagoras( zmid, e0.findmid() )
                #        if zdist < e0.zdist + h2:
                        if zdist < e0.zdist:
                            zangle = angle( e0, zeroangle )
                            dist   = zerodistmeasure( zangle, zdist, h )
                            if dist < e0.zerodist:
                                e0.zerodist = dist
                                e0.zangle   = zangle
                                e0.zdist    = zdist
                                e0.closeto  = zmid
                            p1 = e0.otherp(p0)
                            for e1 in p1.nbedges:
                                zdist = pythagoras( zmid, e1.findmid() )
                    #            if zdist < e1.zdist + h2:
                                if zdist < e1.zdist:
                                    zangle = angle( e1, zeroangle )
                                    dist   = zerodistmeasure( zangle, zdist, h )
                                    if dist < e1.zerodist:
                                        e1.zerodist = dist
                                        e1.zangle   = zangle
                                        e1.zdist    = zdist
                                        e1.closeto  = zmid
##                                    p2 = e1.otherp(p1)                             # Two deep seemed sufficient, else two more levels can be uncommented
##                                    for e2 in p2.nbedges:
##                                        zdist = pythagoras( zmid, e2.findmid() )
##                                        if zdist < e2.zdist + h2:
##                                            zangle = angle( e2, zeroangle )
##                                            dist   = zerodistmeasure( zangle, zdist, h )
##                                            if dist < e2.zerodist:
##                                                e2.zerodist = dist
##                                                e2.zangle   = zangle
##                                                e2.zdist    = zdist
##                                                e2.closeto  = zmid
##                                            p3 = e2.otherp(p2)
##                                            for e3 in p3.nbedges:
##                                                zdist  = pythagoras( zmid, e3.findmid() )
##                                                zangle = angle( e3, zeroangle )
##                                                dist   = zerodistmeasure( zangle, zdist, h )
##                                                if dist < e3.zerodist:
##                                                    e3.zerodist = dist
##                                                    e3.zdist    = zdist
##                                                    e3.zangle   = zangle
##                                                    e3.closeto  = zmid
        return

    def pointrelax(self,relaxthresh):
        # Uses a pointwise mesh relaxation method (does not work)
        pdistarray = self.setuppdistproperties()

        somethingmoved = True
        while somethingmoved:
            somethingmoved = False
            for p in pdistarray[1]:
                if p.relaxpositionv2(relaxthresh):
                    somethingmoved = True
            for p in pdistarray[2]:
                if p.relaxpositionv2(relaxthresh):
                    somethingmoved = True
            for p in pdistarray[3]:
                if p.relaxpositionv2(relaxthresh):
                    somethingmoved = True
        return

    def findshortestpath( self, boundaryplist, boundarytype ):
        # Manages the search for shortest paths for the different cases where the level-set curve crosses the boundary (boundarytype = 1)
        # or when the level-set curve is closed (boundarytype = 2)
        print('---Finding shortest paths')
        if boundarytype == 1:
            pathlist = self.shortestpath( boundaryplist )
            for path in pathlist:
                path[0].iszeropath = True
                path[-1].iszeropath = True
            return pathlist
        elif boundarytype == 2:
            blist = copy(boundaryplist)
            path1 = self.shortestpath( boundaryplist )[0]
            numtokeep = math.floor(len(path1)/4)
            pathtoerase = path1[numtokeep:-numtokeep]
            for p1 in pathtoerase:
                for e in p1.nbedges:
                    e.zerodist = 100000            # Typical edge weights should be much less than 100000 or else increase this value
                    p2 = e.otherp( p1 )
                    for e in p2.nbedges:
                        e.zerodist = 100000
            path2 = self.shortestpath( blist )[0]
            path2[0].iszeropath = True
            path2[-1].iszeropath = True
            return [ path1, path2 ]

    def shortestpath( self, boundaryplist ):
        # Finds the shortest path for the shortest path method using the boundaryplist provided
        zeropathlist    = []
        print('length boundaryplist'+ str(len(boundaryplist)))
        while boundaryplist != []:

            for p in self.points:                # initializing pathlength property
                p.pathdist = math.inf

            startpoints         = boundaryplist.pop(0)
            overallbestpathdist = math.inf
            endplist = False
            for startpoint in startpoints:
                startpoint.path      = [startpoint]
                startpoint.pathdist  = 0
                endp                 = False

                ptoupdate = [startpoint]
                bestpath  = False
                bestdist  = math.inf
                while ptoupdate != []:                  # applying shortest path algorithm until another point in boundaryplist is reached
                    
                    p  = ptoupdate.pop(0)
                    if p.pathdist < bestdist:
                        dist        = p.pathdist
                        path        = p.path
                        for e in p.nbedges:
                            if e.isboundary == False:
                                newdist = dist + e.zerodist
                                newp    = e.otherp(p)
                                if newdist < newp.pathdist and (not newp.iszeropath):
                                    newp.pathdist = newdist
                                    newp.path     = path + [newp]
                                    
                                    ptoupdate     = [point for point in ptoupdate if point is not newp]
                                    if newp.isclosetozero and not newp.isboundary:
                                        i = 0
                                        for point in ptoupdate:
                                            if point.pathdist < newdist:
                                                i += 1
                                        ptoupdate.insert( i, newp )
                                        
                                    if bestpath != False:
                                        if newp in endplist and newdist < bestdist:
                                            bestpath = newp.path
                                            bestdist = newp.pathdist
                                            endp     = newp

                                    elif bestpath == False:
                                        for plist in boundaryplist:
                                            if newp in plist:
                                                endplist = plist
                                                bestpath = newp.path
                                                bestdist = newp.pathdist
                                    

                if bestdist < overallbestpathdist:
                    overallbestpathdist = bestdist
                    overallbestpath     = bestpath
            boundaryplist.remove( endplist )
            zeropathlist.append( overallbestpath )
            for p in overallbestpath[1:-1]:
                p.iszeropath = True
                
        return zeropathlist

    def movetolevelset( self, zeropathlist ):
        # Moves the points on a shortest path to suitable points on zeroedges

        for path in zeropathlist:
            for p1 in path:
                if (not p1.iszeropoint) and (not p1.isboundary):  # ...search for closest zmid of zeroedge in closeto zeroedges
                    
                    candidates      = []
                    othercandidates = []
                    elementstocheck = [ el for el in p1.nbelements if (el.iszeroelement or el.zerotype == 4) ]
                    for elem1 in p1.nbelements:
                        for elem2 in elem1.nbelements:
                            if (elem2.iszeroelement or elem2.zerotype == 4) and (elem2 not in elementstocheck):
                                elementstocheck.append( elem2 )
                    if elementstocheck == []:
                        print('ERROR, no elements to check in movetolevelset')
                        p1.movedtoelem = False
                    else:
                        for el in elementstocheck:
                            if el.iszeroelement:
                                for levelsetline in el.levelsetlines:
                                    candidates.append( [ project( p1, levelsetline ), el, levelsetline ] )
                            elif el.zerotype == 4:
                                for e in el.edges:
                                    if e.islevelset:
                                        candidates.append( [ project( p1, e ), el, e ] )
                        choice             = candidates[0]
                        choicedist         = pythagoras( p1, choice[0] )
                        for candidate in candidates:
                            newdist        = pythagoras( p1, candidate[0] )
                            if newdist < choicedist:
                                choice     = candidate
                                choicedist = newdist
                        p1.move( choice[0], self )
                        p1.iszeropath     = True
                        p1.movedtoelem     = choice[1]
                        p1.movedtolsetline = choice[2]
                        choice[1].pprojected.append( p1 )
                        choice[2].pprojected.append( p1 )
                    
                elif p1.isboundary and (not p1.iszeropoint) and (not p1.iscorner):      # ...search for candidates in nbedges on the boundary
                    candidate = False
                    for e in p1.nbedges:
                        if e.isboundary and e.iszeroedge:
                            elem = e.nbelements[0]
                            for p2 in elem.levelsetlines[0].points:
                                if p2.isboundary:
                                    candidate = [ p2, elem ]
                                    
                    if candidate == False:                        # ...search even further in nbedges on the boundary
                        for e1 in p1.nbedges:
                            if e1.isboundary:
                                p2 = e1.otherp( p1 )
                                if not p2.iscorner:
                                    for e2 in p2.nbedges:
                                        if e2.isboundary:
                                            for p3 in e2.nbelements[0].levelsetlines[0].points:
                                                if p3.isboundary:
                                                    candidate = [ p3, e2.nbelements[0] ]
                    if candidate != False:
                        p1.move( candidate[0], self )
                        p1.iszeropath  = True
                        p1.movedtoelem = candidate[1]
                        candidate[1].pprojected.append( p1 )
                        for levelsetline in candidate[1].levelsetlines:
                            for p in levelsetline.points:
                                if p1 == p:
                                    levelsetline.pprojected.append( p1 )
                    else:
                        print('ERROR, no candidate to move for boundary point in path in movetolevelset')
                        p1.movedtoelem = False
                else:
                    p1.movedtoelem = False

        return

    def redistribute_projectedpoints( self, zeropathlist, zfun ):
        # Applies the redistribution algorithm to redistribute projected points
        # It moves each projected point to the average position of its two neighbouring points and then projects it again
        print('---Redistributing orthogonally projected points')

        for path in zeropathlist:
            if path[0].isboundary or path[-1].isboundary:
                for dummy in range(10):
                    edgestoconsider    = []
                    elementstoconsider = []
                    for p in path:
                        for e in p.nbedges:
                            if e not in edgestoconsider:
                                edgestoconsider.append( e )
                        for elem in p.nbelements:
                            if elem not in elementstoconsider:
                                elementstoconsider.append( elem )
                    for e in edgestoconsider:
                        e.findzero( zfun, self )
                    for elem in elementstoconsider:
                        elem.findzero( zfun, self )
                        
                    halffloor = math.floor( len(path)/2 ) - 1
                    halfceil  = math.ceil( len(path)/2 ) - 1
                    for i in range( 0, halffloor ):
                        p = path[ halffloor - i ]
                        if not p.isboundary:
                            point1 = path[ halffloor - i - 1 ]
                            point2 = path[ halffloor - i + 1 ]
                            p1 = np.array( point1.coordinates )
                            p2 = np.array( point2.coordinates )
                            newcoords = list( ( p1 + p2 ) / 2 )
                            p.move( newcoords, self )
                            self.movetolevelset( [[p]] )
                    for i in range( 0, halfceil ):
                        p = path[ halfceil + i ]
                        if not p.isboundary:
                            point1 = path[ halfceil + i - 1 ]
                            point2 = path[ halfceil + i + 1 ]
                            p1 = np.array( point1.coordinates )
                            p2 = np.array( point2.coordinates )
                            oldcoords = np.array( p.coordinates )
                            newcoords = list( ( p1 + p2 ) / 2 )
                            p.move( newcoords, self )
                            self.movetolevelset( [[p]] )
                for p0 in [ path[0], path[-1] ]:                 # redistribute along the boundaries around boundarypoint of path
                    if p0.isboundary and not p0.iscorner:
                        pointstoconsider1 = [p0]
                        pointstoconsider2 = [p0]
                        firstpoints       = []
                        for p in p0.nbpoints:
                            if p.isboundary and (not p.iscorner) and (not p.iszeropoint):
                                firstpoints.append( p )
                                p.isconsidered = True
                        if len(firstpoints) == 2:
                            pointstoconsider1.append( firstpoints[0] )
                            pointstoconsider2.append( firstpoints[1] )
                            toconsiderlists = ( pointstoconsider1, pointstoconsider2 )
                        elif len(firstpoints) == 1:
                            pointstoconsider1.append( firstpoints[0] )
                            toconsiderlists = ( pointstoconsider1 )
                        else:
                            toconsiderlists = ()
                        for pointstoconsider in toconsiderlists:
                            j = 0
                            cont = True
                            while cont and j < 5:
                                discontinue = 0
                                cont = False
                                for p in pointstoconsider[-1].nbpoints:
                                    if p.isboundary and (p is not pointstoconsider[-2]) and (not p.iscorner) and (not p.iszeropoint) and (not p.isconsidered):
                                        pointstoconsider.append( p )
                                        p.isconsidered = True
                                        cont = True
                                    elif p.isboundary and p.iscorner:
                                        pointstoconsider.append( p )
                                        p.isconsiderd = True
                                        cont = False
                                j += 1

                        for pointstoconsider in toconsiderlists:
                            for dummy in range(10):
                                for i, p in enumerate( pointstoconsider[1:-1] ):
                                    if not p.iszeropoint:
                                        point1    = pointstoconsider[ i ]
                                        point2    = pointstoconsider[ i + 2 ]
                                        p1coords  = np.array( point1.coordinates )
                                        p2coords  = np.array( point2.coordinates )
                                        newcoords = list( ( p1coords + p2coords ) / 2 )
                                        if (p.ismoved == True) or (pythagoras( newcoords, p ) > 0.1*p.hstep):
                                            p.move( newcoords, self )
                            for p in pointstoconsider:
                                p.isconsidered = False
                        

            else:
                for i in range( 1, len(path)-1 ):
                    p = path[i]
                    if not p.isboundary and not p.iszeropoint:
                        point1 = path[ i - 1 ]
                        point2 = path[ ( i + 1 ) ]
                        p1 = np.array( point1.coordinates )
                        p2 = np.array( point2.coordinates )
                        newcoords = list( ( p1 + p2 ) / 2 )
                        p.move( newcoords, self )
                        self.movetolevelset( [[p]] )

    def setupboundaryplist(self):
        # Sets up the boundaryplist to be used in the shortest path algorithm
        # It differentiates between the two cases of a closed or boundary-crossing level-set curve
        pdistarray = self.setuppdistproperties()
        allp       = sum( pdistarray, [] )
        self.setupedgedistproperties()

        boundaryplist = []                          # building list of start and stop points
        for p1 in self.points:                      # ...for shortest path algorithm
            plist        = []
            
            if p1.iszeropoint and p1.isboundary:
                zp = p1.zeropoint
                n = 0
                m = 0
                for e in zp.levelsetedges:
                    if e.isboundary:
                        n += 1
                    else:
                        m += 1
                if n == 0:                          # when no levelsetedges are on the boundary,
                    plist.append(p1)                # ...add also the surrounding boundarypoints
                    p1.ischosenboundary = True
                    if not p1.iscorner:
                        for e in p1.nbedges:
                            if not e.nbelements[0].iszeroelement:
                                p2 = e.otherp(p1)
                                if not p2.iscorner:
                                    plist.append(p2)
                                    p2.ischosenboundary = True
                    boundaryplist.append(plist)
                elif m > 0:                         # when levelsetline is incident with boundary
                    for i in range(m):              # ...but not counting the levelsetedges on boundary
                        boundaryplist.append( [ p1 ] )
            elif len( p1.levelsetedges ) > 2:       # for intersections in the internal part of the grid
                for i in range( len( p1.levelsetedges ) ):
                    boundaryplist.append( [ p1 ] )
        for p1 in self.zeropoints:                  # checking intersections of levelsetlines in the internal of the grid
            if len( p1.levelsetedges ) > 2 and not p1.isboundary:       # ...and moving and adding appropriately
                if not p1.isgrid:
                    if p1.onedge != False:
                        e = p1.onedge
                        if pythagoras( e.points[0], p1 ) < pythagoras( e.points[1], p1 ):
                            p2 = e.points[0]
                        else:
                            p2 = e.points[1]
                        p2.move( p1, self )
                        p2.iszeropoint  = True
                        p2.zeropoint    = p1
                        for i in range( len( p1.levelsetedges ) ):
                            boundaryplist.append( [ p2 ] )
                    elif p1.onelement != False:
                        elem = p1.onelement
                        dist = math.inf
                        for p2 in elem.points:
                            newdist = pythagoras( p2, p1 )
                            if newdist < dist:
                                dist   = newdist
                                choice = p2
                        p2.move( p1, self )
                        p2.iszeropoint  = True
                        p2.zeropoint    = p1
                        for i in range( len( p1.levelsetedges ) ):
                            boundaryplist.append( [ p2 ] )
                    p2.isboundary       = True
                    p2.ischosenboundary = True
                else:
                    gp                  = p1.gridpoint
                    gp.isboundary       = True
                    gp.ischosenboundary = True
                    for i in range( len( p1.levelsetedges ) ):
                        boundaryplist.append( [ gp ] )

        for e1 in self.edges:       # filtering suitable points when levelsetedge has small angle of incidence with boundary
            plist = []
            if e1.iszeroedge and e1.isboundary and not e1.islevelset:
                (p1,p2)           = e1.points
                levelsetline      = e1.nbelements[0].levelsetlines[0]
                levelsetangle     = levelsetline.absangle
                dangle            = angle( levelsetangle, e1.absangle )
                if dangle < np.pi/4:
                    if p1.zdist < p2.zdist and (not p1.iscorner):
                        plist.append( p1 )
                        p1.ischosenboundary = True
                        for p3 in p1.nbpoints:
                            if p3.isboundary and (p3 is not p2) and (not p3.iscorner):
                                plist.append( p3 )
                                p3.ischosenboundary = True
                    elif (not p2.iscorner):
                        plist.append( p2 )
                        p2.ischosenboundary = True
                        for p3 in p2.nbpoints:
                            if p3.isboundary and (p3 is not p1) and (not p3.iscorner):
                                plist.append( p3 )
                                p3.ischosenboundary = True
                    else:
                        print('Error -- no suitable point found to fit to levelsetline due to freedom of movement restrictions on corner points.')
                    boundaryplist.append(plist)
                elif pythagoras( p1, e1.zeropoint ) < pythagoras( p2, e1.zeropoint ) and not p1.iscorner:
                    plist.append( p1 )
                    p1.ischosenboundary = True
                    for p3 in p1.nbpoints:
                        if p3.isboundary and not p3.iscorner:
                            plist.append( p3 )
                            p3.ischosenboundary = True
                    boundaryplist.append(plist)
                elif not p2.iscorner:
                    plist.append( p2 )
                    p2.ischosenboundary = True
                    for p3 in p2.nbpoints:
                        if p3.isboundary and not p3.iscorner:
                            plist.append( p3 )
                            p3.ischosenboundary = True
                    boundaryplist.append(plist)

        if boundaryplist == []:               # building list when zero level-set doesnt cross the boundary of the grid
            boundarytype = 2
            np0 = False
            if self.zpointongrid != []:
                np0     = self.zpointongrid[0]
                newdist = 0
                for p in self.zpointongrid:
                    dist = pythagoras( np0, p )
                    if dist > newdist:
                        np1 = p
                        newdist = dist
                newdist = 0
                for p in self.zpointongrid:
                    dist = pythagoras( np1, p )
                    if dist > newdist:
                        np2 = p
                        newdist = dist
                boundaryplist.append( [ np1 ] )
                np1.ischosenboundary    = True
                if newdist > 0:
                    boundaryplist.append( [ np2 ] )
                    np2.ischosenboundary    = True

            if len(boundaryplist) == 0:        # ...and when there are not enoungh zpointongrids when not crossing boundary
                zdist = math.inf
                for e in self.zeroedges:
                    for p in e.points:
                        dist = abs( p.zdist )
                        if dist < zdist:
                            zdist = dist
                            np0   = p
                boundaryplist.append( [ np0 ] )
                np0.ischosenboundary = True
            if len(boundaryplist) == 1:
                np0     = boundaryplist[0][0]
                newedge = self.zeroedges[0]
                newdist = pythagoras( np0, newedge.zmid )
                for e in self.zeroedges:
                    dist = pythagoras( np0, e.zmid )
                    if dist > newdist:
                        newedge = e
                        newdist = dist
                ( np1, np2 ) = e.points
                if np1.zdist > np2.zdist:
                    boundaryplist.append( [ np1 ] )
                    np1.ischosenboundary = True
                else:
                    boundaryplist.append( [ np2 ] )
                    np2.ischosenboundary = True
        else:
            boundarytype = 1

        self.boundaryplist = boundaryplist

        return boundarytype, boundaryplist

    def finallevelsetinfoupdate( self ):
        # Updates the levelset info of the edges one last time without altering structure (not adding new level-set lines)
        for e in self.edges:
            if e.points[0].iszeropath and e.points[1].iszeropath:
                e.islevelset = True

    def findrelaxpoints( self ):
        # Lists all the points that are allowed to be used in the mesh relaxation method (so no boudary or zeropath points)
        plist       = []
        i = 0
        for p in self.points:
            if not p.isboundary and not p.iszeropath:
                plist.append( p )
                p.isrelaxpoint  = True
                p.relaxindex    = i
                i += 1
        self.relaxpoints = plist
        
        elengths = []
        elist = [e for e in self.edges if not e.isboundary]
        for e in elist:
            elengths.append( e.length )
        self.reflength = np.mean( elengths )
        
        return plist

    def setupmatrix( self ):
        # Sets up the matrix for the fixed-point iteration for the mesh relaxation (diverges)
        if self.relaxpoints == False:
            self.findrelaxpoints()
        plist = self.relaxpoints

        reflength = self.reflength

        elist = [e for e in self.edges if not e.isboundary]
        for e in elist:
            e.springconstant = ( e.length - reflength ) / ( e.length )

        dim = len( plist[0].coordinates )
        n = dim * len( plist )
        A = np.zeros( (n,n) )
        b = np.zeros( n )
        u = np.zeros( n )
        for i,p1 in enumerate(plist):
            for j in range( dim ):
                u[dim*i+j] = p1.coordinates[j]
                p1val = 0
                for e in p1.nbedges:
                    p2 = e.otherp(p1)
                    if p2.isrelaxpoint:
                        A[dim*i+j][dim*p2.relaxindex+j] = e.springconstant
                    else:
                        b[dim*i+j] -= e.springconstant * p2.coordinates[j]
                    p1val -= e.springconstant
                A[dim*i+j][dim*i+j] = p1val
        return A,b,u

    def fixedpoint( self ):
        # Uses fixed-point iteration to relax the mesh (diverges)
        n = 1
        while True:
            print( '---Applying fixed point iteration: ' + str(n) )
            n      += 1
            (A,b,u) = self.setupmatrix()
            Ainv    = np.linalg.inv( A )
            u       = list( Ainv.dot( b ) )

            dim     = len( self.points[0].coordinates )
            moves   = []
            while not u == []:
                move = []
                for dummy in range(dim):
                    move.append( u.pop(0) )
                moves.append( move )

            maxdist = 0
            for i,newcoords in enumerate( moves ):
                p           = self.relaxpoints[i]
                oldcoords   = p.coordinates
                p.move( newcoords, self )
                newdist   = pythagoras( p.coordinates, oldcoords )
                if newdist > maxdist:
                    maxdist = newdist

            print('----------Maxdist = ' + str(maxdist) )

            if maxdist < self.hstep * 0.05 or n > 1:
                break
            
    def eulerforward( self ):
        # Uses an Euler-forward inspired method to relax the mesh (seems to always converge)
        print('---Relaxing points with Euler forward')
        plist = self.findrelaxpoints()
        pdistarray = self.setuppdistproperties()

        cont = True
        while cont:
            cont    = False
            for plist in pdistarray:
                for p1 in [p for p in plist if not p.isboundary and not p.iszeropath]:
                    F = np.array([0,0])
                    n = len( p1.nbpoints )
                    for e in p1.nbedges:
                        p2      = e.otherp(p1)
                        d       = np.subtract( p1.coordinates, p2.coordinates )
                        force   = -d*(1 - self.reflength/e.length)
                        F       = np.add( F, force )
                    p1.force = F/( n * self.reflength )
            
            for plist in pdistarray:
                for p1 in [p for p in plist if not p.isboundary and not p.iszeropath]:
                    displace = self.reflength * 2 * np.arctan( p1.force ) / np.pi
                    displacedist = np.sqrt( displace[0]**2 + displace[1]**2 )
                    if displacedist > self.reflength/100:
                        p1.move( p1.coordinates + displace, self )
                        cont = True
                    elif p1.ismoved and displacedist > self.reflength/1000:
                        p1.move( p1.coordinates + displace, self )
                        cont = True
                
            
  
##########################################################################
##### Quality methods ####################################################
##########################################################################

    def skewness(self):
        # Calculates the mesh-wide skewness measures and returns an array with the measures and their names
        skewness = []
        for el in self.elements:
            skewness.append( el.skewness )
        return [ ['Average skewness',sum(skewness)/len(skewness)],
                 ['Maximum skewness',max(skewness)],
                 ['Standard deviation of skewnesses',np.std(skewness)] ]

    def sizes(self):
        # Calculates the mesh-wide size measures and returns an array with the measures and their names
        sizes = []
        for el in self.elements:
            sizes.append( el.size )
        return [ ['Minimum size',min(sizes)],
                 ['Maximum size',max(sizes)],
                 ['Standard deviation of sizes',np.std(sizes)] ]

##########################################################################
##### Plot methods #######################################################
##########################################################################

    def converttocoordlist(self):
        # Converts the point data to coordinates for plotting purposes
        xcoords = []
        ycoords = []
        for i,p in enumerate(self.points):
            xcoords.append( p.coordinates[0] )
            ycoords.append( p.coordinates[1] )
            p.index = i

        elemindiceslist = []
        for el in self.elements:
            indices = []
            for p in el.points:
                indices.append( p.index )
            elemindiceslist.append( indices )

        return ( [xcoords, ycoords], elemindiceslist )

    def plot(self, show = True):
        # Plots the mesh itself using the coordinate lists with pyplot.triplot
        (points, elements) = self.converttocoordlist()
        triang = tri.Triangulation(points[0], points[1], elements)

        pyplot.figure()            # plot triangulation
        pyplot.gca().set_aspect('equal')
        pyplot.triplot(triang,'k-',lw=1) 
        pyplot.xlabel('$x$', fontsize = 20)
        pyplot.ylabel('$y$', fontsize = 20)
        pyplot.xticks(fontsize = 15)
        pyplot.yticks(fontsize = 15)
        if show:
            pyplot.show()

    def plotz(self, zfun, directory):
        # Plots the level-set function as surfaceplot and as contourplot
        (points, elements) = self.converttocoordlist()
        triang = tri.Triangulation(points[0], points[1], elements)

        z = zfun([points[0],points[1]])
        z = np.array(z)
        
        fig = pyplot.figure(figsize=(10,10))            # plots countours of some function zfun
        pyplot.gca().set_aspect('equal')
        pyplot.tricontourf(triang, z)
        pyplot.colorbar()
        pyplot.tricontour(triang,z,colors='k',levels=[-0.5,0,0.5]) 
        pyplot.xlabel('$x$', fontsize = 20)
        pyplot.ylabel('$y$', fontsize = 20)
        pyplot.xticks(fontsize = 15)
        pyplot.yticks(fontsize = 15)
        
        string = directory + zfun.__name__ + '_contour.png'
        fig.savefig(string, dpi=fig.dpi, bbox_inches = "tight")

        fig = pyplot.figure(figsize=(10,10))      # plot surface
        ax = fig.gca(projection='3d')
        ax.plot_trisurf(triang,z, cmap=pyplot.cm.CMRmap,alpha=0.75) 
        ax.set_xlabel('$x$', fontsize = 20)
        ax.set_ylabel('$y$', fontsize = 20)
        ax.set_zlabel('$z$',fontsize = 20)
        
        string = directory + zfun.__name__ + '_3dsurf.png'
        fig.savefig(string, dpi=fig.dpi, bbox_inches = "tight")


    def plotzeroz(self, savestring, dist = False, mode = 1, plotchosenboundary = False):
        # Plots the mesh along with the zeroelements, edges and points
        # mode defines the plot type, mode 1 shows zeroedges in red, zeroelement in grey etc.
        # mode 2 shows the size mismatch and the skewness in red and green respectively
        # the dist option depicts the distance of each edge as its color, from red to blue
        # plotchosenboundary plots the boundaryplist sets that are used in the shortest path method

        fig = pyplot.figure(figsize=(10,10))
        pyplot.gca().set_aspect('equal') # plot triangulation
        pyplot.xlabel('$x$', fontsize = 20)
        pyplot.ylabel('$y$', fontsize = 20)
        pyplot.xticks(fontsize = 15)
        pyplot.yticks(fontsize = 15)

        
        for el in self.elements:     
            (p1,p2,p3) = el.points
            (x1,y1) = p1.coordinates
            (x2,y2) = p2.coordinates
            (x3,y3) = p3.coordinates
            x = [x1,x2,x3]
            y = [y1,y2,y3]
            if mode == 1:
                if el.iszeroelement:                    # fill the zeroelements light grey
                    pyplot.fill(x,y,'xkcd:light grey')
                elif el.isflipped:
                    pyplot.fill(x,y,'xkcd:light blue')
                elif el.iscut:
                    pyplot.fill(x,y,'xkcd:light yellow')
            elif mode == 2:
                s    = el.skewness
                size = el.size
                h    = el.hstep
                refsize = 4 / len(self.elements)
                t    = 1 - 1 / ( 100*abs(size-refsize)/h + 1 )**(1/4)
                pyplot.fill(x,y,color = ( 1-s, 1-t, max(1-t-s, 0) ) )

        if not mode == 2:
            for e in self.levelsetedges:     # plot the levelsetedges blue
                (p1,p2) = e.points
                (x1,y1) = p1.coordinates
                (x2,y2) = p2.coordinates
                x       = [x1,x2]
                y       = [y1,y2]
                pyplot.plot(x,y,'xkcd:bright blue',linestyle = 'dotted')

        for e in self.edges:            # plots the edges according to type
            (p1,p2) = e.points
            (x1,y1) = p1.coordinates
            (x2,y2) = p2.coordinates
            x       = [x1,x2]
            y       = [y1,y2]
            if e.islevelset:
                if mode == 1:
                    pyplot.plot(x,y,'xkcd:bright green')
                elif mode == 2:
                    pyplot.plot(x,y,'xkcd:bright blue')
            elif not dist and mode == 1:
                if e.iszeroedge:
                    pyplot.plot(x,y,'r-',lw=1.2)
                else:
                    pyplot.plot(x,y,'k-',lw=1)
            elif not dist and mode == 2:
                pyplot.plot(x,y,color = (238/255,233/255,233/255),lw=1)
            elif dist:
                t = 1/(1*e.zerodist + 1)
                pyplot.plot(x,y,color = (t, 0.8*(1-t), 0.8*(1-t)),lw=0.5+t)

        for p in self.points:    # plots the points according to type
            if p.ismoved:
                (x,y) = p.coordinates
                pyplot.plot(x,y,color='orange',marker='o',markersize=3)
            if p.iszeropath:
                (x,y) = p.coordinates
                pyplot.plot(x,y,color='green',marker='o',markersize=5)
        if plotchosenboundary:
            n = len(self.boundaryplist)
            cm = pyplot.get_cmap('gist_rainbow')
            for i,plist in enumerate(self.boundaryplist):
                for p in plist:
                    (x,y) = p.coordinates
                    pyplot.plot(x,y,color=cm(i/n),marker='o',markersize=10)
                

        fig.savefig(savestring, dpi=fig.dpi, bbox_inches = "tight")

            
if __name__ == "__main__":     # for testing purposes
    m = mesh(0.1,0.001)


