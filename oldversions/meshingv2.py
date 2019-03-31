# -*- coding: utf-8 -*-
"""
Created on 

@author: twortelboer, douden
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



#### SOME FUNCTIONS ############################################################################

def create_mesh(x,y,h):
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
    return ( 1 + angle ) * ( 1 + distance/h )

def zzerointerp(p1,p2,zfun):  # finds the zeropoints along zeroedges
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
    for p in e1.points:
        if p in e2.points:
            return p

def angle(a1,a2):
    if isinstance(a1,edge):
        a1 = a1.absangle
    if isinstance(a2,edge):
        a2 = a2.absangle
    dangle = a1 - a2
    return min( abs( dangle ), abs( dangle - np.pi ), abs( dangle + np.pi ) )

def project( point, line ):  # calculates the orthogonal projection of a point on a line, or returns the closest boundary point of the line
    (p1,p2) = line.points
    linevec = np.array( np.subtract( p1.coordinates, p2.coordinates ) )
    pvec    = np.array( np.subtract( point.coordinates, p2.coordinates ) )
    pproj   = ( np.dot( pvec, linevec ) / np.dot( linevec, linevec ) ) * linevec
    newp    = list( np.add( p2.coordinates, pproj ) )
    minx    = min( p1.coordinates[0], p2.coordinates[0] )
    maxx    = max( p1.coordinates[0], p2.coordinates[0] )
    miny    = min( p1.coordinates[1], p2.coordinates[1] )
    maxy    = max( p1.coordinates[1], p2.coordinates[1] )
    if minx < newp[0] < maxx and miny < newp[1] < maxy:
        return newp
    elif pythagoras( point, p1 ) < pythagoras( point, p2 ):
        return p1
    else:
        return p2


#### INITIALIZATION OF point class ################################################################

class point:                           # Builds a point class consisting of coordinates
    def __init__(self,x,y,mesh):
        self.coordinates      = [x,y]
        self.nbedges          = []
        self.nbelements       = []
        self.nbpoints         = []
        self.levelsetedges    = []
        self.iszeropoint      = False
        self.ismoved          = False
        self.isboundary       = False
        self.iszeropath       = False
        self.iscorner         = False
        self.ischosenboundary = False
        self.isconsidered      = False
        self.hstep            = mesh.hstep
        self.nbzerodist       = math.inf
        mesh.points.append(self)

    def __repr__(self):                 # Adds a string representation with some basic info
        #string  = 'Point object with \n'
        string  = 'Coordinates: x = ' + str(self.coordinates[0]) + '\n'
        string += '             y = ' + str(self.coordinates[1]) + '\n'
        #string += 'Number of incident edges = ' + str(len(self.nbedges)) + '\n'
        #string += 'Number of indicent elements = ' + str(len(self.nbedges))
        #string = 'Point'
        return string

    def __eq__(self,other):
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
        for e in self.nbedges:
            e.delete(mesh)
        mesh.points = [p for p in mesh.points if p is not self]
        del self
        return

    def move( self, p, mesh, zfun = False, zpmove = False, updateinfo = True ):     # Performs the necessary operations to move a point to new coordinates p

        if isinstance( p, point ):
            p = p.coordinates
        self.coordinates = p
        self.ismoved     = True

        if zpmove:                    # use this part if the levelsetlines info and zerotypes is important for the rest of the algorithm
            self.iszeropoint = True
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
            for e in self.nbedges:
                e.length = e.calclength()
                if zfun != False:
                    e.findzero( zfun, mesh )
            for el in self.nbelements:
                el.size     = el.calcsize()
                el.skewness = el.calcskew()
                if zfun != False:
                    el.findzero( zfun, mesh )
            return

    def relaxposition(self,relaxthresh):  # does not give good results
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

    def relaxpositionv2(self,relaxthresh,uthresh = 0.01):  # is too complicated, use the relaxpoints for the whole mesh at once instead
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
        self.hstep          = mesh.hstep
        mesh.zeropoints.append(self)

    def delete( self, mesh ):
        for e in self.nbedges:
            e.delete(mesh)
        mesh.zeropoints = [p for p in mesh.zeropoints if p is not self]
        del self
        return


#### INITIALIZATION OF edge class ################################################################

class edge:                        # Builds an edge class consisting of two points
    def __init__(self,p1,p2,mesh):
        self.points     = ( p1, p2 )
        for p in self.points:
            p.nbedges += (self,)
        p1.nbpoints.append(p2)
        p2.nbpoints.append(p1)
        self.nbelements = []
        self.pprojected = []
        self.length     = self.calclength()
        self.islevelset = False
        self.iszeroedge = False
        self.isboundary = False
        self.hstep      = mesh.hstep
        self.absangle   = self.calcangle()
        self.zmid       = self.findmid()
        mesh.edges.append(self)

    def __repr__(self):
        string  = 'Edge object with \n'
        string += 'Number of incident elements = ' + str(len(self.nbelements)) + '\n'
        return string

    def __eq__(self,other):              # Decides if some other line is equal based on point objects
        (s1,s2) = self.points
        if isinstance(other,edge):
            if self is other:
                return True
            else:
                return False
        else:
            (o1,o2) = other
            if ( (s1 == o1) and (s2 == o2) ) or ( (s1 == o2) and (s2 == o1) ):
                return True
            else:
                return False

    def delete(self,mesh,delobject = True):
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
            

    def calclength(self):                # Calculates the length of the line
        return pythagoras( self.points[0], self.points[1] )

    def calcangle(self):
        (dx,dy) = np.subtract( self.points[0].coordinates, self.points[1].coordinates )
        if dx == 0:
            return np.pi/2
        else:
            return np.arctan(dy/dx)

    def findzero( self, zfun, mesh ):             # Finds out whether the zero-levelset line crosses the edge
        (p1, p2) = self.points
        fp1      = zfun( p1.coordinates )
        p1.zdist = fp1
        fp2      = zfun( p2.coordinates )
        p2.zdist = fp2
        if fp1 * fp2 < 0 and abs(fp1) > self.hstep/10 and abs(fp2) > self.hstep/10:
            self.iszeroedge = True
            mesh.zeroedges.append(self)
            (x,y)           = zzerointerp( p1.coordinates, p2.coordinates, zfun )
            self.zeropoint  = zeropoint( x, y, mesh )
            self.zeropoint.onedge = self
            if self.isboundary:
                self.zeropoint.isboundary = True
            self.islevelset = False
            return
        else:
            self.iszeroedge = False
            n=0
            for i,fp in enumerate( (fp1,fp2) ):
                if not abs(fp) > self.hstep/10:
                    self.points[i].iszeropoint = True
                    mesh.zpointongrid.append( self.points[i] )
                    n += 1
                else:
                    self.points[i].iszeropoint = False
            if n == 2:
                self.islevelset = True
                for p in self.points:
                    p.levelsetedges.append( self )
            else:
                self.islevelset = False
            return
        
    def otherp(self,selfp):
        for p in self.points:
            if p is selfp:
                None
            else:
                return p

    def findmid(self):
        return np.array( np.add( self.points[0].coordinates, self.points[1].coordinates ) ) / 2

    def split(self,mesh,delobject = True):              # Creates a new point on the zeropoint and two new lines
        p1 = self.points[0]                             # ... and returns the created elements.
        p2 = self.points[1]
        zp = self.zeropoint.coordinates
        self.zeropoint.delete(mesh)

        np  = point( zp[0], zp[1], mesh )
        ne1 = edge( np, p1, mesh )
        ne2 = edge( np, p2, mesh )

        self.delete(mesh,delobject = delobject)
        
        return [np,ne1,ne2]

class levelsetedge(edge):
    def __init__( self, p1, p2, mesh ):
        self.pprojected = []
        self.nbedges    = []
        self.points     = ( p1, p2 )
        for p in self.points:
            p.levelsetedges.append( self )
        if isinstance( p1, zeropoint ) and isinstance( p2, zeropoint ):
            p1.nbpoints.append(p2)
            p2.nbpoints.append(p1)
        self.length     = self.calclength()
        self.absangle   = self.calcangle()
        self.hstep      = mesh.hstep
        mesh.levelsetedges.append(self)
        self.zmid       = self.findmid()

    def delete(self,mesh,delobject = True):
        (p1,p2) = self.points
        p1.nbpoints = [p for p in p1.nbpoints if p is not p2]
        p2.nbpoints = [p for p in p2.nbpoints if p is not p1]
        for p in self.points:
            p.nbedges = [e for e in p.nbedges if e is not self]
        mesh.levelsetedges = [e for e in mesh.levelsetedges if e is not self]
        if delobject:
            del self
        return



#### INITIALIZATION OF element class ################################################################

class element:                           # Builds an element class consisting of three points and three lines
    def __init__(self, p1, p2, p3, mesh, edges = None):
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
        if isinstance(other,element):
            if self is other:
                return True
            else:
                return False

    def delete(self,mesh, delobject = True):
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
        (p1,p2,p3) = self.points
        (x1,y1) = p1.coordinates
        (x2,y2) = p2.coordinates
        (x3,y3) = p3.coordinates
        return ( (x1+x2+x3)/3, (y1+y2+y3)/3 )

    def calcskew(self):                        # Calculates the skewness of the element
        d1 = self.edges[0].length
        d2 = self.edges[1].length
        d3 = self.edges[2].length
        a1 = np.arccos( (d2**2+d3**2-d1**2) / (2*d2*d3) )
        a2 = np.arccos( (d1**2+d3**2-d2**2) / (2*d1*d3) )
        a3 = np.pi - a1 - a2
        self.angles = (a1,a2,a3)
        amin = min(a1,a2,a3)
        return (np.pi/3 - amin) / (np.pi/3)

    def calcsize(self):                         # Calculates the size of the element
        d1 = self.edges[0].length
        d2 = self.edges[1].length
        d3 = self.edges[2].length
        p  = (d1+d2+d3)/2
        return np.sqrt( p*(p-d1)*(p-d2)*(p-d3) )

    def otherp(self,p1,p2):
        for p in self.points:
            if not ( (p == p1) or (p == p2) ):
                return p
            
    def othere(self,e1,e2):
        for e in self.edges:
            if not ( (e == e1) or (e == e2) ):
                return e

    def findzero( self, zfun, mesh ):
        self.zedges        = []                  # Edges crossing the zero-levelset line
        self.zpointonedge  = []                  # Interpolated zero-levelset point on edge
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
                self.zpoints.append( p )
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
                self.levelsetlines.append( levelsetedge( zp1, zp2, mesh ) )
                return

        elif nzp == 1:
            if nze == 1:                        # ... true when line crosses an edge and a point
                self.iszeroelement = True
                self.zerotype      = 2
                ( zp1, zp2 )       = ( self.zpointonedge[0], self.zpoints[0] )
                self.levelsetlines.append( levelsetedge( zp1, zp2, mesh ) )
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

#### INITIALIZATION OF mesh class ######################################################################

class mesh:

    def __init__(self,hstep,threshold,
                     x = np.array([-1,1,1,-1]),
                     y = np.array([-1,-1,1,1])  ):
        h               = np.array(hstep)
        self.hstep      = hstep
        self.threshold  = threshold
            # check whether saved mesh exists
        file_str = "h"+np.array2string(h)+"_x"+np.array2string(x)+"_y"+np.array2string(y)+".npz"
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

        for p in self.boundarypoints:
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
            
                

    def initzeroedgesinfo( self, zfun, mesh ):
        self.zeroedges     = []
        for e in self.edges:
            e.findzero( zfun, mesh )
        return

    def updatezeroelementinfo( self, zfun, mesh ):
        self.zeropoints    = []
        self.levelsetedges = []
        for el in self.elements:
            el.findzero( zfun, mesh )
        return


#### MESH ADJUSTMENT METHODS ##########################################################################

    def simplegridtozeropoint( self, thresh, mesh, z = False ):

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
                adj[0].move( adj[1], mesh, zfun = z, zpmove = True )
                
            else:
                adj[0].move( adj[1], mesh, zpmove = True )
        return

    def newtriang(self):
        for el in self.elements[:]:
            if el.iszeroelement == True:
                el.cut(self)
        return

    def newtriangwithflip(self):
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
        print('---Setting up zero distance properties')
        for p in self.points:
            p.nbzerodist = math.inf
            p.isclosetozero = False
        pdistarray = [[],[]]
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
##                                for e2 in p1.nbedges:
##                                    p2 = e2.otherp(p1)
##                                    if p2.nbzerodist > 2:
##                                        p2.nbzerodist = 2
##                                        p2.isclosetozero = True
##                                        pdistarray[2].append(p2)
##                                        for e3 in p2.nbedges:
##                                            p3 = e3.otherp(p2)
##                                            if p3.nbzerodist > 3:
##                                                p3.nbzerodist = 3
##                                                p3.isclosetozero = True
##                                                pdistarray[3].append(p3)
        for i in range(len(pdistarray)):
            pdistarray[i] = [ p for p in pdistarray[i] if p.nbzerodist == i ]  # filters out the points that ended up in other distance arrays in pdistarray
        return pdistarray

    def setupedgedistproperties(self):
        """ setupedgedistproperties(self):
        The distance properties of the edges are initialized  using a branching algorithm where some branches will be cut
        off according to how much further away from the current zeroedge they are than from the closest zeroedge. The
        edge distance measures are updated along the branches of neighbouredges of the zeroedges, at most four edges deep."""

        print('---Setting up zero distance properties for edges')
        h  = self.hstep                                     # h2 defines how much further away the edge is allowed to be from
        h2 = 2*h                                            # ...the current zeroedge in comparison with other zeroedges, 
        for e in self.edges:                                # ...otherwise the algorithm is not allowed to continue updating
            e.zerodist = math.inf                           # ...in the current branch of nbedges.
            e.zdist    = math.inf
        for e in self.edges:
            if e.islevelset or e.iszeroedge:
                if e.islevelset:                            # First finding the reference zeroedge properties
                    e.zerodist = 1                          # ...when e.islevelset
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
##                                    p2 = e1.otherp(p1)
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
        print('---Finding shortest paths')
        if boundarytype == 1:
            return self.shortestpath( boundaryplist )
        elif boundarytype == 2:
            blist = copy(boundaryplist)
            path1 = self.shortestpath( boundaryplist )[0]
            numtokeep = math.floor(len(path1)/4)
            pathtoerase = path1[numtokeep:-numtokeep]
            for p1 in pathtoerase:
                for e in p1.nbedges:
                    e.zerodist = 100000
                    p2 = e.otherp( p1 )
                    for e in p2.nbedges:
                        e.zerodist = 100000
            path2 = self.shortestpath( blist )[0]
            return [ path1, path2 ]

    def shortestpath( self, boundaryplist ):
        zeropathlist = []
        while boundaryplist != []:

            for p in self.points:                # initializing pathlength property
                p.pathdist = math.inf

            startpoints         = boundaryplist.pop(0)
            overallbestpathdist = math.inf
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
                                if newdist < newp.pathdist:
                                    newp.pathdist = newdist
                                    newp.path     = path + [newp]
                                    
                                    ptoupdate     = [point for point in ptoupdate if point is not newp]
                                    if newp.isclosetozero and not newp.isboundary:
                                        i = 0
                                        for point in ptoupdate:
                                            if point.pathdist < newdist:
                                                i += 1
                                        ptoupdate.insert( i, newp )
                                        
                                    if bestpath == False:
                                        for plist in boundaryplist:
                                            if newp in plist:
                                                endplist = plist
                                                bestpath = newp.path
                                                bestdist = newp.pathdist
                                    if bestpath != False:
                                        if newp in endplist and newdist < bestdist:
                                            bestpath = newp.path
                                            bestdist = newp.pathdist
                                            endp     = newp

                if bestdist < overallbestpathdist:
                    overallbestpathdist = bestdist
                    overallbestpath     = bestpath
            boundaryplist.remove( endplist )
            zeropathlist.append( overallbestpath )

        for path in zeropathlist:
            for p in path:
                p.iszeropath = True
            
        return zeropathlist

    def movetolevelset( self, zeropathlist ):                     # Moves the points on a shortest path to suitable points on zeroedges
        #print('---Moving to level set line')

        for path in zeropathlist:
            for p1 in path:
                if (not p1.iszeropoint) and (not p1.isboundary):  # ...search for closest zmid of zeroedge in closeto zeroedges
                    
                    candidates      = []
                    othercandidates = []
                    elementstocheck = [ el for el in p1.nbelements if (el.iszeroelement or el.zerotype == 4) ]
                    for elem1 in p1.nbelements:
                        for elem2 in elem1.nbelements:
                            if (elem2.iszeroelement or elem2.zerotype == 4) and elem2 not in elementstocheck:
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
                    if not p.isboundary:
                        point1 = path[ i - 1 ]
                        point2 = path[ ( i + 1 ) ]
                        p1 = np.array( point1.coordinates )
                        p2 = np.array( point2.coordinates )
                        newcoords = list( ( p1 + p2 ) / 2 )
                        p.move( newcoords, self )
                        self.movetolevelset( [[p]] )

    def setupboundaryplist(self):
        pdistarray = self.setuppdistproperties()
        allp       = sum( pdistarray, [] )
        self.setupedgedistproperties()

        boundaryplist = []                      # building list of start and stop points for shortest path algorithm
        for p1 in self.points:
            plist        = []
            p1.isupdated = False
            cont         = True
            for el in p1.nbelements:
                if el.zerotype == 4:
                    cont = False
            if p1.iszeropoint and p1.isboundary and cont:
                plist.append(p1)
                p1.ischosenboundary = True
                if not p1.iscorner:
                    for e in p1.nbedges:
                        if not e.nbelements[0].iszeroelement:
                            p2 = e.otherp(p1)
                            if not p2.iscorner:
                                plist.append(p2)
                                p2.ischosenboundary = True
                boundaryplist.append(plist)
            elif len( p1.levelsetedges ) > 2:
                for i in range( len( p1.levelsetedges ) - 2 ):
                    boundaryplist.append( [ p1 ] )
     #           print('check' + str( len(boundaryplist) ) )
        for p1 in self.zeropoints:
            if len( p1.levelsetedges ) > 2:
                if p1.onedge != False:
                    e = p1.onedge
                    if pythagoras( e.points[0], p1 ) < pythagoras( e.points[1], p1 ):
                        p2 = e.points[0]
                    else:
                        p2 = e.points[1]
                    p2.move( p1, self )
                    for i in range( len( p1.levelsetedges ) - 2 ):
                        boundaryplist.append( [ p2 ] )
                    p2.isboundary = True
   #                 print('check' + str( len(boundaryplist) ) )
                elif p1.onelement != False:
                    elem = p1.onelement
                    dist = math.inf
                    for p2 in elem.points:
                        newdist = pythagoras( p2, p1 )
                        if newdist < dist:
                            dist   = newdist
                            choice = p2
                    p2.move( p1, self )
                    for i in range( len( p1.levelsetedges ) - 2 ):
                        boundaryplist.append( [ p2 ] )
                    p2.isboundary = True
     #               print('check' + str( len(boundaryplist) ) )
                

        for e1 in self.edges:
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
                        print('Error -- not suitable point found to fit to levelsetline due to freedom of movement restrictions on corner points.')
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
                np1.ischosenboundary = True
                if newdist > 0:
                    boundaryplist.append( [ np2 ] )
                    np2.ischosenboundary = True

            if len(boundaryplist) == 0:        # ...and when there are not enoungh zpointongrids
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

#### QUALITY METHODS ###################################################################################

    def skewness(self):
        skewness = []
        for el in self.elements:
            skewness.append( el.skewness )
        #print('---- No skewness entries: ' + str(len(skewness)))
        return [ ['Average skewness',sum(skewness)/len(skewness)],
                 ['Maximum skewness',max(skewness)],
                 ['Standard deviation of skewnesses',np.std(skewness)] ]

    def sizes(self):
        sizes = []
        for el in self.elements:
            sizes.append( el.size )
        #print('---- No sizes entries: ' + str(len(sizes)))
        return [ ['Minimum size',min(sizes)],
                 ['Maximum size',max(sizes)],
                 ['Standard deviation of sizes',np.std(sizes)] ]


#### PLOT METHODS #######################################################################################

    def converttocoordlist(self):
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

    def plotz(self, zfun, directory):    # plots the level-set function as surface and as contours
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


    def plotzeroz(self, savestring, dist = False, mode = 1, plotchosenboundary = False):                 # plots the mesh along with the zeroelements, edges and points
##        (points, elements) = self.converttocoordlist()
##        triang = tri.Triangulation(points[0], points[1], elements) # ...found in findzero

        fig = pyplot.figure(figsize=(10,10))
        pyplot.gca().set_aspect('equal') # plot triangulation
        #pyplot.triplot(triang,'k-',lw=1)  
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
                t = 1 / (e.zerodist)
                pyplot.plot(x,y,color = (t, 0.8*(1-t), 0.8*(1-t)),lw=0.5+t)

        for e in self.levelsetedges:     # plot the levelsetedges blue
            (p1,p2) = e.points
            (x1,y1) = p1.coordinates
            (x2,y2) = p2.coordinates
            x       = [x1,x2]
            y       = [y1,y2]
            pyplot.plot(x,y,'xkcd:bright blue',linestyle = 'dotted')

        for p in self.zeropoints:        # plot the zeropoints black
            (x,y) = p.coordinates
            pyplot.plot(x,y,'ko',markersize = 1.5)
        for p in self.points:
            if p.ismoved:
                (x,y) = p.coordinates
                pyplot.plot(x,y,color='orange',marker='o',markersize=5)
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

            



