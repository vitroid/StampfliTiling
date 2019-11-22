#Put the initial vertex
#subdivide with dodecagon and register new vertices
from math import *
import numpy as np
import yaplotlib as yp
import pairlist as pl
import itertools as it
from FrankKasper import toWater

# dummy funcs
HSB=0

def nofill():
    pass

def stroke(*x):
    pass

def colormode(x):
    pass

def translate(x,y):
    pass

def oval(*a):
    pass

def nostroke(*a):
    pass

def fill(*a):
    pass

def line(*a):
    pass

def fontsize(*x):
    pass


def _mul(mat,vec):
    return ((mat[0][0]*vec[0]+mat[0][1]*vec[1]),
            (mat[1][0]*vec[0]+mat[1][1]*vec[1]))

def _add(vec1,vec2):
    return (vec1[0]+vec2[0],vec1[1]+vec2[1])

def _sub(vec1,vec2):
    v = [0] * len(vec1)
    for i in range(len(vec1)):
        v[i] = vec[1] - vec[2]
    return v

def wrap(c, box):
    s = c.shape[0]
    return c - np.floor( c/box[:s] + 0.5 )*box[:s]

def _sub_pbc(vec1,vec2):
    v = [0] * len(vec1)
    for i in range(len(vec1)):
        v[i] = vec1[i]-vec2[i]
    return wrap(v)


def drawpoly( coord ):
    for i in range(len(coord)):
        line(coord[i-1][0],coord[i-1][1],
             coord[i  ][0],coord[i  ][1])

def subdiv_triangle(v, edges):
    for i in range(3):
        dijx,dijy = edges[v[i]][v[i-1]]
        dikx,diky = edges[v[i]][v[i-2]]
        p = (coord[v[i]][0] + (dijx+dikx)*(0.5 - ratio / sqrt(3.0)),
             coord[v[i]][1] + (dijy+diky)*(0.5 - ratio / sqrt(3.0)))
        coord.append(np.array(p))

def subdiv_square(v, edges):
    for i in range(4):
        #relative diagonal vector of the square
        d = edges[v[i]][v[i-1]] + edges[v[i-1]][v[i-2]]
        p0 = d/sqrt(2)*ratio2
        p = coord[v[i]]+p0
        coord.append(np.array(p))
        p1 = matrix12@p0
        p = coord[v[i]]+p1
        coord.append(np.array(p))
        

def hexagon(center,edge,parity):
    r = edge * ratio * 2.0
    for i in range(6):
        coord.append(center+np.array([r*cos((i*2+parity)*pi/6),r*sin((i*2+parity)*pi/6)]))


def bond(coord,thres,box):
    edges = dict()
    for i in range(len(coord)):
        edges[i] = dict()
    for i in range(len(coord)):
        for j in range(i+1,len(coord)):
            # dx,dy = sub_pbc(coord[j],coord[i])
            d = wrap(coord[j]-coord[i], box)
            #if dx**2+dy**2 < thres**2 * 1.01:
            if d@d < thres**2 * 1.01:
                edges[i][j] = d # (dx,dy)
                edges[j][i] = -d # (-dx,-dy)
    return edges


def dictoflist_add(dol,key,value):
    if key not in dol:
        dol[key] = []
    dol[key].append(value)

def bond3d_fast(coord,thres,box):
    edges = dict()
    for i in range(len(coord)):
        edges[i] = dict()
    ndiv = [int(x) for x in (floor(box[0]/thres),floor(box[1]/thres),floor(box[2]/thres))]
    binw = box[0]/ndiv[0],box[1]/ndiv[1],box[2]/ndiv[2]
    resident = dict()
    for i in range(len(coord)):
        c = coord[i]
        ix,iy,iz = int(floor( c[0] / binw[0] )),int(floor( c[1] / binw[1] )),int(floor( c[2] / binw[2] ))
        if ix < 0:
            ix += ndiv[0]
        if iy < 0:
            iy += ndiv[1]
        if iz < 0:
            iz += ndiv[2]
        dictoflist_add(resident,(ix,iy,iz),i)
    for ix in range(ndiv[0]):
        for iy in range(ndiv[1]):
            for iz in range(ndiv[2]):
                for dx in range(-1,2):
                    for dy in range(-1,2):
                        for dz in range(-1,2):
                            jx = ix + dx
                            jy = iy + dy
                            jz = iz + dz
                            if jx < 0:
                                jx += ndiv[0]
                            if jy < 0:
                                jy += ndiv[1]
                            if jz < 0:
                                jz += ndiv[2]
                            jx %= ndiv[0]
                            jy %= ndiv[1]
                            jz %= ndiv[2]
                            #print jx,jy,jz,resident[(jx,jy,jz)]
                            for i in resident[(ix,iy,iz)]:
                                for j in resident[(jx,jy,jz)]:
                                    if i != j:
                                        xyz = wrap(coord[i]-coord[j], box)
                                        if xyz@xyz < thres**2:
                                            edges[j][i] = xyz
    return edges


def bond3d(coord,thres,box):
    edges = dict()
    for i in range(len(coord)):
        edges[i] = dict()
    for i in range(len(coord)):
        for j in range(len(coord)):
            if i != j:
                xyz = wrap(coord[i]-coord[j], box)
                if xyz@xyz < thres**2:
                    edges[j][i] = xyz
    return edges
            

SAMECOLOR = 1
DIFFCOLOR = -1


def setcolors(en,ec,e,col):
    if e in ec:
        return
    ec[e] = col
    #print e,col
    for i,j,parity in en[e]:
        setcolors(en,ec,(i,j),col * parity)


def findrings(edges):
    #look for rings
    sq = []
    tr = []
    edgeneighbor = dict()
    for i in edges:
        for j in edges[i]:
            for k in edges[i]:
                if i < j and i < k and j < k:
                    if k in edges[j]:
                        #must be compact
                        dx = edges[i][j][0] + edges[j][k][0] + edges[k][i][0]
                        dy = edges[i][j][1] + edges[j][k][1] + edges[k][i][1]
                        if abs(dx)<0.001 and abs(dy) < 0.001:
                            tr.append((i,j,k))
                            dictoflist_add(edgeneighbor,(i,j),(i,k,SAMECOLOR))
                            dictoflist_add(edgeneighbor,(i,j),(j,k,SAMECOLOR))
                            dictoflist_add(edgeneighbor,(i,k),(i,j,SAMECOLOR))
                            dictoflist_add(edgeneighbor,(i,k),(j,k,SAMECOLOR))
                            dictoflist_add(edgeneighbor,(j,k),(i,j,SAMECOLOR))
                            dictoflist_add(edgeneighbor,(j,k),(i,k,SAMECOLOR))
                    else:
                        for l in edges[j]:
                            if l in edges[k]:
                                if l > i:
                                    if l not in edges[i]:
                                        #must be compact
                                        dx = edges[i][j][0] + edges[j][l][0] + edges[l][k][0] + edges[k][i][0]
                                        dy = edges[i][j][1] + edges[j][l][1] + edges[l][k][1] + edges[k][i][1]
                                        if abs(dx)<0.001 and abs(dy) < 0.001:
                                            sq.append((i,j,l,k))
                                            kk,kl = k,l
                                            if k > l:
                                                kk,kl = l,k
                                            jj,jl = j,l
                                            if j > l:
                                                jj,jl = l,j
                                            dictoflist_add(edgeneighbor,(i,j),(i,k,DIFFCOLOR))
                                            dictoflist_add(edgeneighbor,(i,k),(i,j,DIFFCOLOR))
                                            dictoflist_add(edgeneighbor,(i,k),(kk,kl,DIFFCOLOR))
                                            dictoflist_add(edgeneighbor,(kk,kl),(i,k,DIFFCOLOR))
                                            dictoflist_add(edgeneighbor,(kk,kl),(jj,jl,DIFFCOLOR))
                                            dictoflist_add(edgeneighbor,(jj,jl),(kk,kl,DIFFCOLOR))
                                            dictoflist_add(edgeneighbor,(i,j),(jj,jl,DIFFCOLOR))
                                            dictoflist_add(edgeneighbor,(jj,jl),(i,j,DIFFCOLOR))
                                            dictoflist_add(edgeneighbor,(i,j),(kk,kl,SAMECOLOR))
                                            dictoflist_add(edgeneighbor,(i,k),(jj,jl,SAMECOLOR))
                                            dictoflist_add(edgeneighbor,(kk,kl),(i,j,SAMECOLOR))
                                            dictoflist_add(edgeneighbor,(jj,jl),(i,k,SAMECOLOR))
    return tr,sq,edgeneighbor
 

def inflate(coord,edgelen,edges,triangles,squares,depth):
    clen = len(coord)
    if depth == 0:
        for i in range(clen):
            hexagon(coord[i],edgelen,0) #ast number is the direction of hexagon
    else:
        for i in range(clen):
            hexagon(coord[i],edgelen,i%2) #ast number is the direction of hexagon
    for triangle in triangles:
        subdiv_triangle(triangle, edges)
    for square in squares:
        subdiv_square(square, edges)
    return coord


blue = 1
red  = -1


def onelayer(coord,edgelen,edges,tr,sq,edgeneighbor,box):
    """
    2D layer to 3D layer
    """
    atoms = []
    i = 0
    for c in coord:
        stroke(0)
        nofill()
        xy = wrap(c, box)
        x,y = xy
        oval( x-4,y-4,8,8 )
        nostroke()
        fill(0)
        oval( x-2,y-2,4,4 )
        #text("%s" % i,c[0],c[1] )
        atoms.append((x,y,box[2]/4))
        atoms.append((x,y,box[2]*3/4))
        i += 1
    #return
    edgecolors = dict()
    i = list(edges.keys())[0]
    j = list(edges[i].keys())[0]
    if i > j:
        i,j = j,i
    setcolors(edgeneighbor,edgecolors,(i,j),red)
    for i in edges:
        for j in edges[i]:
            if i < j:
                dx,dy = edges[i][j]
                #ij = (coord[i][0]+dx*0.5,
                #      coord[i][1]+dy*0.5)
                ij = coord[i] + edges[i][j]*0.5
                x,y = wrap(ij, box)
                if edgecolors[(i,j)] == red:
                    stroke(0)
                    nofill()
                    oval( x-3,y-3,6,6 )
                    atoms.append((x,y,0.0))
                    stroke(0,1,1)
                else:
                    nostroke()
                    fill(0)
                    oval( x-3,y-3,6,6 )
                    atoms.append((x,y,box[2]/2))
                    stroke(0.666,1,1)
                x,y = wrap(coord[i], box)
                line(x,y,x+dx,y+dy)
    for i,j,k in tr:
        #dijx,dijy = edges[i][j]
        #dikx,diky = edges[i][k]
        #center = (coord[i][0]+(dijx+dikx)/3,
        #          coord[i][1]+(dijy+diky)/3)
        center = coord[i] + (edges[i][j]+edges[i][k])/3
        x,y = wrap(center, box)
        if edgecolors[(i,j)] == red:
            nostroke()
            fill(0)
            z = box[2]/2
        else:            
            stroke(0)
            nofill()
            z = 0.0
        oval(x-3,y-3,6,6)
        atoms.append((x,y,z))
        #text("(%s,%s,%s)" % (i,j,k), center[0]+5,center[1]+5)
        #print (i,j,k),center
    for i,j,l,k in sq:
        p = 3.0/8.0
        q = 1.0/8.0
        #dijx,dijy = edges[i][j]
        #dikx,diky = edges[i][k]
        #dilx,dily = dijx+dikx, dijy+diky
        #ij = (coord[i][0] + dijx*p + dikx*q + dilx*q,
        #      coord[i][1] + dijy*p + diky*q + dily*q)
        #ik = (coord[i][0] + dijx*q + dikx*p + dilx*q,
        #      coord[i][1] + dijy*q + diky*p + dily*q)
        #jl = (coord[i][0] + dijx*p + dikx*q + dilx*p,
        #      coord[i][1] + dijy*p + diky*q + dily*p)
        #kl = (coord[i][0] + dijx*q + dikx*p + dilx*p,
        #      coord[i][1] + dijy*q + diky*p + dily*p)
        ij = coord[i] + edges[i][j]*p + edges[i][k]*q + (edges[i][j]+edges[i][k])*q
        ik = coord[i] + edges[i][j]*q + edges[i][k]*p + (edges[i][j]+edges[i][k])*q
        jl = coord[i] + edges[i][j]*p + edges[i][k]*q + (edges[i][j]+edges[i][k])*p
        kl = coord[i] + edges[i][j]*q + edges[i][k]*p + (edges[i][j]+edges[i][k])*p
        if edgecolors[(i,j)] == red:
            nostroke()
            fill(0)
            z = box[2] / 2
        else:            
            stroke(0)
            nofill()
            z = 0.0
        x,y = wrap(ij, box)
        oval(x-3,y-3,6,6)
        atoms.append((x,y,z))
        x,y = wrap(kl, box)
        oval(x-3,y-3,6,6)
        atoms.append((x,y,z))
        if edgecolors[(i,j)] == blue:
            nostroke()
            fill(0)
            z = box[2] / 2
        else:            
            stroke(0)
            nofill()
            z = 0.0
        x,y = wrap(ik, box)
        oval(x-3,y-3,6,6)
        atoms.append((x,y,z))
        x,y = wrap(jl, box)
        oval(x-3,y-3,6,6)
        atoms.append((x,y,z))
        #text("%s" % len(atoms), x,y)
    return np.array(atoms)


def tetrahedra(atoms,edges):
    tets = dict()
    for i in range(len(atoms)):
        for j in edges[i]:
            if i < j:
                for k in edges[i]:
                    if j < k:
                        for l in edges[i]:
                            if k < l:
                                if k in edges[j] and l in edges[j] and l in edges[k]:
                                    tets[i,j,k,l] = atoms[i] + (edges[i][j] + edges[i][k] + edges[i][l])/4.0
    return tets
                                        

def draw3d(atoms,edges):
    nofill()
    hue = 0.0
    for i in edges:
        print(i,len(edges[i]),edges[i].keys())
        for j in edges[i]:
            dx,dy,dz = edges[i][j]
            stroke(hue,1,1)
            line(atoms[i][0]+atoms[i][2]*0.03,atoms[i][1]+atoms[i][2]*0.09,
                 atoms[i][0]+dx+(atoms[i][2]+dz)*0.03,atoms[i][1]+dy+(atoms[i][2]+dz)*0.09)
            hue += (sqrt(5.0)-1.0)/2
            hue %= 1.0
    nostroke()
    i = 0
    for x,y,z in atoms:
        x += z*0.03
        y += z*0.09
        fill(0,1,1)
        oval(x-5,y-5,10,10)
        fill(0,1,0)
        text("%s" % i,x,y)
        i += 1



triangle_top = 0
triangle_left = 1
triangle_right = 2

theta = pi*15.0/180.0
ratio = 0.5 * tan(theta)
ratio2 = 0.5 / cos(theta)
matrix12 = ((cos(theta*2),sin(theta*2)),
          (-sin(theta*2),cos(theta*2)))



def argparser(arg):
    global coord, cell, waters, bondlen, density, cages
    gglen = 1.0 / (2 + 2*sqrt(3))

    p = gglen * sqrt(3)
    box2 = np.array([1.0, 1.0])
    boxz = box2[0] / sqrt(2) /23.18 * 24.29 / 2.0
    # Archimede tiling
    coord = [np.array([x,y]) for x,y in [(gglen,0.0),
                                         (gglen+2*p,0.0),
                                         (gglen+p, gglen),
                                         (0.0, p),
                                         (p,gglen+p),
                                         (2*gglen+p,gglen+p),
                                         (0.0,p+2*gglen),
                                         (gglen+p,gglen+2*p)]]
        
    edgelen = 2*gglen
    # connect nearest vertices
    edges = bond(coord,edgelen, box2)
    # connect to make shapes
    triangles,squares,edgeneighbor = findrings(edges)
    
    # yaplot size
    R = 0.02
    
    # inflation
    if arg == "1":
        coord = inflate(coord,edgelen,edges,triangles,squares,1)
        edgelen *= 2*ratio
        gglen  *= 2*ratio
        boxz *= (2*ratio)
        R *= 2*ratio
        edges = bond(coord,edgelen,box2)
        triangles,squares,edgeneighbor = findrings(edges)
    
    # print
    # print(triangles)
    #print(squares)
    #print(edgeneighbor)
    #print(coord)
    
    #2d to 3d
    box3 = np.array([1.0, 1.0, boxz])
    
    atoms = onelayer(coord,edgelen,edges,triangles,squares,edgeneighbor,box3)
    # double the layer (2 layers are minimum to define cages)
    N = atoms.shape[0]
    a = np.zeros([N*2,3])
    a[:N] = atoms
    a[N:] = atoms + np.array([0, 0, box3[2]])
    atoms = a
    box3[2] *= 2

    # use toWater
    #atoms /= box3
    #atoms -= np.floor(atoms+0.5)
    #waters = np.array([w for w in toWater(atoms, np.diag(box3), tolerance=1.7)])
    ## for GenIce
    #box3 *= 2.76 / (gglen*0.5)
    #coord="relative"
    #cell = np.diag(box3)
    #bondlen = 3
    #density=0.8
    #cages = ""
    #for atom in atoms:
    #    cages += "{0} {1} {2} {3}\n".format("cage", *atom)
    #return

    atoms -= np.floor(atoms/box3)*box3
    
    
    edges = [dict() for i in range(len(atoms))]
    for i,j,d in pl.pairs_iter(atoms/box3, gglen*2*sqrt(19.0/48.0)*1.01, np.diag(box3)):
        d = wrap(atoms[j] - atoms[i], box3)
        edges[i][j] = d
        edges[j][i] = -d
    tets = tetrahedra(atoms,edges)
        
    tetc = np.array(list(tets.values()))
    tetc -= np.floor(tetc/box3)*box3

    # reduce cell size
    # The original (level 0) structure is sqrt2 x sqrt2 of TS1 because it was generated from cage layout and it was minimal unit cell
    # Here we reduce the cell by clipping the central part of the cell
    #newbox = box3 / np.array([2**0.5, 2**0.5, 1.0])
    #N = tetc.shape[0]
    #tetc2 = np.zeros([N*2,3])
    #tetc2[:N] = tetc
    #tetc2[N:] = tetc - np.array([0.0, box3[1], 0.0])
    ## rotate 45 degree
    #tetc2 = tetc2 @ np.array([[2**0.5/2,2**0.5/2, 0.0], [-(2**0.5)/2, 2**0.5/2, 0.0], [0.0, 0.0, 1.0]])
    #tetc2 -= np.array([1e-4,1e-4,0])
    #N2 = 0
    #newtetc = []
    #for i in range(N*2):
    #    if 0 <= tetc2[i,0] < newbox[0] and 0 <= tetc2[i,1] < newbox[1]:
    #        N2 += 1
    #        newtetc.append(tetc2[i])
    #assert N2 == N/2
    #box3 = newbox
    #tetc = np.array(newtetc)

    box3 *= 2.76 / (gglen*0.5)


    # for GenIce
    waters = tetc * 2.76 / (gglen*0.5)
    coord="absolute"
    cell = np.diag(box3)
    bondlen = 3
    density=0.8

argparser("0")
