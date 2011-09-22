#Put the initial vertex
#subdivide with dodecagon and register new vertices
from math import *

blue = 1
red  = -1

triangle_top = 0
triangle_left = 1
triangle_right = 2

theta = pi*15.0/180.0
ratio = 0.5 * tan(theta)
ratio2 = 0.5 / cos(theta)
matrix12 = ((cos(theta*2),sin(theta*2)),
          (-sin(theta*2),cos(theta*2)))

def mul(mat,vec):
    return ((mat[0][0]*vec[0]+mat[0][1]*vec[1]),
            (mat[1][0]*vec[0]+mat[1][1]*vec[1]))

def add(vec1,vec2):
    return (vec1[0]+vec2[0],vec1[1]+vec2[1])

def sub(vec1,vec2):
    v = [0] * len(vec1)
    for i in range(len(vec1)):
        v[i] = vec[1] - vec[2]
    return v

def wrap(c):
    v = [0] * len(c)
    for i in range(len(c)):
        v[i] =  c[i]-floor(c[i]/box[i]+0.5)*box[i]
    return v

def sub_pbc(vec1,vec2):
    v = [0] * len(vec1)
    for i in range(len(vec1)):
        v[i] = vec1[i]-vec2[i]
    return wrap(v)


def drawpoly( coord ):
    for i in range(len(coord)):
        line(coord[i-1][0],coord[i-1][1],
             coord[i  ][0],coord[i  ][1])

def subdiv_triangle(v):
    for i in range(3):
        dijx,dijy = edges[v[i]][v[i-1]]
        dikx,diky = edges[v[i]][v[i-2]]
        p = (coord[v[i]][0] + (dijx+dikx)*(0.5 - ratio / sqrt(3.0)),
             coord[v[i]][1] + (dijy+diky)*(0.5 - ratio / sqrt(3.0)))
        coord.append(p)

def subdiv_square(v):
    for i in range(4):
        #relative diagonal vector of the square
        d = (edges[v[i]][v[i-1]][0]+edges[v[i-1]][v[i-2]][0],
                     edges[v[i]][v[i-1]][1]+edges[v[i-1]][v[i-2]][1])
        p0 = (d[0]/sqrt(2)*ratio2,
              d[1]/sqrt(2)*ratio2)
        p = add(coord[v[i]],p0)
        coord.append(p)
        p1 = mul(matrix12,p0)
        p = add(coord[v[i]],p1)
        coord.append(p)
        

def hexagon(center,edge,parity):
    r = edge * ratio * 2.0
    for i in range(6):
        coord.append(add(center,(r*cos((i*2+parity)*pi/6),r*sin((i*2+parity)*pi/6))))


def bond(coord,thres):
    edges = dict()
    for i in range(len(coord)):
        edges[i] = dict()
    for i in range(len(coord)):
        for j in range(i+1,len(coord)):
            dx,dy = sub_pbc(coord[j],coord[i])
            if dx**2+dy**2 < thres**2 * 1.01:
                edges[i][j] = (dx,dy)
                edges[j][i] = (-dx,-dy)
    return edges


def dictoflist_add(dol,key,value):
    if not dol.has_key(key):
        dol[key] = []
    dol[key].append(value)

def bond3d_fast(coord,thres):
    edges = dict()
    for i in range(len(coord)):
        edges[i] = dict()
    ndiv = floor(box[0]/thres),floor(box[1]/thres),floor(box[2]/thres)
    ndiv = map(int,ndiv)
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
                                        x,y,z =sub_pbc(coord[i],coord[j])
                                        if x**2 + y**2 + z**2 < thres**2:
                                            edges[j][i] = x,y,z
    return edges


def bond3d(coord,thres):
    edges = dict()
    for i in range(len(coord)):
        edges[i] = dict()
    for i in range(len(coord)):
        for j in range(len(coord)):
            if i != j:
                x,y,z = sub_pbc(coord[i],coord[j])
                if x**2 + y**2 + z**2 < thres**2:
                    edges[j][i] = x,y,z
    return edges
            

SAMECOLOR = 1
DIFFCOLOR = -1


def setcolors(en,ec,e,col):
    if ec.has_key(e):
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
                    if edges[j].has_key(k):
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
                            if edges[k].has_key(l):
                                if l > i:
                                    if not edges[i].has_key(l):
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
 

def inflate(coord,edge,triangles,squares,depth):
    clen = len(coord)
    if depth == 0:
        for i in range(clen):
            hexagon(coord[i],edge,0) #ast number is the direction of hexagon
    else:
        for i in range(clen):
            hexagon(coord[i],edge,0) #i%2#ast number is the direction of hexagon
    for triangle in triangles:
        subdiv_triangle(triangle)
    for square in squares:
        subdiv_square(square)
    return coord

def onelayer(coord,edge,tr,sq,edgeneighbor):
    atoms = []
    i = 0
    for c in coord:
        x,y = wrap(c)        
        atoms.append((x,y,box[2]/4))
        atoms.append((x,y,box[2]*3/4))
        i += 1
    #return
    edgecolors = dict()
    i = edges.keys()[0]
    j = edges[i].keys()[0]
    if i > j:
        i,j = j,i
    setcolors(edgeneighbor,edgecolors,(i,j),red)
    for i in edges:
        for j in edges[i]:
            if i < j:
                dx,dy = edges[i][j]
                ij = (coord[i][0]+dx*0.5,
                      coord[i][1]+dy*0.5)
                x,y = wrap(ij)
                if edgecolors[(i,j)] == red:
                    atoms.append((x,y,0.0))
                else:
                    atoms.append((x,y,box[2]/2))
                x,y = wrap(coord[i])
    for i,j,k in tr:
        dijx,dijy = edges[i][j]
        dikx,diky = edges[i][k]
        center = (coord[i][0]+(dijx+dikx)/3,
                  coord[i][1]+(dijy+diky)/3)
        x,y = wrap(center)
        if edgecolors[(i,j)] == red:
            z = box[2]/2
        else:            
            z = 0.0
        atoms.append((x,y,z))
    for i,j,l,k in sq:
        p = 3.0/8.0
        q = 1.0/8.0
        dijx,dijy = edges[i][j]
        dikx,diky = edges[i][k]
        dilx,dily = dijx+dikx, dijy+diky
        ij = (coord[i][0] + dijx*p + dikx*q + dilx*q,
              coord[i][1] + dijy*p + diky*q + dily*q)
        ik = (coord[i][0] + dijx*q + dikx*p + dilx*q,
              coord[i][1] + dijy*q + diky*p + dily*q)
        jl = (coord[i][0] + dijx*p + dikx*q + dilx*p,
              coord[i][1] + dijy*p + diky*q + dily*p)
        kl = (coord[i][0] + dijx*q + dikx*p + dilx*p,
              coord[i][1] + dijy*q + diky*p + dily*p)
        if edgecolors[(i,j)] == red:
            z = box[2] / 2
        else:            
            z = 0.0
        x,y = wrap(ij)
        atoms.append((x,y,z))
        x,y = wrap(kl)
        atoms.append((x,y,z))
        if edgecolors[(i,j)] == blue:
            z = box[2] / 2
        else:            
            z = 0.0
        x,y = wrap(ik)
        atoms.append((x,y,z))
        x,y = wrap(jl)
        atoms.append((x,y,z))
    return atoms


def tetrahedra(atoms,edges):
    tets = dict()
    for i in range(len(atoms)):
        for j in edges[i]:
            if i < j:
                for k in edges[i]:
                    if j < k:
                        for l in edges[i]:
                            if k < l:
                                if (edges[j].has_key(k) and
                                    edges[j].has_key(l) and
                                    edges[k].has_key(l) ):
                                    dx = (edges[i][j][0]+edges[i][k][0]+edges[i][l][0])/4.0
                                    dy = (edges[i][j][1]+edges[i][k][1]+edges[i][l][1])/4.0
                                    dz = (edges[i][j][2]+edges[i][k][2]+edges[i][l][2])/4.0
                                    tets[(i,j,k,l)] = (atoms[i][0]+dx,
                                                       atoms[i][1]+dy,
                                                       atoms[i][2]+dz)
    return tets
                                        


zoom = 23.18*sqrt(2)
p = zoom * sqrt(3)
box = [2*zoom + 2*p, 2*zoom + 2*p, 0.0]
coord = [(zoom,0.0), (zoom+2*p,0.0), (zoom+p, zoom),
         (0.0, p), (p,zoom+p), (2*zoom+p,zoom+p),
         (0.0,p+2*zoom), (zoom+p,zoom+2*p)]
        
edge = 2*zoom
box[2] = box[0] / sqrt(2) /23.18 * 24.29 / 2.0
print box[2],edge

edges = bond(coord,edge)
triangles,squares,edgeneighbor = findrings(edges)

coord = inflate(coord,edge,triangles,squares,1)
edge *= 2*ratio
box[2] *= 2*ratio
edges = bond(coord,edge)
triangles,squares,edgeneighbor = findrings(edges)

#coord = inflate(coord,edge,triangles,squares,1)
#edge *= 2*ratio
#edges = bond(coord,edge)
#triangles,squares,edgeneighbor = findrings(edges)

thres = edge*sqrt(19.0/48.0)*1.01
atoms = onelayer(coord,edge,triangles,squares,edgeneighbor)
natoms = len(atoms)
#double the layers
for i in range(natoms):
    atoms.append((atoms[i][0],atoms[i][1],atoms[i][2]+box[2]))
#    atoms.append((atoms[i][0],atoms[i][1],atoms[i][2]+box[2]*2))
box[2] *= 2 #3
edges = bond3d_fast(atoms,thres)
tets = tetrahedra(atoms,edges)
#print "#box"
#print box[0],box[1],box[2]
#print "#atoms"
#print len(atoms)
#for i in range(len(atoms)):
#    print "G",atoms[i][0],atoms[i][1],atoms[i][2]
#print len(tets)
#for tet in tets:
#    c = tets[tet]
#    print "W",c[0],c[1],c[2]
l=2*zoom*(1+sqrt(3))/(2*(2+sqrt(3)))
translate(zoom*5,zoom*5)
nofill()
stroke(0)
rect(-zoom,-l,zoom*2,l*2)
print "@BOX3"
print zoom*2,l*2,box[2]
print "@AR3A"
s = ""
n = 0
for i in range(len(atoms)):
    oval(atoms[i][0]-2,atoms[i][1]-2,4,4)
    if -zoom-0.8 < atoms[i][0] < zoom-0.8 and -l-0.8 < atoms[i][1] < l-0.8:
        s += "%s %s %s\n" % (atoms[i][0],atoms[i][1],atoms[i][2])
        n += 1
print n
print s,


    