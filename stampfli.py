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
    return (vec1[0]-vec2[0],vec1[1]-vec2[1])

def wrap(c):
    return (c[0]-floor(c[0]/box[0]+0.5)*box[0],
            c[1]-floor(c[1]/box[1]+0.5)*box[1])

def sub_pbc(vec1,vec2):
    dx = vec1[0]-vec2[0]
    dy = vec1[1]-vec2[1]
    return wrap((dx,dy))


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
            #hexagon(coord[i],edge,i%2) #ast number is the direction of hexagon
            hexagon(coord[i],edge,0) #ast number is the direction of hexagon
    for triangle in triangles:
        subdiv_triangle(triangle)
    for square in squares:
        subdiv_square(square)
    return coord

def draw(coord,edge,tr,sq,edgeneighbor):
    i = 0
    for c in coord:
        stroke(0)
        nofill()
        x,y = wrap(c)        
        oval( x-4,y-4,8,8 )
        nostroke()
        fill(0)
        oval( x-2,y-2,4,4 )
        #text("%s" % i,c[0],c[1] )
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
                    stroke(0)
                    nofill()
                    oval( x-3,y-3,6,6 )
                    stroke(0,1,1)
                else:
                    nostroke()
                    fill(0)
                    oval( x-3,y-3,6,6 )
                    stroke(0.666,1,1)
                x,y = wrap(coord[i])
                line(x,y,x+dx,y+dy)
    for i,j,k in tr:
        dijx,dijy = edges[i][j]
        dikx,diky = edges[i][k]
        center = (coord[i][0]+(dijx+dikx)/3,
                  coord[i][1]+(dijy+diky)/3)
        x,y = wrap(center)
        if edgecolors[(i,j)] == red:
            nostroke()
            fill(0)
        else:            
            stroke(0)
            nofill()

        oval(x-3,y-3,6,6)
        #text("(%s,%s,%s)" % (i,j,k), center[0]+5,center[1]+5)
        #print (i,j,k),center
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
            nostroke()
            fill(0)
        else:            
            stroke(0)
            nofill()
        x,y = wrap(ij)
        oval(x-3,y-3,6,6)
        x,y = wrap(kl)
        oval(x-3,y-3,6,6)
        if edgecolors[(i,j)] == blue:
            nostroke()
            fill(0)
        else:            
            stroke(0)
            nofill()
        x,y = wrap(ik)
        oval(x-3,y-3,6,6)
        x,y = wrap(jl)
        oval(x-3,y-3,6,6)
        #text("(%s,%s,%s,%s)" % (i,j,l,k), ij[0],ij[1])




nofill()
stroke(0)
colormode(HSB)
#coord = [(100.0,100.0),(500.0,100.0),(300.0,450.0)]
#subdiv_triangle((0,1,2))
zoom = 100.0
translate(5*zoom,5*zoom)
p = zoom * sqrt(3)
box = (2*zoom + 2*p, 2*zoom + 2*p)
coord = [(zoom,0.0), (zoom+2*p,0.0), (zoom+p, zoom),
         (0.0, p), (p,zoom+p), (2*zoom+p,zoom+p),
         (0.0,p+2*zoom), (zoom+p,zoom+2*p)]
        
edge = 2*zoom

edges = bond(coord,edge)
triangles,squares,edgeneighbor = findrings(edges)

coord = inflate(coord,edge,triangles,squares,1)
edge *= 2*ratio
edges = bond(coord,edge)
triangles,squares,edgeneighbor = findrings(edges)

#coord = inflate(coord,edge,triangles,squares,1)
#edge *= 2*ratio
#edges = bond(coord,edge)
#triangles,squares,edgeneighbor = findrings(edges)

l=2*zoom*(1+sqrt(3))/(2*(2+sqrt(3)))
rect(-zoom,-l,2*zoom,l*2)
draw(coord,edge,triangles,squares,edgeneighbor)

print len(coord)
