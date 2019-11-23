#!/usr/bin/env python

def dictoflist_add(dol,key,value):
    if not dol.has_key(key):
        dol[key] = []
    dol[key].append(value)

def sub(vec1,vec2):
    v = [0] * len(vec1)
    for i in range(len(vec1)):
        v[i] = vec[1] - vec[2]
    return v

def wrap(c,box):
    v = [0] * len(c)
    for i in range(len(c)):
        v[i] =  c[i]-floor(c[i]/box[i]+0.5)*box[i]
    return v

def sub_pbc(vec1,vec2,box):
    v = [0] * len(vec1)
    for i in range(len(vec1)):
        v[i] = vec1[i]-vec2[i]
    return wrap(v,box)


from math import *

def bond3d_fast(coord,thres,box):
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
                            if resident.has_key((ix,iy,iz)) and resident.has_key((jx,jy,jz)):
                                for i in resident[(ix,iy,iz)]:
                                    for j in resident[(jx,jy,jz)]:
                                        if i != j:
                                            x,y,z =sub_pbc(coord[i],coord[j],box)
                                        if x**2 + y**2 + z**2 < thres**2:
                                            edges[j][i] = x,y,z
    return edges

import sys

#read water (type==1) only
natom = 0
box = [0.0] * 3
coord = []
while True:
    line = sys.stdin.readline()
    if line == "":
        break
#    columns = line.split()
    if line == "ITEM: BOX BOUNDS\n":
        for i in range(3):
            line = sys.stdin.readline()
            columns = line.split()
            columns = map(float,columns)
            box[i] = columns[1] - columns[0]
    if line == "ITEM: NUMBER OF ATOMS\n":
        line = sys.stdin.readline()
        natom = int(line)
    if line == "ITEM: ATOMS\n":
        for i in range(natom):
            line = sys.stdin.readline()
            columns = line.split()
            columns[0:2] = map(int,columns[0:2])
            columns[2:5] = map(float,columns[2:5])
            if columns[1] == 1:
                coord.append(columns[2:5])
                for i in range(3):
                    coord[-1][i] *= box[i]

edges = bond3d_fast(coord,3.0,box)
for i in edges:
    print len(edges[i])
