"""
Approximant of a dodecagonal quasicrystal designed by Stampfli.
"""

# Put the initial vertex
# subdivide with dodecagon and register new vertices
import numpy as np
import pairlist as pl
from genice2.FrankKasper import toWater
import genice2.lattices
from logging import getLogger
from math import pi, sin, cos, tan, sqrt

desc = {
    "ref": {"stampfli": "Stampfli 1986"},
    "usage": __doc__,
    "brief": "Approximant of Stampfli's dodecagonal quasicrystal",
    "test": ("[1]", "[2]",)
}


def wrap(c, box):
    s = c.shape[0]
    return c - np.floor(c / box[:s] + 0.5) * box[:s]


def subdiv_triangle(v, edges, coord):
    for i in range(3):
        dijx, dijy = edges[v[i]][v[i - 1]]
        dikx, diky = edges[v[i]][v[i - 2]]
        p = (coord[v[i]][0] + (dijx + dikx) * (0.5 - ratio / sqrt(3.0)),
             coord[v[i]][1] + (dijy + diky) * (0.5 - ratio / sqrt(3.0)))
        coord.append(np.array(p))


def subdiv_square(v, edges, coord):
    for i in range(4):
        # relative diagonal vector of the square
        d = edges[v[i]][v[i - 1]] + edges[v[i - 1]][v[i - 2]]
        p0 = d / sqrt(2) * ratio2
        p = coord[v[i]] + p0
        coord.append(np.array(p))
        p1 = matrix12 @ p0
        p = coord[v[i]] + p1
        coord.append(np.array(p))


def hexagon(center, edge, parity, coord):
    r = edge * ratio * 2.0
    for i in range(6):
        coord.append(
            center + np.array([r * cos((i * 2 + parity) * pi / 6), r * sin((i * 2 + parity) * pi / 6)]))


def bond(coord, thres, box):
    edges = dict()
    for i in range(len(coord)):
        edges[i] = dict()
    for i in range(len(coord)):
        for j in range(i + 1, len(coord)):
            # dx,dy = sub_pbc(coord[j],coord[i])
            d = wrap(coord[j] - coord[i], box)
            # if dx**2+dy**2 < thres**2 * 1.01:
            if d @ d < thres**2 * 1.01:
                edges[i][j] = d  # (dx,dy)
                edges[j][i] = -d  # (-dx,-dy)
    return edges


def dictoflist_add(dol, key, value):
    if key not in dol:
        dol[key] = []
    dol[key].append(value)


def bond3d_fast(coord, thres, box):
    edges = dict()
    for i in range(len(coord)):
        edges[i] = dict()
    ndiv = [
        int(x) for x in (
            floor(
                box[0] /
                thres),
            floor(
                box[1] /
                thres),
            floor(
                box[2] /
                thres))]
    binw = box[0] / ndiv[0], box[1] / ndiv[1], box[2] / ndiv[2]
    resident = dict()
    for i in range(len(coord)):
        c = coord[i]
        ix, iy, iz = int(floor(c[0] /
                               binw[0])), int(floor(c[1] /
                                                    binw[1])), int(floor(c[2] /
                                                                         binw[2]))
        if ix < 0:
            ix += ndiv[0]
        if iy < 0:
            iy += ndiv[1]
        if iz < 0:
            iz += ndiv[2]
        dictoflist_add(resident, (ix, iy, iz), i)
    for ix in range(ndiv[0]):
        for iy in range(ndiv[1]):
            for iz in range(ndiv[2]):
                for dx in range(-1, 2):
                    for dy in range(-1, 2):
                        for dz in range(-1, 2):
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
                            # print jx,jy,jz,resident[(jx,jy,jz)]
                            for i in resident[(ix, iy, iz)]:
                                for j in resident[(jx, jy, jz)]:
                                    if i != j:
                                        xyz = wrap(coord[i] - coord[j], box)
                                        if xyz @ xyz < thres**2:
                                            edges[j][i] = xyz
    return edges


def bond3d(coord, thres, box):
    edges = dict()
    for i in range(len(coord)):
        edges[i] = dict()
    for i in range(len(coord)):
        for j in range(len(coord)):
            if i != j:
                xyz = wrap(coord[i] - coord[j], box)
                if xyz @ xyz < thres**2:
                    edges[j][i] = xyz
    return edges


SAMECOLOR = 1
DIFFCOLOR = -1



def setcolors(en, ec, e, col):
    queue = []
    queue.append([en, e, col])

    while len(queue) > 0:
        en, e, col = queue.pop(0)
        if e in ec:
            continue
        ec[e] = col
        # print e,col
        for i, j, parity in en[e]:
            queue.append([en, (i, j), col * parity])



def _setcolors(en, ec, e, col):
    if e in ec:
        return
    ec[e] = col
    # print e,col
    for i, j, parity in en[e]:
        _setcolors(en, ec, (i, j), col * parity)


def findrings(edges):
    # look for rings
    sq = []
    tr = []
    edgeneighbor = dict()
    for i in edges:
        for j in edges[i]:
            for k in edges[i]:
                if i < j and i < k and j < k:
                    if k in edges[j]:
                        # must be compact
                        dx = edges[i][j][0] + edges[j][k][0] + edges[k][i][0]
                        dy = edges[i][j][1] + edges[j][k][1] + edges[k][i][1]
                        if abs(dx) < 0.001 and abs(dy) < 0.001:
                            tr.append((i, j, k))
                            dictoflist_add(
                                edgeneighbor, (i, j), (i, k, SAMECOLOR))
                            dictoflist_add(
                                edgeneighbor, (i, j), (j, k, SAMECOLOR))
                            dictoflist_add(
                                edgeneighbor, (i, k), (i, j, SAMECOLOR))
                            dictoflist_add(
                                edgeneighbor, (i, k), (j, k, SAMECOLOR))
                            dictoflist_add(
                                edgeneighbor, (j, k), (i, j, SAMECOLOR))
                            dictoflist_add(
                                edgeneighbor, (j, k), (i, k, SAMECOLOR))
                    else:
                        for l in edges[j]:
                            if l in edges[k]:
                                if l > i:
                                    if l not in edges[i]:
                                        # must be compact
                                        dx = edges[i][j][0] + edges[j][l][0] + \
                                            edges[l][k][0] + edges[k][i][0]
                                        dy = edges[i][j][1] + edges[j][l][1] + \
                                            edges[l][k][1] + edges[k][i][1]
                                        if abs(dx) < 0.001 and abs(dy) < 0.001:
                                            sq.append((i, j, l, k))
                                            kk, kl = k, l
                                            if k > l:
                                                kk, kl = l, k
                                            jj, jl = j, l
                                            if j > l:
                                                jj, jl = l, j
                                            dictoflist_add(
                                                edgeneighbor, (i, j), (i, k, DIFFCOLOR))
                                            dictoflist_add(
                                                edgeneighbor, (i, k), (i, j, DIFFCOLOR))
                                            dictoflist_add(
                                                edgeneighbor, (i, k), (kk, kl, DIFFCOLOR))
                                            dictoflist_add(
                                                edgeneighbor, (kk, kl), (i, k, DIFFCOLOR))
                                            dictoflist_add(
                                                edgeneighbor, (kk, kl), (jj, jl, DIFFCOLOR))
                                            dictoflist_add(
                                                edgeneighbor, (jj, jl), (kk, kl, DIFFCOLOR))
                                            dictoflist_add(
                                                edgeneighbor, (i, j), (jj, jl, DIFFCOLOR))
                                            dictoflist_add(
                                                edgeneighbor, (jj, jl), (i, j, DIFFCOLOR))
                                            dictoflist_add(
                                                edgeneighbor, (i, j), (kk, kl, SAMECOLOR))
                                            dictoflist_add(
                                                edgeneighbor, (i, k), (jj, jl, SAMECOLOR))
                                            dictoflist_add(
                                                edgeneighbor, (kk, kl), (i, j, SAMECOLOR))
                                            dictoflist_add(
                                                edgeneighbor, (jj, jl), (i, k, SAMECOLOR))
    return tr, sq, edgeneighbor


def inflate(coord, edgelen, edges, triangles, squares, depth):
    clen = len(coord)
    if depth == 0:
        for i in range(clen):
            # ast number is the direction of hexagon
            hexagon(coord[i], edgelen, 0, coord)
    else:
        for i in range(clen):
            hexagon(
                coord[i], edgelen, i %
                2, coord)  # ast number is the direction of hexagon
    for triangle in triangles:
        subdiv_triangle(triangle, edges, coord)
    for square in squares:
        subdiv_square(square, edges, coord)
    return coord


blue = 1
red = -1


def onelayer(coord, edgelen, edges, tr, sq, edgeneighbor, box):
    """
    2D layer to 3D layer
    """
    ratoms = []
    for c in coord:
        xy = c / box[:2]
        xy -= np.floor(xy) 
        x, y = xy
        ratoms.append((x, y, 1 / 4))
        ratoms.append((x, y, 3 / 4))

    edgecolors = dict()
    i = list(edges.keys())[0]
    j = list(edges[i].keys())[0]
    if i > j:
        i, j = j, i
    setcolors(edgeneighbor, edgecolors, (i, j), red)
    for i in edges:
        for j in edges[i]:
            if i < j:
                dx, dy = edges[i][j]
                # ij = (coord[i][0]+dx*0.5,
                #      coord[i][1]+dy*0.5)
                ij = coord[i] + edges[i][j] * 0.5
                xy = ij / box[:2]
                xy -= np.floor(xy)
                x, y = xy
                if edgecolors[i, j] == red:
                    ratoms.append((x, y, 0))
                else:
                    ratoms.append((x, y, 1 / 2))
    for i, j, k in tr:
        center = coord[i] + (edges[i][j] + edges[i][k]) / 3
        xy = center / box[:2]
        xy -= np.floor(xy)
        x, y = xy
        if edgecolors[i, j] == red:
            z = 1 / 2
        else:
            z = 0
        ratoms.append((x, y, z))
    for i, j, l, k in sq:
        p = 3.0 / 8.0
        q = 1.0 / 8.0
        ij = coord[i] + edges[i][j] * p + edges[i][k] * \
            q + (edges[i][j] + edges[i][k]) * q
        ik = coord[i] + edges[i][j] * q + edges[i][k] * \
            p + (edges[i][j] + edges[i][k]) * q
        jl = coord[i] + edges[i][j] * p + edges[i][k] * \
            q + (edges[i][j] + edges[i][k]) * p
        kl = coord[i] + edges[i][j] * q + edges[i][k] * \
            p + (edges[i][j] + edges[i][k]) * p
        if edgecolors[(i, j)] == red:
            z = 1 / 2
        else:
            z = 0.0
        xy = ij / box[:2]
        xy -= np.floor(xy)
        x, y = xy
        ratoms.append((x, y, z))
        xy = kl / box[:2]
        xy -= np.floor(xy)
        x, y = xy
        ratoms.append((x, y, z))
        if edgecolors[i, j] == blue:
            z = 1 / 2
        else:
            z = 0.0
        xy = ik / box[:2]
        xy -= np.floor(xy)
        x, y = xy
        ratoms.append((x, y, z))
        xy = jl / box[:2]
        xy -= np.floor(xy)
        x, y = xy
        ratoms.append((x, y, z))
    return np.array(ratoms) * box


def tetrahedra(atoms, edges):
    tets = dict()
    for i in range(len(atoms)):
        for j in edges[i]:
            if i < j:
                for k in edges[i]:
                    if j < k:
                        for l in edges[i]:
                            if k < l:
                                if k in edges[j] and l in edges[j] and l in edges[k]:
                                    tets[i, j, k, l] = atoms[i] + \
                                        (edges[i][j] + edges[i][k] + edges[i][l]) / 4.0
    return tets




triangle_top = 0
triangle_left = 1
triangle_right = 2

theta = pi * 15.0 / 180.0
ratio = 0.5 * tan(theta)
ratio2 = 0.5 / cos(theta)
matrix12 = ((cos(theta * 2), sin(theta * 2)),
            (-sin(theta * 2), cos(theta * 2)))


class Lattice(genice2.lattices.Lattice):
    def __init__(self, **kwargs):
        logger = getLogger()
        assert len(kwargs) > 0, desc["usage"]

        self.gen = 1
        for k, v in kwargs.items():
            if k == 'generation':
                self.gen = int(v)
            elif v:  # in case only the char string is given
                self.gen = int(k)
            else:
                logger.error(f"Unknown option for stampfli plugin: {k}={v}")

        gglen = 1.0 / (2 + 2 * sqrt(3))

        p = gglen * sqrt(3)
        box2 = np.array([1.0, 1.0])
        boxz = box2[0] / sqrt(2) / 23.18 * 24.29 / 2.0
        # Archimede tiling
        coord = [np.array([x, y]) for x, y in [(gglen, 0.0),
                                            (gglen + 2 * p, 0.0),
                                            (gglen + p, gglen),
                                            (0.0, p),
                                            (p, gglen + p),
                                            (2 * gglen + p, gglen + p),
                                            (0.0, p + 2 * gglen),
                                            (gglen + p, gglen + 2 * p)]]

        edgelen = 2 * gglen
        # connect nearest vertices
        edges = bond(coord, edgelen, box2)
        # connect to make shapes
        triangles, squares, edgeneighbor = findrings(edges)

        # inflation
        for gen in range(self.gen):
            coord = inflate(coord, edgelen, edges, triangles, squares, 1)
            edgelen *= 2 * ratio
            gglen *= 2 * ratio
            boxz *= (2 * ratio)
            # R *= 2 * ratio
            edges = bond(coord, edgelen, box2)
            triangles, squares, edgeneighbor = findrings(edges)

        # 2d to 3d
        box3 = np.array([1.0, 1.0, boxz])

        atoms = onelayer(
            coord,
            edgelen,
            edges,
            triangles,
            squares,
            edgeneighbor,
            box3)
        # double the layer (2 layers are minimum to define cages)
        atoms = np.vstack([atoms, atoms + np.array([0, 0, box3[2]])])
        box3[2] *= 2

        atoms -= np.floor(atoms / box3) * box3

        edges = [dict() for i in range(len(atoms))]
        for i, j, d in pl.pairs_iter(
                atoms / box3, gglen * 2 * sqrt(19.0 / 48.0) * 1.01, np.diag(box3)):
            d = wrap(atoms[j] - atoms[i], box3)
            edges[i][j] = d
            edges[j][i] = -d
        tets = tetrahedra(atoms, edges)

        tetc = np.array(list(tets.values()))
        tetc -= np.floor(tetc / box3) * box3

        box3 *= 2.76 / (gglen * 0.5)

        # for GenIce
        self.waters = tetc * 2.76 / (gglen * 0.5)
        self.coord = "absolute"
        self.cell = np.diag(box3)
        self.bondlen = 3
        self.density = 0.8


