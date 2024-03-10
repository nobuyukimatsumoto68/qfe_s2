#!/usr/bin/env python3

import os
from mpmath import mp
import numpy as np
from quat import Quat
from grp_o3 import GrpElemO3

mp.dps = 30

sqrt1_2 = 1 / mp.sqrt(2)
sqrt3_4 = mp.sqrt(3) / 2
sqrt2_3 = mp.sqrt(2) / mp.sqrt(3)
sqrt1_3 = 1 / mp.sqrt(3)

C309 = mp.cos(2 * mp.pi / 5)
C425 = mp.sqrt(1 / 2 + mp.sqrt(5) / 10) / 2
C587 = mp.cos(3 * mp.pi / 10)
C688 = mp.sqrt(1 + 2 / mp.sqrt(5)) / 2
C809 = mp.cos(mp.pi / 5)


def FixElem(g):
    # apply sign convention to avoid double counting degenerate quaternions
    w = mp.chop(g.q.w)
    x = mp.chop(g.q.x)
    y = mp.chop(g.q.y)
    z = mp.chop(g.q.z)

    s = 1
    if w < 0:
        s = -1
    elif w == 0 and x < 0:
        s = -1
    elif w == 0 and x == 0 and y < 0:
        s = -1
    elif w == 0 and x == 0 and y == 0 and z < 0:
        s = -1

    q = Quat(w * s, x * s, y * s, z * s)
    return GrpElemO3(q, g.inv)


def ElemName(g):
    # make a string to identify this group element
    w = mp.nstr(mp.chop(g.q.w), 20)
    x = mp.nstr(mp.chop(g.q.x), 20)
    y = mp.nstr(mp.chop(g.q.y), 20)
    z = mp.nstr(mp.chop(g.q.z), 20)

    return ','.join([str(g.inv), w, x, y, z])


def VertexName(v):
    # make a string to identify this vertex
    x = mp.nstr(mp.chop(v[0]), 20)
    y = mp.nstr(mp.chop(v[1]), 20)
    z = mp.nstr(mp.chop(v[2]), 20)

    return ','.join([x, y, z])


def GenerateGroup(q):

    print(f'Generating group o3q{q}...')

    # list of group generators
    gen = []
    if q == 3:
        gen.append(GrpElemO3(Quat(0.5, 0, 0, sqrt3_4), 1))
        gen.append(GrpElemO3(Quat(0, sqrt2_3, 0, sqrt1_3), 1))
        gen.append(GrpElemO3(Quat(0, 0, 1.0, 0), -1))
    elif q == 4:
        gen.append(GrpElemO3(Quat(sqrt1_2, sqrt1_2, 0, 0), 1))
        gen.append(GrpElemO3(Quat(0.5, 0.5, 0.5, 0.5), -1))
    elif q == 5:
        gen.append(GrpElemO3(Quat(C809, 0, 0, C587), 1))
        gen.append(GrpElemO3(Quat(0.5, C425, C309, C688), -1))

    # start with the identity element only
    G = [GrpElemO3.Identity()]
    elem_names = {ElemName(G[0]): 1}

    while (True):
        # new elements
        new_ones = []

        # multiply all current elements times all generators
        for g in G:
            for h in gen:
                # generate a new group element
                gh = FixElem(g * h)

                # skip if it's not new
                hash_string = ElemName(gh)
                if hash_string in elem_names:
                    continue

                # add it to the list of new ones
                new_ones.append(gh)
                elem_names[hash_string] = 1

        # exit when no new elements were generated
        if (len(new_ones) == 0):
            break

        # add new elements to the list
        G += new_ones

    print(f'Group has {len(G)} elements')

    # generate a set of vertices
    vertices = []
    vertex_names = {}

    z_hat = np.array([0.0, 0.0, 1.0])

    for g in G:
        v = g * z_hat
        hash_string = VertexName(v)
        if hash_string in vertex_names:
            continue

        vertices.append(v)
        vertex_names[hash_string] = 1

    # determine the edge length
    edge_length = 2.0
    for s in range(1, len(vertices)):
        length = np.linalg.norm(vertices[s] - vertices[0])
        if length < edge_length:
            edge_length = length

    # generate a set of faces
    faces = []
    for s1 in range(len(vertices)):
        for s2 in range(s1 + 1, len(vertices)):
            length12 = np.linalg.norm(vertices[s1] - vertices[s2])
            if not mp.almosteq(length12, edge_length):
                continue
            for s3 in range(s2 + 1, len(vertices)):
                length13 = np.linalg.norm(vertices[s1] - vertices[s3])
                if not mp.almosteq(length13, edge_length):
                    continue
                length23 = np.linalg.norm(vertices[s2] - vertices[s3])
                if not mp.almosteq(length23, edge_length):
                    continue

                faces.append([s1, s2, s3])

    print(f'Polyhedron has {len(vertices)} vertices and {len(faces)} faces')
    print('Writing group data to file...\n')

    # write group elements to a data file
    file = open(f'elem/o3q{q}.dat', 'w')
    for i, g in enumerate(G):
        w = mp.chop(g.q.w)
        x = mp.chop(g.q.x)
        y = mp.chop(g.q.y)
        z = mp.chop(g.q.z)
        file.write('%d %+d ' % (i, g.inv))
        file.write('%+.18f %+.18f %+.18f %+.18f\n' % (w, x, y, z))
    file.close()

    # write polyhedron lattice to a data file
    file = open(f'lattice/o3q{q}.dat', 'w')
    file.write('begin_sites\n')
    file.write(f'n_sites {len(vertices)}\n')
    for i, v in enumerate(vertices):
        phi = 0.0
        cos_theta = mp.chop(v[2])
        if mp.almosteq(cos_theta, 1.0):
            theta = 0.0
        elif mp.almosteq(cos_theta, 1.0):
            theta = mp.pi
        else:
            theta = mp.acos(cos_theta)
            phi = mp.atan2(mp.chop(v[1]), mp.chop(v[0]))

        theta = mp.chop(theta)
        phi = mp.chop(phi)
        file.write('%d 1.0 0 %+.18f %+.18f\n' % (i, theta, phi))
    file.write('end_sites\n')

    file.write('begin_faces\n')
    file.write(f'n_faces {len(faces)}\n')
    for i, f in enumerate(faces):
        file.write('%d 1.0 %d %d %d\n' % (i, f[0], f[1], f[2]))
    file.write('end_faces\n')
    file.close()


if __name__ == '__main__':
    os.makedirs('elem', exist_ok=True)
    os.makedirs('lattice', exist_ok=True)
    os.makedirs('site_g', exist_ok=True)
    GenerateGroup(3)
    GenerateGroup(4)
    GenerateGroup(5)
