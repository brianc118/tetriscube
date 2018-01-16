"""
Solver for the tetris cube problem

"""

from datetime import datetime
from pprint import pprint
from functools import reduce
from cubeviewer import *
import json
import pulp
import math

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

probfile = 'standardcube.json'

def mmin(m, i, j):
    return [r[:j] + r[j + 1:] for r in (m[:i] + m[i + 1:])]

def det(m):
    """ Returns determinant of matrix """
    if len(m) == 2:
        return m[0][0] * m[1][1] - m[0][1] * m[1][0]
    d = 0
    for c in range(len(m)):
        d += ((-1) ** c) * m[0][c] * det(mmin(m, 0, c))
    return d

def permutations(l):
    """ Returns permutations of elements of list l """
    if len(l) <= 0:
        raise ValueError("List to permute is empty")
    elif len(l) == 1:
        return [[l[0]]]
    else:
        p = []
        for i in range(len(l)):
            for sp in permutations(l[:i] + l[i + 1:]):
                p.append([l[i]] + sp)
        return p

def getmat(p):
    """ Returns transformation corresponding to permutation list p """
    n = len(p)
    if min(p) < 0 or max(p) >= n:
        raise ValueError("Invalid permutation list")
    m = []
    for i in range(n):
        m.append([1 if j == p[i] else 0 for j in range(n)])
    return m


def point2int(point, dimensions):
    x = 0
    n = math.ceil(math.log2(max(dimensions)))
    for i in point:
        x <<= n
        x |= i
    return x

def int2point(x, dimensions):
    n = math.ceil(math.log2(max(dimensions)))
    ndim = len(dimensions)
    point = []
    for i in range(ndim):
        point.append(x & (2 ** n - 1))
        x >>= n
    return list(reversed(point))

def orthrotations(points, remove_identical=False):
    """
    Returns a list of all orthogonal simple rotations of points.

    """
    def points2mutable(points):
        com = [sum([p[i] for p in points]) / len(points) for i in range(len(points[0]))]
        # make com = origin
        return tuple(sorted([tuple([p[i] - com[i] for i in range(len(points[0]))]) for p in points]))

    if len(points) <= 0:
        raise ValueError("No points given")
    else:
        ans = []
        n = len(points[0])
        perms = permutations(list(range(n)))
        ids = set()

        for i in range(2**n):
            bmap = [1 if i >> k & 1 else -1 for k in range(0, n)]
            for o in perms:
                if det(getmat(o)) * reduce(lambda x, y: x * y, bmap) != 1:
                    continue
                base = []
                for p in points:
                    base.append([0] * n)
                    for j in range(n):
                        base[-1][j] = p[o[j]] * bmap[j]
                id = points2mutable(base)
                if id in ids:
                    continue
                ids.add(id)
                ans.append(base)
        return ans


# enumerate piece where piece is an n-orthotope
def enumpiece(piece, dimensions, rotate=True, fcoord=lambda x: x, fpiece=lambda x: x):
    if len(piece) <= 0:
        raise ValueError("Piece has no volume")

    n = len(piece[0])  # no. of dimensions
    if len(dimensions) != n:
        raise ValueError("Bad dimensions")

    enum = []
    if rotate:
        # n-D simple rotation about plane
        for p in orthrotations(piece, remove_identical=True):
            enum.extend(enumpiece(p, dimensions, rotate=False, fcoord=fcoord, fpiece=fpiece))
        return enum

    pmin, pmax = [], []

    for i in range(n):
        pmin.append(min(piece, key=lambda x: x[i])[i])
        pmax.append(max(piece, key=lambda x: x[i])[i])

    mincoord, maxcoord = [0] * n, [dimensions[i] for i in range(n)]

    def nestedshift(d, clones, shift=n * [0]):
        if d < 0:
            raise ValueError("Invalid recursion depth d")
        elif d == 0:
            newpiece = []
            for p in piece:
                tempcoord = [p[i] + shift[i] for i in range(n)]
                newpiece.append(fcoord(tempcoord))
            clones.append(fpiece(newpiece))
        else:
            d -= 1
            for i in range(mincoord[d] - pmin[d], maxcoord[d] - pmax[d]):
                newshift = shift[:]
                newshift[d] = i
                nestedshift(d, clones=clones, shift=newshift)

    nestedshift(d=n, clones=enum)
    return enum


with open(probfile) as pjson:
    pd = json.load(pjson)
    npieces = len(pd['pieces'])

# create universe set
U = enumpiece([[0, 0, 0]], pd['dimensions'], fcoord=lambda x: point2int(x, pd['dimensions']), fpiece=lambda x: x[0]) + list(range(-1, -1-npieces, -1))

# create piece enumerations
pieces = []
for i, piece in enumerate(pd['pieces']):
    if i == 0:
        pieces += enumpiece(piece, pd['dimensions'], rotate=False, fcoord=lambda x: point2int(x, pd['dimensions']), fpiece=lambda x: x + [-i-1])
    pieces += enumpiece(piece, pd['dimensions'], fcoord=lambda x: point2int(x, pd['dimensions']), fpiece=lambda x: x + [-i-1])

x = pulp.LpVariable.dicts('Choice', (list(range(len(pieces)))), 0, 1, pulp.LpInteger)

prob = pulp.LpProblem(pd['name'], pulp.LpMinimize)

prob += 0, 'Arbitrary objective'

# Constraint: require exactly 1 piece to cover each point in space
for u in U:
    prob += pulp.lpSum([x[i] if u in pieces[i] else 0 for i in range(len(pieces))]) == 1, ''

prob.writeLP("tetriscube.lp")

tetriscubeout = open('tetriscubeout.txt','w')

nsolutions = 0
maxsolutions = 10
solutions = []
alpha = 0.4
piecetocolor = [tuple(list(plt.cm.get_cmap('jet', npieces)(i)[:-1]) + [alpha]) for i in range(npieces)]

while True:
    prob.solve()
    if pulp.LpStatus[prob.status] == 'Optimal':
        solutions.append([])
        nsolutions += 1
        for i in range(len(pieces)):
            if pulp.value(x[i]) != 0:
                solutions[-1].append((i, [int2point(id, pd['dimensions']) for id in pieces[i][:-1]]))
        print("\r{} Solutions found = {}".format(datetime.now(), nsolutions), sep=' ', end='')
        colors = []
        ax=plt.figure().gca(projection='3d')
        for i, (j, piece) in enumerate(solutions[-1]):
            viewpiece(piece, pd['dimensions'], ax=ax, color=piecetocolor[-(pieces[j][-1]+1)], export='solution_{}_step_{}'.format(nsolutions, i))

        if nsolutions >= maxsolutions:
            print('\nReached target number of solutions')
            break
        
        # add constraint to avoid finding solution again
        prob += pulp.lpSum([x[i] for i, piece in solutions[-1]]) <= len(solutions[-1]) - 1, ''
    else:
        print('\nNo more solutions')
        break

tetriscubeout.close()