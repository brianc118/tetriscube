"""
Solver for the tetris cube problem

"""

from datetime import datetime
from pprint import pprint
from cubeviewer import *
from matrix import *
from tetrismath import *
import json
import pulp

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

probfile = 'standardcube.json'
answerfile = 'tetriscubeout.txt'
maxsolutions = 10
alpha = 0.4  # alpha for output

def enumpiece(piece, dimensions, rotate=True, fcoord=lambda x: x, fpiece=lambda x: x):
    """ enumerate piece where piece is an n-orthotope. """
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
U = enumpiece([[0, 0, 0]], pd['dimensions'], fcoord=lambda x: point2int(
    x, pd['dimensions']), fpiece=lambda x: x[0]) + list(range(-1, -1 - npieces, -1))

# create piece enumerations
pieces = []
for i, piece in enumerate(pd['pieces']):
    if i == 0:
        pieces += enumpiece(piece, pd['dimensions'], rotate=False, fcoord=lambda x: point2int(
            x, pd['dimensions']), fpiece=lambda x: x + [-i - 1])
    pieces += enumpiece(piece, pd['dimensions'], fcoord=lambda x: point2int(
        x, pd['dimensions']), fpiece=lambda x: x + [-i - 1])

# create problem variables
x = pulp.LpVariable.dicts(
    'Choice', (list(range(len(pieces)))), 0, 1, pulp.LpInteger)

prob = pulp.LpProblem(pd['name'], pulp.LpMinimize)

prob += 0, 'Arbitrary objective'

# Constraint: require exactly 1 piece to cover each point in space
for u in U:
    prob += pulp.lpSum([x[i] if u in pieces[i] else 0 for i in range(len(pieces))]) == 1, ''

prob.writeLP("tetriscube.lp")

tetriscubeout = open(answerfile, 'w')

nsolutions = 0
solutions = []
piecetocolor = [tuple(list(plt.cm.get_cmap('jet', npieces)(i)[:-1]) + [alpha]) for i in range(npieces)]

while True:
    prob.solve()
    if pulp.LpStatus[prob.status] == 'Optimal':
        solutions.append([])
        nsolutions += 1
        for i in range(len(pieces)):
            if pulp.value(x[i]) != 0:
                solutions[-1].append((i, [int2point(id, pd['dimensions']) for id in pieces[i][:-1]]))
        print('\r{} Solutions found = {}'.format(
            datetime.now(), nsolutions), sep=' ', end='')
        colors = []
        ax = plt.figure().gca(projection='3d')
        for i, (j, piece) in enumerate(solutions[-1]):
            viewpiece(piece, pd['dimensions'], ax=ax, color=piecetocolor[
                      -(pieces[j][-1] + 1)], export='solution_{}_step_{}'.format(nsolutions, i))

        if nsolutions >= maxsolutions:
            print('\nReached target number of solutions')
            break

        # add constraint to avoid finding solution again
        prob += pulp.lpSum([x[i] for i, piece in solutions[-1]]) <= len(solutions[-1]) - 1, ''
    else:
        print('\nNo more solutions')
        break

tetriscubeout.close()
