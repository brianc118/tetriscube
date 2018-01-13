"""
Solver for the tetris cube problem

"""

from pprint import pprint
from functools import reduce
import json
import pulp

probfile = 'standardcube.json'

def mmin(m,i,j):
    return [r[:j] + r[j+1:] for r in (m[:i]+m[i+1:])]

def det(m):
    if len(m) == 2:
        return m[0][0] * m[1][1] - m[0][1] * m[1][0]
    d = 0
    for c in range(len(m)):
        d += ((-1) ** c) * m[0][c] * det(mmin(m, 0, c))
    return d

def permutations(l):
    """ 
    Returns permutations of elements of list l
    
    """
    if len(l) <= 0:
        raise ValueError("List to permute is empty")
    elif len(l) == 1:
        return [[l[0]]]
    else:
        p = []
        for i in range(len(l)):
            for sp in permutations(l[:i] + l[i+1:]):
                p.append([l[i]] + sp)
        return p

def getmat(p):
    """
    Returns transformation corresponding to permutation list p

    """
    n = len(p)
    if min(p) < 0 or max(p) >= n:
        raise ValueError("Invalid permutation list")
    m = []
    for i in range(n):
        m.append([1 if j == p[i] else 0 for j in range(n)])
    return m

def orthrotations(points):
    """
    Returns a list of all orthogonal simple rotations of points.

    """
    if len(points) <= 0:
        raise ValueError("No points given")
    else:
        ans = []
        n = len(points[0])
        perms = permutations(list(range(n)))
        
        for i in range(2**n):
            bmap = [1 if i >> k & 1 else -1 for k in range(0, n)]
            # bmap = [1 if digit=='1' else -1 for digit in bin(i)[2:]]
            for o in perms:
                if det(getmat(o)) * reduce(lambda x, y: x * y, bmap) != 1:
                    continue
                base = []
                for p in points:
                    base.append([0] * n)
                    for j in range(n):
                        base[-1][j] = p[o[j]] * bmap[j]
                    # print(bmap, o, base)
                ans.append(base)
        return ans


# enumerate piece where piece is an n-orthotope
# rotations assume NO symmetry. i.e. every rotation results in unique polytope
def enumpiece(piece, dimensions, rotate=True):
    if len(piece) <= 0:
        raise ValueError("Piece has no volume")

    n = len(piece[0])  # no. of dimensions
    if len(dimensions) != n:
        raise ValueError("Bad dimensions")

    enum = []
    if rotate:
        # n-D simple rotation about plane
        for p in orthrotations(piece):
            enum.extend(enumpiece(p, dimensions, rotate=False))
        return enum
    
    pmin, pmax = [], []

    for i in range(n):
        pmin.append(min(piece, key=lambda x: x[i])[i])
        pmax.append(max(piece, key=lambda x: x[i])[i])

    mincoord, maxcoord = [0] * n, [dimensions[i] for i in range(n)]

    def nestedshift(d, clones, shift=n*[0]):
        # print(d, shift)
        if d < 0:
            raise ValueError("Invalid recursion depth d")
        elif d == 0:
            clones.append([[v for v in p] for p in piece])
            for p in clones[-1]:
                for i in range(n):
                    p[i] += shift[i]
        else:
            d -= 1
            for i in range(mincoord[d]-pmin[d], maxcoord[d]-pmax[d]):
                newshift = shift[:]
                newshift[d] = i
                nestedshift(d, clones=clones, shift=newshift)
    
    nestedshift(d=n, clones=enum)
    return enum

"""
Now we set up the problem.

"""

with open(probfile) as pjson:
    pd = json.load(pjson)
    npieces = len(pd['pieces'])

prob = pulp.LpProblem(pd['name'], LpMinimize)

# set objective function to minimise pieces used
prob += sum(ans)

# set constraints


pprint(enumpiece(pd['pieces'][0], pd['dimensions']))
# pprint(orthrotations([[1,2,3]]))
