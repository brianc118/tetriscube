import math
from matrix import *
from functools import reduce

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