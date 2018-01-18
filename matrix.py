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