import numpy as np
from itertools import product

class Point:
    def __init__(self, y, err):
        self.y=y
        assert err>0
        self.err=err
        self.dx = 0

    def is_conflict(self, other):
        assert isinstance(other, Point)
        if other.dx != self.dx:
            return False
        if abs(self.y - other.y)<(self.err + other.err):
            return True
        else:
            return False

def solve(_data):
    for l in _data:
        assert len(l)<len(_data)

    solution = np.zeros(len(_data))


class test_iter:
    def __init__(self, n):
        self.n = n

    def __iter__(self):
        return self

    def next(self):
        pass

def test(sol, constraints, verbosse = False):
    for index, contraint in enumerate(constraints):
        if sol[index] in map(lambda i: sol[i], contraint):
            if verbosse:
                print "failed ", index, sol[index], contraint
            return False
    return True

def gen_puzzle(n):
    result_not = [set() for i in range(n)]
    for i in range(n):

        a = range(n)
        a.remove(i)
        result_not[i] |= set(np.random.choice(a, np.random.randint(0,len(a)-1), replace=False))
        for j in result_not[i]:
            result_not[j].add(i)
    # for i, s in enumerate(result):
    #     for j in s:
    #         result[j].add(i)


    for i, s in enumerate(result_not):
        eq = set(range(n)) - s - set([i])
        print "{0}: != {1}  ".format(i, " ".join(map(str, s)), " ".join(map(str, eq)))

    best_sum = None

    best_sol = None

    for s in list(product(*([range(n-3)]*(n)))):
        if test(s, result_not):
            v = len(set(s))
            if best_sum is None:
                best_sum = v
                best_sol = s
            else:
                if v<best_sum:
                    best_sum = v
                    best_sol = s
    print "random searched best solution: ", best_sol

    return result_not


def solve(l):
    sol = range(len(l))

    # print l
    # print l

    l = [set(i) for i in l]

    # l = sorted(l, key=lambda x: -len(x))

    def replace(old, new):
        for i, _set_ in enumerate(l):
            assert isinstance(_set_, set)
            if old in _set_:
                _set_.remove(old)
                _set_.add(new)

    for i, const in enumerate(l):
        best = min(set(range(len(l))) - set(const))
        replace(i, best)
        sol[i] = best

    return sol





puzzle = gen_puzzle(8)
sol = solve(puzzle)

print "Algorithmic sol. :             ",  sol

print "test ", test(sol, puzzle, True)

