import timeit


def one():
    def two():
        return 1
    return two()


def three():
    return four()


def four():
    return 1


print timeit.timeit('one()', setup="from __main__ import one")
print timeit.timeit('three()', setup="from __main__ import three, four")
