import math
from sys import float_info

import integral

print(integral.__doc__)


def test_distribution_interpolation():

    def f(x):
        return 1. / (math.exp(x / math.pi) + 1)

    xs = list(range(-5, 5))
    ys = list(map(f, xs))

    xxs = [x / 5. for x in range(-25, 21)]
    yys = map(f, xxs)

    assert [abs(integral.distribution_interpolation(xs, ys, x) - y) < 10 * float_info.epsilon
            for x, y in zip(xs, ys)]
    assert [abs(integral.distribution_interpolation(xs, ys, x) - y) < 10 * float_info.epsilon
            for x, y in zip(xxs, yys)]
    print("Distribution interpolation is just fine")


test_distribution_interpolation()
print(integral.D([2., 1., 1., 2.], [2, 1, 1, 2], [0, 0, 0, 0], 1, 1, [0, 1, 2, 3], [-1, -1, 1, 1]))
print(integral.Db([0., 1., 1., 0.], [0, 1, 1, 0], [0, 0, 0, 0], 1, 1, [0, 1, 2, 3], [-1, -1, 1, 1]))
