from library.SM import particles

from . import eps


def oscillation_pattern_test():
    map = particles.leptons.oscillations_map()

    flavours = ('electron', 'muon', 'tau')

    for f1 in flavours:
        assert abs(sum(map[(f1, f2)] for f2 in flavours) - 1.) < eps
