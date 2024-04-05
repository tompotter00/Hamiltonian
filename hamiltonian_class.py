class Hamilton:  # 哈密顿算符中的项
    def __init__(self, hami_c, hami_c_sgn, hami_h):
        self.c = hami_c
        self.c_sgn = hami_c_sgn
        self.h = hami_h


class Conf:
    def __init__(self, hamilton_operator, lattice_length):
        self.hamilton_operator = hamilton_operator
        self.lattice_length = lattice_length


def init_hamilton(conf):
    hamilton_operator = [Hamilton(c, c_sgn, h) for c, c_sgn, h in conf.hamilton_operator]
    total_hamilton_operator = []
    for hamilton_item in hamilton_operator:
        for j in range(conf.M):
            total_hamilton_operator.append(Hamilton([(j + k) % conf.M for k in hamilton_item.c],
                                                    hamilton_item.c_sgn,
                                                    hamilton_item.h))
    return total_hamilton_operator
