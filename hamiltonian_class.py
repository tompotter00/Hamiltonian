import matplotlib.pyplot as plt
import matplotlib


class Hamilton:  # 哈密顿算符中的项
    def __init__(self, hami_c, hami_c_sgn, hami_h):
        self.c = hami_c
        self.c_sgn = hami_c_sgn
        self.h = hami_h

    # 增加了长度和迭代器方法
    def __len__(self):
        if len(self.c_sgn) != len(self.c):
            raise ValueError('Hamilton must have same number of elements')
        return len(self.c)

    def dim(self):
        return len(self.c[0])

    def __iter__(self):
        if len(self.c_sgn) != len(self.c):
            raise ValueError('Hamilton must have same number of elements')
        return iter([self.c[i], self.c_sgn[i]] for i in range(self.__len__()))


class Conf:
    def __init__(self, hamilton_operator, lattice_length, is_spin=False):
        self.hamilton_operator = hamilton_operator
        self.lattice_length = lattice_length
        self.is_spin = is_spin

    def dim(self):
        if self.lattice_length is None:
            return 0
        elif isinstance(self.lattice_length, int):
            return 1
        else:
            return len(self.lattice_length)


def hamiltonian2conf(hamilton_operator, is_spin=False):
    op = [[i.c, i.c_sgn, i.h] for i in hamilton_operator]
    return Conf(op, None, is_spin)


def multi_range(lattice_length: list, is_spin=False):
    if is_spin:
        return [[0] + i for i in multi_range(lattice_length)]
    if len(lattice_length) == 1:
        return [[i] for i in range(lattice_length[0])]
    tp = multi_range(lattice_length[1:])
    ans = []
    for i in range(lattice_length[0]):
        for j in tp:
            ans.append([i] + j)
    return ans


def init_hamilton(conf):
    hamilton_operator = [Hamilton(c, c_sgn, h) for c, c_sgn, h in conf.hamilton_operator]
    total_hamilton_operator = []
    for hamilton_item in hamilton_operator:
        if conf.dim() == 0:
            total_hamilton_operator.append(hamilton_item)
        elif conf.dim() == 1:
            for j in range(conf.lattice_length):
                if conf.is_spin:
                    j_spin = [0, j]
                    total_hamilton_operator.append(Hamilton(
                        [[(j_spin[i] + k[i]) % conf.lattice_length for i in range(conf.dim() + 1)] for k in
                         hamilton_item.c],
                        hamilton_item.c_sgn, hamilton_item.h))
                else:
                    j_spin = j
                    total_hamilton_operator.append(
                        Hamilton([(j_spin + k) % conf.lattice_length for k in hamilton_item.c], hamilton_item.c_sgn,
                                 hamilton_item.h))
        else:
            if conf.is_spin:
                lattice_length_spin = [2] + conf.lattice_length
            else:
                lattice_length_spin = conf.lattice_length
            for j in multi_range(conf.lattice_length, conf.is_spin):
                total_hamilton_operator.append(Hamilton(
                    [[(j[i] + k[i]) % lattice_length_spin[i] for i in range(conf.dim() + (1 if conf.is_spin else 0))]
                     for
                     k in hamilton_item.c], hamilton_item.c_sgn, hamilton_item.h))
    return total_hamilton_operator


# DONE func:显示latex的Hamilton算符
def subscript(x, cnt=0, is_sum=True, is_spin=False):
    if is_spin:
        if isinstance(x, int):
            return r'\uparrow ' if x else r'\downarrow '
        else:
            return subscript(x[1:], cnt, is_sum, False) + (r'\uparrow ' if x[0] else r'\downarrow ')
    it_list = ['i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's']
    if x == 0:
        if is_sum:
            return it_list[cnt % len(it_list)]
        else:
            return '0'
    if isinstance(x, int):
        if is_sum:
            return f'i+{x}' if x > 0 else f'i{x}'
        else:
            return f'{x}'
    elif hasattr(x, '__iter__'):
        ans = ''
        for i in x:
            if ans and not is_spin:
                ans += ','
            ans += subscript(i, cnt, is_sum, False)
            cnt += 1
        return ans
    else:
        raise TypeError('x must be an integer or an iterable')


def latex_hamilton(conf, print_it=False, fig_size=None):
    # for print work properly:
    # install:
    # cm-super texlive texlive-latex-extra texlive-latex-recommended dvipng
    if fig_size is None:
        fig_size = [15, 3]
    hamilton_operator = [Hamilton(c, c_sgn, h) for c, c_sgn, h in conf.hamilton_operator]
    # hamilton_operator = init_hamilton(conf)
    hamilton_operator_latex = ''
    for hamilton_item in hamilton_operator:
        hamilton_item_latex = ''
        add_op = bool(hamilton_operator_latex)
        for i in hamilton_item:
            hamilton_item_latex += 'c{}_{{{}}}'.format(r'^\dagger' if i[1] == 1 else '',
                                                       subscript(i[0], is_sum=(conf.lattice_length is not None),
                                                                 is_spin=conf.is_spin))
        # latex of Coefficient
        if isinstance(hamilton_item.h, (int, float, complex)):
            if hamilton_item.h == 1:
                hamilton_coefficient_latex = ''
            else:
                hamilton_coefficient_latex = str(hamilton_item.h)
            if not isinstance(hamilton_item.h, complex):
                if hamilton_item.h < 0:
                    add_op = False
            else:
                hamilton_coefficient_latex = '\\underline{{{}}}'.format(hamilton_coefficient_latex)
                # underline for complex number
        elif callable(hamilton_item.h):
            if hamilton_item.h.__doc__:
                hamilton_coefficient_latex = hamilton_item.h.__doc__ + ' '
            else:
                hamilton_coefficient_latex = hamilton_item.h.__name__
        else:
            raise TypeError('h mast be number or function')
        if add_op:
            hamilton_operator_latex += '+'
        hamilton_operator_latex += hamilton_coefficient_latex
        hamilton_operator_latex += hamilton_item_latex

    # latex of sum_op
    if conf.lattice_length is None:
        sum_op = ''
    elif conf.dim() == 1:
        sum_op = r'\sum_i^{{{}}}'.format(conf.lattice_length)
    else:
        sum_op = r'\sum_{{{}}}^{{{}}}'.format(subscript([0 for _ in range(len(conf.lattice_length))]),
                                              str(conf.lattice_length)[1:-1])

    hamilton_operator_latex = sum_op + hamilton_operator_latex
    if print_it:
        matplotlib.rcParams['text.usetex'] = True
        fig = plt.figure(figsize=fig_size)
        ax = fig.add_axes((0., 0., 1., 1.))
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        plt.text(0.5, 0.5, f'${hamilton_operator_latex}$', fontsize=64,
                 verticalalignment='center', horizontalalignment='center',
                 fontdict=dict(family='monospace'))
        # plt.text(0., 0., hamilton_operator_latex)
        # fig.tight_layout()
        plt.show()
    return hamilton_operator_latex


if __name__ == '__main__':
    def beta():
        """\\beta (t)"""
        return 1


    #
    hmt = [
        [[[0, 0, 0], [0, 0, 0]], [0, 0], 1],
        [[[1, 0, 0], [0, 1, 0]], [0, 1], 2j],
        [[[0, 1, 0], [1, 2, 0]], [1, 1], beta],
    ]
    cf = Conf(hmt, [2, 2], is_spin=True)
    print(latex_hamilton(cf, print_it=True, fig_size=[20, 3]))
    print(latex_hamilton(hamiltonian2conf(init_hamilton(cf), is_spin=cf.is_spin), print_it=True, fig_size=[60, 3]))

    # hmt = [
    #     [[[0, 0], [0, 0]], [0, 0], 1],
    #     [[[1, 0], [0, 1]], [0, 1], 2j],
    #     [[[0, 1], [1, 2]], [1, 1], beta],
    # ]
    # cf = Conf(hmt, 2, is_spin=True)
    # print(latex_hamilton(cf, print_it=True))
    # print(latex_hamilton(hamiltonian2conf(init_hamilton(cf), is_spin=cf.is_spin), print_it=True, fig_size=[50, 3]))

    # hmt = [
    #     [[0, 0], [0, 0], 1],
    #     [[1, 0], [0, 1], 3],
    #     # [[[0, 9], [0, 0]], [0, 1], 1],
    # ]
    # cf = Conf(hmt, None,is_spin=True)
    # print(latex_hamilton(cf, print_it=True))
    # print(latex_hamilton(hamiltonian2conf(init_hamilton(cf),is_spin=cf.is_spin), print_it=True, fig_size=[50, 3]))
