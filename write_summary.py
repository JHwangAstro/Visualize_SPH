

def get_header(titles, keys, run, base):

    # write header
    header = base

    for title, k in zip(titles, keys):

        # grab a sample run
        if isinstance(run[k], list):
            for i, y in enumerate(run[k]):
                header = header + '\t' + title + str(i+1)
                if len(title) < 3: header = header + '\t'

        else:
            header = header + '\t' + title
            if len(title) < 3: header = header + '\t'

    return header




def write_system_parameters(fname, runs, params):

    f = open(fname, 'w')

    titles, keys = zip(*params)

    # write header
    header = get_header(titles, keys, runs[0][1], 'run')

    f.write(header + '\n')


    # write the values
    for index, run in runs:

        f.write(str(index) + '\t')

        for k in keys:

            if isinstance(run[k], list):
                for x in run[k]:
                    f.write(str(x) + '\t')

            else:
                f.write(str(run[k]) + '\t')

        f.write('\n')

    f.close()




def write_system_parameters_formatted(fname, runs, params, base = 'run', sort_key = 'status', paper = 0):

    f = open(fname, 'w')

    titles, keys, conv, rules = zip(*params)

    header = get_header(titles, keys, runs[0][1], base)

    f.write(header + '\n')

    # write the values
    # sort by sort key (status)
    prefix = ''
    suffix = '\t'
    final_suffix = ''
    if paper:
        prefix = '$ '
        suffix = '\t$ & $\t'
        final_suffix = '\t$ & \\\\'

    # rules for writing values
    def write_value(_f, value, _conv, _rule, _suffix = suffix):

        if _rule == 'string':
            _f.write(value + _suffix)
        elif _rule == 'sci':
            _f.write(str("{0:.2E}".format(value/_conv)).replace('E-05','\\times10^{-5}').replace('E-06','\\times10^{-6}').replace('E-07','\\times10^{-7}') + _suffix)
        elif _rule == 'small':
            _f.write(str("{0:.4f}".format(value/_conv)) + _suffix)
        else:
            _f.write(str("{0:.3f}".format(value/_conv)) + _suffix)


    i = 0
    for index, run in sorted(runs, key = lambda _run: _run[1][sort_key]):

        i = i + 1
        if paper:
            f.write(prefix + str(i) + suffix)
        else:
            f.write(prefix + str(index) + suffix)

        for j, (k, c, rule) in enumerate(zip(keys, conv, rules)):

            if isinstance(run[k], list):
                if j+1 < len(keys):
                    for x in run[k]:
                        write_value(f, x, c, rule)
                else:
                    for h, x in enumerate(run[k]):
                        if h+1 < len(run[k]):
                            write_value(f, x, c, rule)
                        else:
                            write_value(f, x, c, rule, final_suffix)

            else:
                if j+1 < len(keys):
                    write_value(f, run[k], c, rule)
                else:
                    write_value(f, run[k], c, rule, final_suffix)

        f.write('\n')

    f.close()




'''
def write_system_parameters_formatted_extra(fname, runs):

    me_to_g = 5.9722*10.**27.
    ms_to_g = 1.99*10.**33.
    ms_to_me = ms_to_g/me_to_g

    keys = ['m', 'mgas', 'mc', 'rhoc', 'a', 'e', 'd', 'period_ratio', 'd_ratio', 'd_min', 'e_c', 'status']
    strings = ['status']
    not_attribute = ['mc', 'rhoc']
    convert = ['m']

    f = open(fname, 'w')
    f.write('run\tm1\t\tm2\t\tmgas1\tmgas2\tmc1\tmc2\trhoc1\trhoc2\ta1\t\ta2\t\te1\t\te2\t\trho1\trho2\tP2/P1\td2/d1\tdmin\tec\t\tstatus\n')

    for index, run in sorted(runs, key = lambda k: k[1]['status']):
        mcore = [m*ms_to_me*(1.-mgas) for m, mgas in zip(run['m'], run['mgas'])]
        f.write(str(index) + '\t')
        for k in keys:
            if k not in not_attribute:
                if isinstance(run[k], list):
                    for x in run[k]:
                        if k in strings: f.write(x + '\t')
                        if k in convert: f.write(str("{0:.3f}".format(x*ms_to_me)) + '\t')
                        if k not in strings and k not in convert: f.write(str("{0:.3f}".format(x)) + '\t')
                else:
                    if k in strings: f.write(run[k] + '\t')
                    if k not in strings: f.write(str("{0:.3f}".format(run[k])) + '\t')

            else:
                if k == 'mc':
                    for _mcore in mcore:
                        f.write(str("{0:.3f}".format(_mcore)) + '\t')

                if k == 'rhoc':
                    for _mcore in mcore:
                        rhoc = float(mass_radius_tables.GetDensity(33, _mcore))
                        f.write(str("{0:.3f}".format(rhoc)) + '\t')

        f.write('\n')
    f.close()
'''





