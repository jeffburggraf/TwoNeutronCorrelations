from collections import OrderedDict
import numpy as np

Data = \
    OrderedDict([((24, 76), [((0.6, 0.6), (0.39, 0.06)), ((1.8, 0.6), (0.65, 0.07)), ((1.8, 1.8), (1.43, 0.08)),
                             ((3.0, 0.6), (0.97, 0.13)), ((3.0, 1.8), (1.6, 0.15)), ((3.0, 3.0), (1.82, 0.24)),
                             ((4.2, 0.6), (0.78, 0.17)), ((4.2, 1.8), (1.58, 0.21)), ((4.2, 3.0), (2.04, 0.41)),
                             ((4.2, 4.2), (2.63, 0.57)), ((5.4, 0.6), (0.64, 0.21)), ((5.4, 1.8), (1.03, 0.24)),
                             ((5.4, 3.0), (0.67, 0.32)), ((5.4, 4.2), (1.8, 0.66)), ((5.4, 5.4), (2.63, 1.05))]), (
                 (76, 128), [((0.6, 0.6), (0.4, 0.06)), ((1.8, 0.6), (0.7, 0.07)), ((1.8, 1.8), (0.91, 0.07)),
                             ((3.0, 0.6), (0.78, 0.12)), ((3.0, 1.8), (1.51, 0.15)), ((3.0, 3.0), (1.91, 0.29)),
                             ((4.2, 0.6), (0.72, 0.16)), ((4.2, 1.8), (0.91, 0.18)), ((4.2, 3.0), (0.73, 0.3)),
                             ((4.2, 4.2), (0.25, 0.31)), ((5.4, 0.6), (0.22, 0.17)), ((5.4, 1.8), (0.39, 0.18)),
                             ((5.4, 3.0), (1.0, 0.4)), ((5.4, 4.2), (1.02, 0.64)), ((5.4, 5.4), (2.74, 1.09))]), (
                 (128, 180), [((0.6, 0.6), (0.45, 0.07)), ((1.8, 0.6), (0.77, 0.08)), ((1.8, 1.8), (1.55, 0.09)),
                              ((3.0, 0.6), (1.2, 0.15)), ((3.0, 1.8), (1.83, 0.16)), ((3.0, 3.0), (2.17, 0.29)),
                              ((4.2, 0.6), (0.81, 0.18)), ((4.2, 1.8), (2.12, 0.25)), ((4.2, 3.0), (1.95, 0.39)),
                              ((4.2, 4.2), (2.47, 0.74)), ((5.4, 0.6), (0.49, 0.21)), ((5.4, 1.8), (1.42, 0.28)),
                              ((5.4, 3.0), (0.86, 0.39)), ((5.4, 4.2), (0.64, 0.47)), ((5.4, 5.4), (0.42, 0.66))])])

table_index = 0


print('[[2n corr theta VS n-n energies|go_back]]')
for theta_range,data  in Data.items():
    n_entries = len(data)
    l = int(0.5*(-1 + np.sqrt(1 + 8*(n_entries))))
    bin_width = ((6.) / l)

    if table_index == 0:
        print('\n[[File:2NCorr th vs nn {0}x{1}.png| 600px]]\n'.format(l, len(Data)))

    theta_range = '{0} - {1}'.format(*theta_range)
    table_heading = """{{| class="wikitable" style="width:30%"\
    \n|-\n! colspan = {n_entries}| theta range: {theta_range}\n|-\n""".format(n_entries=l+1, theta_range=theta_range,
                                                                 pos = ['Left','right'][table_index%2])

    table = table_heading

    bins = list(sorted(set([d[0][0] for d in data])))

    I = 0

    table += '!Energy range!!' + ' !! '.join(['{0:.1f} - {1:.1f}'.format(b-bin_width/2., b+ bin_width/2.) for b in bins]) + '\n' + '|-\n'
    # table += '! rowspan = {0} !! Energy !!'.format(l+1) + ' !! '.join(['{:.2f}'.format(b) for b in bins]) + '\n' + '|-\n'
    for i in range(l):
        row = []
        for j in range(l):
            if j>i:
                row += [' -- ']
                continue
            else:
                value = data[I][1][0]
                err = data[I][1][1]
                row += ['{0:.2f} +/- {1:.2f}'.format(value,err)]
                I += 1
        b = bins[i]
        table += '!{0:.1f} - {1:.1f}\n|'.format(b-bin_width/2., b+ bin_width/2.) + '||'.join(row) + '\n|-\n'

    table_index += 1


    table += "|}"

    print (table)

print('\n\n=={0} energy bins, {1} theta bins=='.format(l, table_index))
