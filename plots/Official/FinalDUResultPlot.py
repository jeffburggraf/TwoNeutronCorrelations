import pickle
import numpy as np
import matplotlib as mpl
from numpy.polynomial.legendre import legval

fast = False


if not fast:
    font = {'family':'sans-serif',
            'sans-serif':['Helvetica'],
            'size': 15}
    mpl.rc('font', **font)
    mpl.rc("savefig", dpi=400)

mpl.rc('text', usetex=True)

f = open('/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/plots/Official/FinalDUResultData.pickle', 'rb')

data_dict = pickle.load(f, encoding='latin1')

x_data =data_dict['x_data']
print(x_data)

def Legendre(cos_X, a, b, c):
    return a*1 + b*0.5*(3*cos_X**2 - 1) + c*0.5*(5*cos_X**3 - 3*cos_X)


def get_fit(x, params, max_n):
    x = np.linspace(min(x), max(x), 5*len(x))
    # y = [Legendre(np.cos(3.1415*xxx/180), *params) for xxx in x]
    y = [legval(np.cos(3.1415*xxx/180), params[0:max_n]) for xxx in x]
    return x, y

markers =['^','x']
e =1

colors = ['orange','blue']

axs= []
_max = 0

import matplotlib.pyplot as plt

for index, (cut, _dict) in enumerate(data_dict['cuts'].items()):
    if index == 0:
        ax = plt.subplot(2, 1, 1)
        axs.append(ax)
    elif index == 2:
        ax = plt.subplot(2, 1, 2)
        axs.append(ax)
    ax.set_xticks(np.arange(0,180+30,30))


    cut = '${0}$'.format(cut)

    y_err = _dict['y_err']
    y_data = _dict['y_data']
    params = _dict['params']
    ratio = '{0:.2f}'.format(params[2]/params[0])

    x_fit, y_fit = get_fit(x_data, params, 4)
    i = int(0.9*len(x_fit))
    # plt.text(x_fit[i], -0.15 + y_fit[i], '$P_{{0}} = {1:.2f} ; r ={0}$'.format(ratio, params[0]))
    label_info = '\\small{{$p_{{0}} = {1:.1f}, \\frac{{p_{{2}}}}{{p_{{0}}}} ={0}$}}'.format(ratio, params[0])

    label = '{0}, {1}'.format(cut, label_info)

    ax.plot(x_fit, y_fit, ls='--',c=colors[index%2])
    _x_data = np.array(x_data) + (-1)**(index)
    ax.errorbar(_x_data, y_data,yerr=y_err, label=label, marker=markers[index%2], markersize=4, fmt='o', capsize=3, elinewidth=.6,markeredgewidth=.6, c=colors[index%2])

    if max(y_data)>_max:
        _max = max(y_data)

    l = ax.legend(fontsize=13, bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode="expand", borderaxespad=0.)
    ax.set_xlabel('$\\theta_{nn}$ \small{[degrees]}')
    ax.set_ylabel('correlation \small{[arb. units]}')

    print('Number of trues for {0}: {1}'.format(cut, _dict['n_trues']))
    print('Number of accidentals for {0}: {1}'.format(cut, _dict['n_accidentals']))

    if index%2 != 0:
        for _i_, text in enumerate(l.get_texts()):
            text.set_color(colors[_i_%2])

for ax in axs:
    ax.set_ylim(0,1.1*_max)
    ax.set_yticks(np.arange(0, int(1.1*_max)+2,2))


plt.subplots_adjust(hspace=0.71)
plt.show()
