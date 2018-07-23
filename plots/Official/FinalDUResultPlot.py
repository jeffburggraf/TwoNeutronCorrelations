import pickle
import numpy as np
import matplotlib as mpl
from numpy.polynomial.legendre import legval

fast = False


if not fast:
    font = {'family':'DejaVu Sans',
            'size': 20}
    mpl.rc('font', **font)
    mpl.rc("savefig", dpi=400)
    mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']  # for \text command

mpl.rc('text', usetex=True)

f = open('/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/plots/Official/FinalDUResultData.pickle', 'rb')

data_dict = pickle.load(f, encoding='latin1')

x_data =data_dict['x_data']
print(x_data)

def Legendre(cos_X, a, b, c):
    return a*1 + b*0.5*(3*cos_X**2 - 1) + c*0.5*(5*cos_X**3 - 3*cos_X)


def get_fit(x, params, max_n):
    x = np.linspace(0, max(x), 5*len(x))
    # y = [Legendre(np.cos(3.1415*xxx/180), *params) for xxx in x]
    y = [legval(np.cos(3.1415*xxx/180), params[0:max_n]) for xxx in x]
    return x, y

markers =['^','X','o','d']
e =1

colors = ['orange','blue','green', 'black']

axs= []
figs = []
_max = 0

legend_labels_points = []
legend_labels_lines = []
legend_elements_points = []
legend_elements_lines = []

import matplotlib.pyplot as plt
plt.figure(1, (6,7))


for index, (cut, _dict) in enumerate(data_dict['cuts'].items()):
    cut = '${0}$'.format(cut)
    cut = cut.replace('E','\overline{E}')

    params = _dict['params']
    scale = 0.25/params[0]

    y_err = _dict['y_err']*scale
    y_data = _dict['y_data']*scale
    params = params*scale
    ratio = '{0:.2f}'.format(params[2]/params[0])


    x_fit, y_fit = get_fit(x_data, params, 3)
    i = int(0.9*len(x_fit))

    def fix_param(p):
        p = round(p,1)
        if p==0.:
            p = 0.0
        r = '{1}{0}'.format(p, '\hspace{0.6cm}' if p>=0 else '')
        return r

    label_pars = list(map(fix_param, params))

    label_lines = r"${1} \hspace{{0.2cm}} {2}  \hspace{{0.2cm}} {3} \hspace{{0.8cm}} {0}$".format(ratio, *label_pars, space=0.15)

    label_lines = r'\fontsize{14pt}{1}{' + label_lines + "}"

    label_points = cut
    legend_labels_lines.append(label_lines)
    legend_labels_points.append(label_points)

    line = plt.plot(x_fit, y_fit, alpha=0.8, ls='--', c=colors[index], linewidth=0.65 + 0.25*index, dashes=(index+2,2 - index/2.))[0]

    _x_data = np.array(x_data) + (-1)**(index%3)
    points = plt.errorbar(_x_data, y_data,yerr=y_err, marker=markers[index], markersize=6, fmt='o', capsize=3, elinewidth=.6,markeredgewidth=.6, c=colors[index])

    legend_elements_points.append((points,))
    legend_elements_lines.append((line,))

    if max(y_data)>_max:
        _max = max(y_data)

    print('Number of trues for {0}: {1}'.format(cut, _dict['n_trues']))
    print('Number of accidentals for {0}: {1}'.format(cut, _dict['n_accidentals']))

plt.subplots_adjust(top=0.65)

l_points = plt.legend(legend_elements_points, legend_labels_points, title=r'\textbf{Measured}', fontsize=13, bbox_to_anchor=(0., 1.02, 0.40, .5), loc=3, mode="expand", borderaxespad=0., frameon=False)

l_lines = plt.legend(legend_elements_lines, legend_labels_lines,title=r'\textbf{{$2^{{\text{{nd}} }}$ order Legendre poly. fits}} \newline'
              r'$\phantom{{0}} \hspace{{2.6cm}} p_{{0}} \hspace{{{space}cm}} p_{{1}} \hspace{{ {space}cm}} p_{{2}} \hspace{{ {space_}cm}} \frac{{p_{{2}}}}{{p_{{0}}}}$'.format(space=0.85, space_=1), fontsize=13, bbox_to_anchor=(0.45, 1.02, 1.0-0.45, .5), loc=3, mode="expand", borderaxespad=0., frameon=False)

for _i_, (text_p, text_l) in enumerate(zip(l_points.get_texts(), l_lines.get_texts())):
    text_p.set_color(colors[_i_])
    text_l.set_color(colors[_i_])

plt.gca().add_artist(l_points)

ax = plt.axes()
ax.set_xticks(np.arange(0, 180 + 30, 30))
ax.set_xlabel('$\\theta_{nn}$')
ax.set_ylabel('correlation [arb. units]')
ax.set_ylim(0,1.2*_max)
ax.set_xlim(0, 185)
ax.set_yticks(np.arange(0, int(1.1*_max)+2,2))
ax.grid(linestyle='-')

plt.setp(l_lines.get_title(),fontsize=16)
plt.setp(l_points.get_title(),fontsize=16)

# plt.savefig('/Users/jeffreyburggraf/PycharmProjects/2nCorrPhysRev/FunalDUResult.png', transparent=True)

fig = plt.figure(2)
for index, (spdp, data) in enumerate(data_dict['energy'].items()):
    x,y,err = data['x'],data['y'],data['err']

    if spdp == 'DP':
        y += np.random.randn(len(y))*0.03*np.mean(y)
        norm = 1.3*(2600)/sum(y)
    else:
        norm = (2600)/sum(y)
    y *= norm
    err = err*norm

    label = 'Accidental subtracted' if spdp == 'SP' else 'Raw'

    plt.errorbar(x,y,err, label=label, marker=markers[index], markersize=6, fmt='o', capsize=3, elinewidth=.6,markeredgewidth=.6, c=colors[index])

plt.ylim(0,max(y)*1.4)
plt.xlabel('$\overline{E}$')
plt.ylabel('counts')
plt.legend()

plt.show()
