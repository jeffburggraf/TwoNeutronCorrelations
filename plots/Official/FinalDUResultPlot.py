import pickle
import numpy as np
import matplotlib as mpl
from numpy.polynomial.legendre import legval
from scipy import optimize
from numpy.polynomial.legendre import legfit
fast = False
from lmfit import Model
from lmfit import Parameters
import re


if not fast:
    font = {'family':'DejaVu Sans',
            'size': 20}
    mpl.rc('font', **font)
    mpl.rc("savefig", dpi=400)
    mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']  # for \text command

mpl.rc('text', usetex=True)

f = open('/Users/jeffreyburggraf/PycharmProjects/TwoNeutronCorrelations/plots/Official/FinalDUResultData.pickle', 'rb')

data_dict = pickle.load(f, encoding='latin1')


def Legendre(x, p0, p1, p2):
    ps = [p0,p1,p2]
    r =legval(np.cos(3.1415*x/180), ps)
    return r


def get_fit(x, params):
    x = np.linspace(0, max(x), 5*len(x))
    y = Legendre(x,*params)
    assert x.shape == y.shape, '{0}, {1}, {2}'.format(x.shape, y.shape, y)
    return x, y

markers =['^','X','o','d']
e =1
line_styles = ['--','-.','--','-.']
line_styles = ['-', '--', '-.', ':']
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

##########################
Scale_to_true = False
######################

x_data =data_dict['x_data']

x_err = [data_dict['theta_bin_width']/2.]*len(x_data)


ParamErrors = []
for index, (cut, _dict) in enumerate(data_dict['cuts'].items()):
    if index ==0:
        init_norm = _dict['params'][0]/_dict['n_trues']

    if Scale_to_true:
        scale = init_norm*_dict['n_trues']
    else:
        scale = 1

    cut = '${0}$'.format(cut)
    cut = cut.replace('E','\overline{E}')
    y_err = _dict['y_err'] *scale
    y_data = _dict['y_data'] *scale


    model = Model(Legendre)

    init_params = legfit(np.cos(3.1415*x_data/180), y_data, 2,w = 1.0/y_err)



    params = Parameters()
    params.add('p0', init_params[0], min=0, max=2)
    params.add('p1', init_params[1], min=-2, max=2)
    params.add('p2', init_params[2], min=-2, max=2)
    result = model.fit(y_data, params, x=x_data, weights=1.0/y_err)

    _m = re.findall(r'p[0-9]: .+ \+/- ([0-9]+\.[0-9]+)', result.fit_report())
    par_errors = (list(map(float,_m)))

    params = init_params
    print('params',params)
    print('par errors', par_errors)

    x_fit, y_fit = get_fit(x_data,params)
    i = int(0.9*len(x_fit))

    ratio = '{0:.2f}'.format(params[2] / params[0])

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

    line = plt.plot(x_fit, y_fit, alpha=0.8, linestyle=line_styles[index], c=colors[index], linewidth=0.65 + 0.25*index)[0] #dashes=(index+2,2 - index/2.))[0]

    _x_data = np.array(x_data) + 0.7*(-1)**(index%2)
    points = plt.errorbar(_x_data, y_data,yerr=y_err,xerr=x_err, marker=markers[index], markersize=8, fmt='o', capsize=3, elinewidth=.6,markeredgewidth=.6, c=colors[index])

    legend_elements_points.append((points,))
    legend_elements_lines.append((line,))

    if max(y_fit)>_max:
        _max = max(y_fit)
    ratio_err = np.sqrt(((par_errors[0]*params[2])**2 + (params[0]*par_errors[2])**2)/params[0]**4)

    errs= []
    for err in list(par_errors) + [ratio_err]:
        errs.append(float('{:.2f}'.format(err)))

    ParamErrors.append(errs)

    print('Number of trues for {0}: {1}'.format(cut, _dict['n_trues']))
    print('Number of accidentals for {0}: {1}'.format(cut, _dict['n_accidentals']))

print('Parameter errors:\n{}'.format(np.array(ParamErrors)))

plt.subplots_adjust(top=0.65, left=0.16)

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
ax.set_ylabel(r'$(n\text{-}n_{\text{corr.}})/(n\text{-}n_{\text{uncorr.}})$')
ax.set_ylim(0,1.2*_max)
ax.set_xlim(0, 185)
ax.set_yticks(np.linspace(0, int(1.2*_max),4))
ax.grid(linestyle='-')

plt.setp(l_lines.get_title(),fontsize=16)
plt.setp(l_points.get_title(),fontsize=16)

fig = plt.figure(2)

DP_data = data_dict['energy']['DP']
SP_data = data_dict['energy']['SP']
for index, (spdp, data) in enumerate({"DP":DP_data, "SP":SP_data}.items()):
    x,y,err = data['x'],data['y'],data['err']

    if spdp == 'DP':
        y += np.random.randn(len(y))*0.03*np.mean(y)
        norm = 1.3*(2600)/sum(y)
        at_most = max(y*norm)
        print(at_most)
    else:
        norm = 0.6*at_most/max(y)
        print(at_most)
    y *= norm
    err = err*norm

    label = 'Accidental subtracted' if spdp == 'SP' else 'Raw'

    plt.errorbar(x,y,err, label=label, marker=markers[index], markersize=6, fmt='o', capsize=3, elinewidth=.6,markeredgewidth=.6, c=colors[index])

plt.ylim(0,max(y)*1.4)
plt.xlabel('$\overline{E}$')
plt.ylabel('counts')
plt.legend()

plt.show()
