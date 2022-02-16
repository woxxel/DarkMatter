import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import math
import pandas as pd
import numpy as np
from numpy.ma import masked_array
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from pythonCode.network import network

from darkMatter import darkMatter

# def plot_interactive(
def plot_fins(x_arr,y_arr,gamma,chi,regions):

    levs = range(100)
    plt_para = {
        'cb_plotted': False,
        'bnw': mcolors.LinearSegmentedColormap.from_list(name='black_n_white',colors=[(0,(0,0,0)),(1,(1,1,1))],N=len(levs)-1),
        'heat': mcolors.LinearSegmentedColormap.from_list(name='heat',colors=[(0,'y'),(0.5,'r'),(1,(0.1,0.1,0.5))],N=len(levs)-1),
        'bnw_regions': mcolors.LinearSegmentedColormap.from_list(name='black_n_white_regions',colors=[(0,(1,1,1)),(1,(0.8,0.8,0.8))],N=3),
        'ax_label': [],
        'const_label': []
    }

    # set_labels(results[i],para_order,plt_para)

    mask_inconsistent = (regions == 3)
    mask_no_peak = (regions == 2)
    mask_implausible = (regions == 1)

    mask_dark_matter = (gamma**2 < 1)

    plot_gamma = masked_array(gamma**2,mask_inconsistent + mask_no_peak)#masked_array(gamma,gamma>0)
    plot_chi = masked_array(chi,mask_inconsistent + mask_no_peak + mask_dark_matter)
    plot_regions = masked_array(regions,np.invert(mask_inconsistent + mask_no_peak))
    plot_implausible = masked_array(regions,np.invert(mask_implausible))

    plt_para['heat'].set_bad('k',0.)
    norm = mcolors.Normalize(vmin=0, vmax=5, clip=False)
    mapper = cm.ScalarMappable(norm=norm, cmap=plt_para['heat'])
    plot_chi = mapper.to_rgba(plot_chi)

    plt_para['bnw'].set_bad('k',0.)
    norm = mcolors.Normalize(vmin=0, vmax=1, clip=False)
    mapper = cm.ScalarMappable(norm=norm, cmap=plt_para['bnw'])
    plot_gamma = mapper.to_rgba(plot_gamma)

    plot_regions = mapper.to_rgba(1-plot_regions/3.)

    hovertemplate = "" + \
    "<b>alpha: %{customdata[0]}</b> " + \
    "<b>nu: %{customdata[1]}</b>"

    fig = go.Figure()
    fig.add_traces(go.Image(z=plot_gamma, colormodel='rgba',zmax=[2,2,2,1]))
    fig.add_traces(go.Image(z=plot_chi, colormodel='rgba',zmax=[1,1,1,1]))
    fig.add_traces(go.Image(z=plot_regions, colormodel='rgba',zmin=[-4,-4.,-4.,0],zmax=[0.5,0.5,0.5,1],customdata=[results['alpha_0'],results['rateWnt']]))
    # fig.add_traces(go.Image(z=plot_implausible, colormodel='rgba',zmax=[1,1,1,1]))
    fig.update_yaxes(autorange=True)

    # fig.update_traces(hovertemplate = hovertemplate)
    return fig

steps=100
rateWnt=[0,20]
alpha_0=[0,0.2]
tau_G=[0.005]
eps=[0.]
eta=[0.]
n=[0]
J=-1.
Npop=1
drive=0
save=0
file_format='png'
rerun=False
compile=False

    ## load data
options = {
    'Npop': Npop,
    'order': ['rateWnt','alpha_0','tau_G','n','eta','eps'],
    'rateWnt': rateWnt,
    'alpha_0': alpha_0,
    'tau_G': tau_G,
    'eps': eps,
    'eta': eta,
    'n': n,
    'drive': drive,
    'mode_stats': 0,
    'J': J
}

results = darkMatter(steps=steps,options=options,rerun=rerun,compile=compile)

# plt.rcParams['font.family'] = ['Tahoma','Verdana']
# plt.rcParams['font.size'] = 12
# plt.rcParams['xtick.labelsize'] = 12
# plt.rcParams['ytick.labelsize'] = 12

fig = plot_fins(results[options['order'][0]],results[options['order'][1]],results['gamma'][0,...],results['chi'][0,...],results['regions'][0,...])


print('now plotting')
app = dash.Dash(__name__)

app.layout = html.Div(children=[
    html.Div(id='graph-container',children=[
        dcc.Graph(
            id='sharkfin',
            figure=fig,
            style={'height':'600px','width': '60%'}
        ),
        html.Div(
            id='graphs',
            style={'width': '40%'},
            children=[
                dcc.Graph(
                    id='distribution',
                ),
                html.Div(
                    children=[
                        dcc.Graph(
                            id='gamma_x',
                            style={'height':'400px','width': '50%'}
                        ),
                        dcc.Graph(
                            id='gamma_y',
                            style={'height':'400px','width': '50%'}
                        )
                    ],
                    style={'display':'flex'}
                )
            ]
        )],
        style={'width':'95%','display':'flex','justifyContent':'start'}
    )
])
print('done?')

@app.callback(
    Output('distribution', 'figure'),
    Output('gamma_x', 'figure'),
    Output('gamma_y', 'figure'),
    Input('sharkfin', 'hoverData')
)
def onHover(hoverData):
    # print(hoverData)
    x=hoverData['points'][0]['x']
    y=hoverData['points'][0]['y']
    alpha=results['alpha_0'][y]
    net = network(alpha_0=alpha)
    steps = 1001
    rate_max = net.rate_max()
    nu = np.linspace(0,rate_max,steps)
    # nu = np.linspace(0,1,steps)
    # nu = np.logspace(-4,1.5,steps)
    rate_ratio = nu/rate_max
    gamma = results['gamma'][0,y,x]
    chi = results['chi'][0,y,x]
    delta = net.delta(nu,results['q'][0,y,x])

    distr = gamma/(rate_max*np.sqrt(-math.pi*np.log(rate_ratio)))*np.exp(-(delta**2) / 2)*rate_ratio**(gamma**2-1)*np.cosh(gamma*delta*np.sqrt(-2*np.log(rate_ratio)))
    # print(distr)
    fig = px.line(x=nu,y=distr)
    fig.update_layout(yaxis_range=[0,1],xaxis_range=[0,rate_max])


    fig_gamma_x = px.line(x=results['rateWnt'],y=results['gamma'][0,y,:])
    fig_gamma_x.update_layout(yaxis_range=[0,3],xaxis_range=[0,rate_max])

    fig_gamma_y = px.line(x=results['alpha_0'],y=results['gamma'][0,:,x])
    fig_gamma_y.update_layout(yaxis_range=[0,3],xaxis_range=[0,0.2])

    return fig,fig_gamma_x,fig_gamma_y

# @app.callback(
#     Output('distribution', 'figure'),
#     Input('sharkfin', 'hoverData'),
#     Input('sharkfin', 'clickData')
# )
# def onHover(hoverData,clickData):
#     print('in')
#     print(hoverData)
#     print(clickData)
#
#     # fig = px.Plot([1,1],[0,1])
#     return fig

if __name__ == '__main__':
    app.run_server(debug=True)


def set_labels(results,para_order,plt_para):

    for key in para_order:

        if key=='eps':
            plt_para['ax_label'].append(r'$\displaystyle \varepsilon$')
        elif key == 'eta':
            plt_para['ax_label'].append(r'$\displaystyle \eta$')
        elif key == 'n':
            plt_para['ax_label'].append(r'$\displaystyle b$')
        elif key == 'alpha_0':
            plt_para['ax_label'].append(r'$\displaystyle \alpha_0$')
        elif key == 'tau_G':
            plt_para['ax_label'].append(r'$\displaystyle \tau_G\,[ms]$')
        elif key == 'rateWnt':
            plt_para['ax_label'].append(r'$\displaystyle \bar{\nu}\,[Hz]$')

        if len(results[key]) == 1:
            plt_para['const_label'].append(plt_para['ax_label'][-1][:-1] + ' = %g$'%results[key][0])
            break
