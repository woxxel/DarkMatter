import dash
from dash import html

from scipy.stats import *
import numpy as np

import sys
from pathlib import Path
# root_dir = os.path.dirname(os.path.abspath(''))

root_dir = Path(__file__).parent.parent
if not root_dir in sys.path: sys.path.insert(0,str(root_dir))

# import DM_theory
from inference.transform_meta_to_bio import *

from dash import dcc
from dash_extensions import EventListener
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go


# Initialize the Dash app
app = dash.Dash(__name__)
app.config.suppress_callback_exceptions = True

fig_1population_surface = go.Figure()
fig_1population_linear = go.Figure()

steps = 101
gamma_arr = np.linspace(1,4,steps)
delta_arr = np.linspace(0,10,steps)

gamma,delta = np.meshgrid(gamma_arr,delta_arr)
# nu_max = 25.


# Add a slider to control the value of nu_max
app.layout = html.Div(children=[
    html.H1(children='Hello Dash'),
    html.Div(children='''Dash: A web application framework for Python.'''),
    html.Div(children=[
        dcc.Graph(id='graph_1population_surface', style={'width': '50%', 'height': '800px'},figure=fig_1population_surface),
        dcc.Graph(id='graph_1population_linear',figure=fig_1population_linear)#, style={'width': '50%', 'height': '800px'}),
    ], style={'display': 'flex', 'flex-direction': 'row'}),
    dcc.Slider(
        id='nu_max_slider',
        min=10.,
        max=100.,
        step=1.,
        value=25.,
    )
])
nu_max = 25.
nu_bar = get_nu_bar(gamma, delta, nu_max=nu_max)
alpha_0 = get_alpha_0(gamma, delta, nu_max=nu_max, nP=1)[..., 0]

fig_1population_surface.add_trace(go.Surface(
    z=nu_bar,
    x=gamma_arr,
    y=delta_arr,
    colorscale=[[0, "grey"], [1, "grey"]],
    opacity=0.4,
))
fig_1population_surface.add_trace(go.Surface(
    z=np.where(np.isfinite(alpha_0), nu_bar, np.nan),
    x=gamma_arr,
    y=delta_arr,
    surfacecolor=alpha_0,
    colorscale='fall',
    colorbar=dict(title='Alpha 0 (Finite)')
))


# fig_1population_linear.
nu_max_array = np.linspace(10, 100, steps)

fig_1population_linear.add_trace(go.Scatter(
    x=get_tau_I(nu_max_array),
    y=nu_max_array,
    mode='lines+markers',
    name='Tau_i vs Nu_max'
))

# Add a marker at the current value of nu_max and tau_I
# current_tau_I = get_tau_I(nu_max)

fig_1population_linear.add_trace(go.Scatter(
    x=[get_tau_I(nu_max)],
    y=[nu_max],
    mode='markers',
    marker=dict(color='red', size=10),
    name='Current Nu_max'
))


fig_1population_surface.update_layout(
    title='Solutions to the 1 population model',
    scene=dict(
        xaxis=dict(title='Gamma Axis'),
        yaxis=dict(title='Delta Axis'),
        zaxis=dict(title='Nu Bar')
    )
)


@app.callback(
    Output('graph_1population_surface', 'figure'),
    Input('nu_max_slider', 'value'),
    State('graph_1population_surface', 'figure'),
    prevent_initial_call=True
)
def update_surface(nu_max, fig):

    nu_bar = get_nu_bar(gamma, delta, nu_max=nu_max)
    alpha_0 = get_alpha_0(gamma, delta, nu_max=nu_max, nP=1)[..., 0]

    updated_fig = go.Figure(fig)
    # print(updated_fig['data'][1])
    updated_fig['data'][0].z = nu_bar
    # nu_new = np.where(np.isfinite(alpha_0), nu_bar, np.nan)
    # print('bla')
    # print(nu_new)
    updated_fig['data'][1].z = np.where(np.isfinite(alpha_0), nu_bar, np.nan)
    updated_fig['data'][1]['surfacecolor'] = alpha_0

    z_max = max(50, 1.1 * np.nanmax(nu_bar))
    updated_fig.update_layout(
        scene=dict(
            zaxis=dict(range=[0, z_max])
        )
    )

    return updated_fig


@app.callback(
    Output('graph_1population_linear', 'figure'),
    Input('nu_max_slider', 'value'),
    State('graph_1population_linear', 'figure'),
    prevent_initial_call=True
)
def update_linear(nu_max,fig):

    updated_fig = go.Figure(fig)
    # Update the linear plot with the current value of tau_I and nu_max
    tau_I = get_tau_I(nu_max)
    updated_fig['data'][1].x = [tau_I]
    updated_fig['data'][1].y = [nu_max]

    return updated_fig



fig_1population_surface.update_layout(
    title='Solutions to the 1 population model',
    xaxis=dict(title='Gamma Axis'),
    yaxis=dict(title='Delta Axis'),
    xaxis2=dict(title='Gamma Axis'),
    yaxis2=dict(title='Delta Axis'),
    grid=dict(rows=1, columns=1, pattern='independent')
)



# Run the app
if __name__ == '__main__':
    app.run_server(debug=True)