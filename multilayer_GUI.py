import numpy as np
import matplotlib.pyplot as plt

from darkMatter import darkMatter
from general.plot_statistics import *
from general.utils import set_plot_params

import PySimpleGUI as sg

session_list_column = [
    [
        sg.Text("Network layout"),
    ],
    [
        sg.Text('Layers'),
        sg.In(size=(4,1), enable_events=True, key="-LAYER-NUMBER-"),
        sg.Text('Populations'),
        sg.In(size=(4,1), enable_events=True, key="-POPULATION-NUMBER-"),
        sg.In(size=(4,1), enable_events=True, key="-POPULATION-NUMBER-"),
        # sg.FolderBrowse(),
    ],
    [
        sg.Graph(
            canvas_size=(1200, 1200),
            graph_bottom_left=(0, 0),
            graph_top_right=(100, 100),
            key="-GRAPH-",
            change_submits=True,  # mouse click events
            background_color='lightblue',
            drag_submits=True)
    ]
]

processing_column = [
    [
        sg.Button('Add layer',key="-ADD-LAYER-"),
        sg.Button('Add population',key="-ADD-POPULATION-")
    ],
    [sg.Text("In here comes some input fields and processing status updates")],
    # [sg.Text(size=(40,1), key="-TOUT-")],
]

layout = [
    [
        sg.Column(session_list_column),
        sg.VSeparator(),
        sg.Column(processing_column),
    ]
]

window = sg.Window("Data pipeline", layout)

layer_margins = (5,2)
layer_sz = (100,20)

population_margins = (2,2)
population_sz = (30,layer_sz[1]-2*population_margins[1])


parameters = {
    'L': 0,
    'P': [],
    'S': [[]],
}

def addLayer(paras):
    print('adding a new layer')
    l = paras['L']
    paras['L'] += 1
    start_point = (
        layer_margins[0],
        l*layer_sz[1]+layer_margins[1]
    )
    end_point = (
        100-layer_margins[0],
        (l+1)*layer_sz[1]-layer_margins[1]
    )
    graph.draw_rectangle(start_point, end_point,fill_color='green', line_color='red')

    paras['P'].append(0)
    addPopulation(paras,l)
    addPopulation(paras,l)


def addPopulation(paras,l):
    print(f"adding population to layer {l}")
    paras['P'][l] += 1
    p = paras['P'][l]
    start_point = (
        layer_margins[0]+(p-1)*population_sz[0]+population_margins[0],
        l*layer_sz[1]+layer_margins[1]+population_margins[1]
    )
    end_point = (
        p*population_sz[0]-population_margins[0],
        (l+1)*layer_sz[1]-layer_margins[1]-population_margins[1]
    )

    graph.draw_rectangle(start_point, end_point,fill_color='black', line_color='red')

    paras['S'].append(0)


def addSynapse(paras,l,p):
    print(f"adding population to layer {l}")
    paras['P'][l] += 1
    p = paras['P'][l]
    start_point = (
        layer_margins[0]+(p-1)*population_sz[0]+population_margins[0],
        l*layer_sz[1]+layer_margins[1]+population_margins[1]
    )
    end_point = (
        p*population_sz[0]-population_margins[0],
        (l+1)*layer_sz[1]-layer_margins[1]-population_margins[1]
    )

    graph.draw_rectangle(start_point, end_point,fill_color='black', line_color='red')

    paras['S'].append(0)



# get the graph element for ease of use later
graph = window["-GRAPH-"]

erasing = False
while True:
    event,values = window.read()

    if event=="Exit" or event==sg.WIN_CLOSED or event is None:
        break

    if event == "-GRAPH-":
        x, y = values["-GRAPH-"]
        if not erasing:
            print(f"(x,y)=({x},{y})")
            start_point = (x, y)
            # dragging = True
            drag_figures = graph.get_figures_at_location((x,y))
            print(drag_figures)
            lastxy = x, y
        else:
            end_point = (x, y)

    if event=="-ADD-LAYER-":
        addLayer(parameters)

    if event=="-ADD-POPULATION-":
        addPopulation(parameters,0)


window.close()
