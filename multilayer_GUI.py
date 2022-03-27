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
        sg.Button('Add population',key="-ADD-POPULATION-"),
        sg.Button('Add synapse',key="-ADD-SYNAPSE-")
    ],
    [
        sg.Button('Delete stuff',key="-DELETE-ELEMENT-")
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

synapse_margins = (1,1)
synapse_sz = (10,10)



class Elements():

    graph_registry = {}

    paras = {
        'L': 0,
        'P': [],
        'S': [],
    }

    def register_graph(self,id,lvl):
        print(f"registering id {id} for level {lvl}")

    def addLayer(self):
        print('adding a new layer')
        l = self.paras['L']
        self.paras['L'] += 1
        start_point = (
            layer_margins[0],
            l*layer_sz[1]+layer_margins[1]
        )
        end_point = (
            100-layer_margins[0],
            (l+1)*layer_sz[1]-layer_margins[1]
        )
        id = graph.draw_rectangle(start_point, end_point,fill_color='green', line_color='red')

        self.register_graph(id,'L')
        self.paras['P'].append(0)
        self.paras['S'].append([])
        self.addPopulation(l)
        self.addPopulation(l)


    def addPopulation(self,l):
        print(f"adding population to layer {l}")
        p = self.paras['P'][l]
        self.paras['P'][l] += 1
        start_point = (
            layer_margins[0]+p*population_sz[0]+population_margins[0],
            l*layer_sz[1]+layer_margins[1]+population_margins[1]
        )
        end_point = (
            layer_margins[0]+(p+1)*population_sz[0]-population_margins[0],
            (l+1)*layer_sz[1]-layer_margins[1]-population_margins[1]
        )

        id = graph.draw_rectangle(start_point, end_point,fill_color='black', line_color='red')

        self.register_graph(id,'P')

        self.paras['S'][l].append(0)
        self.addSynapse(l,p)


    def addSynapse(self,l,p):
        print(f"adding synapse to population {p} in layer {l}")
        # s_idx = paras['P'][l] + p
        s = self.paras['S'][l][p]
        print('Synapse', s)
        self.paras['S'][l][p] += 1
        start_point = (
            layer_margins[0]+p*population_sz[0]+population_margins[0] + s*synapse_sz[0]+synapse_margins[0],
            l*layer_sz[1]+layer_margins[1]+population_margins[1]+synapse_margins[1]
        )
        end_point = (
            layer_margins[0]+p*population_sz[0]+population_margins[0]+(s+1)*synapse_sz[0]-synapse_margins[0],
            l*layer_sz[1]+layer_margins[1]+population_margins[1]+synapse_sz[1]-synapse_margins[1]
        )

        id = graph.draw_rectangle(start_point, end_point,fill_color='yellow', line_color='red')

        self.register_graph(id,'S')

els = Elements()

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
        els.addLayer()
    if event=="-ADD-POPULATION-":
        els.addPopulation(0)
    if event=="-ADD-SYNAPSE-":
        els.addSynapse(0,0)
    if event=="-DELETE-ELEMENT-":
        for figure in drag_figures:
            graph.delete_figure(figure)


window.close()
