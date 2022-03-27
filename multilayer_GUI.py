import numpy as np
import matplotlib.colors as mcol
import matplotlib.pyplot as plt
import pandas as pd

from darkMatter import darkMatter
from general.plot_statistics import *
from general.utils import set_plot_params

import PySimpleGUI as sg

L=1
J_l = np.ones((L,L))
np.fill_diagonal(J_l,0)
print(J_l)

default_options = {
    'L': {
        'eps': 1./np.sqrt(2),
        'eta': 0.9,
        'J0_l': J_l,
        'kappa': 1.,
    },
    'P': {
        # population level parameters
        'I_ext': 1,
        'rateWnt': 1.,
        'alpha_0': 0.02,
        'tau_M': 0.01,
        'J0': 1.,
    },
    'S': {
        # psp level parameters
        'tau_I': [0.01,0.005,0.2],
        'tau_n': 0.,
        'tau_norm': 1.,
    },
}

naming = {
    'L': 'layer',
    'P': 'population',
    'S': 'synapse'
}

options = {

    'L': {},
    'P': {},
    'S': {},

    # 'order': ['tau_G','n','alpha_0','rateWnt','eta','eps'],
    'mode_stats': 0,
    'mode_calc': 0,
    'simulation': {
        # for each iteration parameter, specify (layer,population,psp)-tuple
        # specify -1 if a level of hierarchy is non-applicable
        # specify 'None' if should be applied to all candidates
        'rateWnt': [0.,20.],
        'alpha_0': [0.,0.2],

        'sim_prim': [0,-1,0],       # when population parameters are iterated, specify population number(s) (empty = all)
        'sim_sec': [0,-1,0],     # when synaptic timeconstants are iterated, specify number within population
    }
}

session_list_column = [
    [
        sg.Text("Network layout"),
    ],
    [
        sg.Text('Layers'),
        sg.In(size=(4,1), enable_events=True, key="-LAYER-NUMBER-"),
        sg.Text('Populations'),
        sg.In(size=(4,1), enable_events=True, key="-POPULATION-NUMBER-"),
        # sg.In(size=(4,1), enable_events=True, key="-POPULATION-NUMBER-"),
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

para_elements = {}
for lvl in ['L','P','S']:
    para_elements[lvl] = []
    for key in default_options[lvl]:
        print(key)
        para_elements[lvl].extend([
            sg.Text(key),
            sg.In(size=(4,1),enable_events=True, key=f'-PARA-{key}-'),
        ])

processing_column = [
    [
        sg.Button('Add layer',key="-ADD-L-"),
        sg.Button('Add population',key="-ADD-P-"),
        sg.Button('Add synapse',key="-ADD-S-")
    ],
    [
        sg.Button('Delete stuff',key="-DELETE-ELEMENT-")
    ],
]

for lvl in ['L','P','S']:
    # options = [
    #     [
    #         sg.Text(f'Some elements to change {naming[lvl]} parameters'),
    #     ],
    #     para_elements[lvl]
    # ]
    processing_column.extend([
        [sg.HSeparator()],
        [
            sg.Column([
                [sg.Text(f'Some elements to change {naming[lvl]} parameters')],
                para_elements[lvl]
                ],
                key=f'-{lvl}-OPTIONS-',visible=False)
        ],
    ])

layout = [
    [
        sg.Column(session_list_column),
        sg.VSeparator(),
        sg.Column(processing_column,key='-OPTIONS-'),
    ]
]

window = sg.Window("Data pipeline", layout)

lvl_paras = {
    'L': {
        'margins': (5,2),
        'size': (100,20),
        'col': (mcol.to_hex([0,1,0]),mcol.to_hex([0,0.8,0]))
    },
    'P': {
        'margins': (2,2),
        'size': (30,16),
        'col': (mcol.to_hex([0,0,0]),mcol.to_hex([0.2,0.2,0.2]))
    },
    'S': {
        'margins': (1,1),
        'size': (10,10),
        'col': (mcol.to_hex([0,0,1]),mcol.to_hex([0,0,0.8]))
    },
}



class Elements():

    register = pd.DataFrame(index=['graph_id'],columns=['lvl','L','P','S'])

    paras = {
        'L': 0,
        'P': [],
        'S': [],
    }

    picked = {}

    def __init__(self,window,graph):
        self.graph = graph
        self.window = window

    def register_graph(self,f_id,ids,lvl):
        print(f"registering id {f_id} for level {lvl}")
        self.register.loc[f_id] = [lvl, ids[0], ids[1] if len(ids)>1 else None, ids[2] if len(ids)>2 else None]

        # self.register[lvl]

    def pickObject(self,f_id):

        for lvl in ['L','P','S']:
            window[f'-{lvl}-OPTIONS-'].update(visible=False)

        if bool(self.picked):
            print(self.picked)
            if 'L' in self.picked.keys():
                graph_id = np.where((self.register['lvl']=='L') & (self.register['L']==self.picked['L']))[0]
                # print('layer id:',graph_id)
                if len(graph_id):
                    self.graph.Widget.itemconfig(graph_id[0], fill=lvl_paras['L']['col'][0])

            if 'P' in self.picked.keys():
                graph_id = np.where((self.register['lvl']=='P') & (self.register['L']==self.picked['L']) & (self.register['P']==self.picked['P']))[0]
                # print('population id:',graph_id)
                if len(graph_id):
                    self.graph.Widget.itemconfig(graph_id[0], fill=lvl_paras['P']['col'][0])

            if 'S' in self.picked.keys():
                graph_id = np.where((self.register['lvl']=='S') & (self.register['L']==self.picked['L']) & (self.register['P']==self.picked['P']) & (self.register['S']==self.picked['S']))[0]
                # print('synapse id:',graph_id)
                if len(graph_id):
                    self.graph.Widget.itemconfig(graph_id[0], fill=lvl_paras['S']['col'][0])

        for id in f_id:
            lvl = self.register.loc[id]['lvl']
            self.picked[lvl] = self.register.loc[id][lvl]

            window[f'-{lvl}-OPTIONS-'].update(visible=True)

            self.graph.Widget.itemconfig(id, fill=lvl_paras[lvl]['col'][1])




    # def delete


    def addLayer(self):
        # print('adding a new layer')
        l = self.paras['L']
        if l>=4:
            print("No adding of further layers possible")
            return

        self.paras['L'] += 1
        start_point = (
            lvl_paras['L']['margins'][0],
            l*lvl_paras['L']['size'][1]+lvl_paras['L']['margins'][1]
        )
        end_point = (
            start_point[0] + lvl_paras['L']['size'][0] - lvl_paras['L']['margins'][0],
            start_point[1] + lvl_paras['L']['size'][1] - lvl_paras['L']['margins'][1]
        )
        id = graph.draw_rectangle(start_point, end_point,fill_color=lvl_paras['L']['col'][0], line_color='red')

        self.register_graph(id,[l],'L')
        self.paras['P'].append(0)
        self.paras['S'].append([])
        self.picked['L'] = l
        self.addPopulation(l)
        self.addPopulation(l)


    def addPopulation(self,l):

        p = self.paras['P'][l]
        if p>=3:
            print("No adding of further populations possible")
            return
        # print(f"adding population to layer {l}")
        self.paras['P'][l] += 1
        start_point = (
            lvl_paras['L']['margins'][0]+p*lvl_paras['P']['size'][0]+lvl_paras['P']['margins'][0],
            l*lvl_paras['L']['size'][1]+lvl_paras['L']['margins'][1]+lvl_paras['P']['margins'][1]
        )
        end_point = (
            start_point[0] + lvl_paras['P']['size'][0] - lvl_paras['P']['margins'][0],
            start_point[1] + lvl_paras['P']['size'][1] - lvl_paras['P']['margins'][1]
        )

        id = graph.draw_rectangle(start_point, end_point,fill_color=lvl_paras['P']['col'][0], line_color='red', line_width=5)

        self.register_graph(id,[l,p],'P')

        self.paras['S'][l].append(0)

        self.picked['P'] = p
        self.addSynapse(l,p)


    def addSynapse(self,l,p):

        s = self.paras['S'][l][p]
        if s>=2:
            print("No adding of further synapses possible")
            return

        # print(f"adding synapse to population {p} in layer {l}")
        # s_idx = paras['P'][l] + p
        print('Synapse', s)
        self.paras['S'][l][p] += 1
        start_point = (
            lvl_paras['L']['margins'][0]+p*lvl_paras['P']['size'][0]+lvl_paras['P']['margins'][0] + s*lvl_paras['S']['size'][0]+lvl_paras['S']['margins'][0],
            l*lvl_paras['L']['size'][1]+lvl_paras['L']['margins'][1]+lvl_paras['P']['margins'][1]+lvl_paras['S']['margins'][1]
        )
        end_point = (
            start_point[0] + lvl_paras['S']['size'][0] - lvl_paras['S']['margins'][0],
            start_point[1] + lvl_paras['S']['size'][1] - lvl_paras['S']['margins'][1]
        )

        id = graph.draw_rectangle(start_point, end_point,fill_color=lvl_paras['S']['col'][0], line_color='red')

        self.register_graph(id,[l,p,s],'S')
        self.picked['S'] = s

graph = window["-GRAPH-"]
els = Elements(window, graph)

# get the graph element for ease of use later


# picked = None
while True:
    event,values = window.read()

    if event=="Exit" or event==sg.WIN_CLOSED or event is None:
        break

    if event == "-GRAPH-":
        x, y = values["-GRAPH-"]
        drag_figures = graph.get_figures_at_location((x,y))

        if len(drag_figures)>0:
            # print(f"(x,y)=({x},{y})")
            # start_point = (x, y)
            # dragging = True
            # print(drag_figures)
            els.pickObject(drag_figures)
            # lastxy = x, y
        else:
            els.picked = {}
            # end_point = (x, y)

    if event=="-ADD-L-":
        els.addLayer()
    if event=="-ADD-P-":
        els.addPopulation(els.picked['L'])
    if event=="-ADD-S-":
        els.addSynapse(els.picked['L'],els.picked['P'])
    if event=="-DELETE-ELEMENT-":
        window.extend_layout(window['-OPTIONS-'],[[sg.Button('Delete stuff',key="-DELETE-ELEMENT-")]])
        # for figure in drag_figures:
            # graph.delete_figure(figure)


window.close()
