import numpy as np
import matplotlib.colors as mcol
import matplotlib.pyplot as plt
import pandas as pd

from darkMatter import darkMatter
from general.plot_statistics import *
from general.utils import set_plot_params

import PySimpleGUI as sg

# --------------- necessary stuff for plotting --------------------
from matplotlib.widgets  import RectangleSelector
import matplotlib.figure as figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

def draw_figure_w_toolbar(canvas, fig, canvas_toolbar):
    if canvas.children:
        for child in canvas.winfo_children():
            child.destroy()
    if canvas_toolbar.children:
        for child in canvas_toolbar.winfo_children():
            child.destroy()
    figure_canvas_agg = FigureCanvasTkAgg(fig, master=canvas)
    figure_canvas_agg.draw()
    toolbar = Toolbar(figure_canvas_agg, canvas_toolbar)
    toolbar.update()
    figure_canvas_agg.get_tk_widget().pack(side='right', fill='both', expand=1)

class Toolbar(NavigationToolbar2Tk):
    def __init__(self, *args, **kwargs):
        super(Toolbar, self).__init__(*args, **kwargs)


default_params = {
    'L': {
        # layer level parameters
        'eps': 1./np.sqrt(2),
        'eta': 0.9,
        'kappa': 1.,
    },
    'P': {
        # population level parameters
        'I_ext': 1,
        'rateWnt': 1.,
        'alpha_0': 0.02,
        'tau_M': 0.01,
        'tau_n': 0.,
        'J0': 1.,
    },
    'S': {
        # psp level parameters
        'tau_I': 0.01,
        'tau_norm': 1.,
    },
    'J': {
        'J0_l': np.zeros(0),
        'in': 1.,
        'out': 1.,
    }
}

naming = {
    'L': 'layer',
    'P': 'population',
    'S': 'synapse',
    'J': 'weight'
}

options = {
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
        sg.Graph(
            canvas_size=(600, 600),
            graph_bottom_left=(0, 0),
            graph_top_right=(100, 100),
            key="-GRAPH-",
            change_submits=True,  # mouse click events
            background_color='lightblue',
            drag_submits=True)
    ]
]

sim_list = []
para_elements = {}
for lvl in ['L','P','S']:
    para_elements[lvl] = []
    for key in default_params[lvl]:
        # print(key)
        para_elements[lvl].extend([
            sg.Text(key),
            sg.In(size=(5,1),enable_events=True, key=f'-PARA-{lvl}-{key}-'),
        ])
        sim_list.append(key)

para_elements['J'] = [
    sg.Text('J_in',key=f'-PARA-TEXT-J-IN-'),
    sg.In(size=(5,1),enable_events=True, key=f'-PARA-J-in-'),
    sg.Text('J_out',key=f'-PARA-TEXT-J-OUT-'),
    sg.In(size=(5,1),enable_events=True, key=f'-PARA-J-out-'),
]

processing_column = [
    [
        sg.Text('Simulation variable 1'),
        sg.Combo(sim_list, enable_events=True, readonly=True, key='-SIM-PARA-1-KEY-'),
        sg.In(size=(4,1), enable_events=True, key="-SIM-PARA-1-LOWER-"),
        sg.Text('-'),
        sg.In(size=(4,1), enable_events=True, key="-SIM-PARA-1-UPPER-"),
        sg.Text('L'),
        sg.In(size=(2,1), enable_events=True, key="-SIM-PARA-1-L-"),
        sg.Text('P'),
        sg.In(size=(2,1), enable_events=True, key="-SIM-PARA-1-P-"),
        sg.Text('S'),
        sg.In(size=(2,1), enable_events=True, key="-SIM-PARA-1-S-"),
        sg.Button('Choose',key="-SIM-CHOOSE-1-")
    ],
    [
        sg.Text('Simulation variable 2'),
        sg.Combo(sim_list, enable_events=True, readonly=True, key='-SIM-PARA-2-KEY-'),
        sg.In(size=(4,1), enable_events=True, key="-SIM-PARA-2-LOWER-"),
        sg.Text('-'),
        sg.In(size=(4,1), enable_events=True, key="-SIM-PARA-2-UPPER-"),
        sg.Text('L'),
        sg.In(size=(2,1), enable_events=True, key="-SIM-PARA-2-L-"),
        sg.Text('P'),
        sg.In(size=(2,1), enable_events=True, key="-SIM-PARA-2-P-"),
        sg.Text('S'),
        sg.In(size=(2,1), enable_events=True, key="-SIM-PARA-2-S-"),
        sg.Button('Choose',key="-SIM-CHOOSE-2-")
    ],
    [sg.HSeparator()],
    [
        sg.Button('Add layer',key="-ADD-L-"),
        sg.Button('Add population',key="-ADD-P-"),
        sg.Button('Add synapse',key="-ADD-S-")
    ],
    [
        sg.Button('Delete stuff',key="-DELETE-ELEMENT-")
    ],
]

for lvl in ['L','P','S','J']:
    processing_column.extend([
        [sg.HSeparator()],
        [
            sg.Column([
                    [sg.Text(f'{naming[lvl]} parameters',font=('bold',18))],
                    para_elements[lvl]
                ],
                key=f'-{lvl}-OPTIONS-',visible=False)
        ],
    ])
processing_column.extend(
    [
        [sg.Canvas(key='controls_cv')],
        [sg.Canvas(key='fig_cv',
                   # it's important that you set this size
                   size=(1200, 600)
                   )]
    ],
    # background_color='#DAE0E6',
    # pad=(0, 0)
)

layout = [
    [
        sg.Column(session_list_column),
        sg.VSeparator(),
        sg.Column(processing_column,key='-OPTIONS-'),
    ],
    [
        sg.Button('Run simulation',key='-RUN-'),
        sg.Text('Steps:'),
        sg.In(size=(5,1), enable_events=True, key='-SIM-STEPS-')
    ],
]

font = ("Arial", 16)
window = sg.Window("Multilayer network analysis", layout,font=font, finalize=True)

lvl_paras = {
    'L': {
        'margins': (5,5),
        'size': (90,25),
        'col': (mcol.to_hex([0,1,0]),mcol.to_hex([0,0.8,0]))
    },
    'P': {
        'margins': (2,2),
        'size': (25,16),
        'col': (mcol.to_hex([0,0,0]),mcol.to_hex([0.2,0.2,0.2]))
    },
    'S': {
        'margins': (1,1),
        'size': (10,10),
        'col': (mcol.to_hex([0,0,1]),mcol.to_hex([0,0,0.8]))
    },
    'J': {
        'margins': (60,None),
        'size': (2,None),
        'col': (mcol.to_hex([1,0,0]),mcol.to_hex([0.8,0,0]))
    },
}


def get_element_location(lvl,ids):

    if lvl=='L':
        l,p,s = ids
        start_point = (
            lvl_paras['L']['margins'][0],
            l*lvl_paras['L']['size'][1]+lvl_paras['L']['margins'][1]
        )
        end_point = (
            start_point[0] + lvl_paras['L']['size'][0] - lvl_paras['L']['margins'][0],
            start_point[1] + lvl_paras['L']['size'][1] - lvl_paras['L']['margins'][1]
        )
    elif lvl=='P':
        l,p,s = ids
        start_point = (
            lvl_paras['L']['margins'][0]+p*lvl_paras['P']['size'][0]+lvl_paras['P']['margins'][0],
            l*lvl_paras['L']['size'][1]+lvl_paras['L']['margins'][1]+lvl_paras['P']['margins'][1]
        )
        end_point = (
            start_point[0] + lvl_paras['P']['size'][0] - lvl_paras['P']['margins'][0],
            start_point[1] + lvl_paras['P']['size'][1] - lvl_paras['P']['margins'][1]
        )
    elif lvl=='S':
        l,p,s = ids
        start_point = (
            lvl_paras['L']['margins'][0]+p*lvl_paras['P']['size'][0]+lvl_paras['P']['margins'][0] + s*lvl_paras['S']['size'][0]+lvl_paras['S']['margins'][0],
            l*lvl_paras['L']['size'][1]+lvl_paras['L']['margins'][1]+lvl_paras['P']['margins'][1]+lvl_paras['S']['margins'][1]
        )
        end_point = (
            start_point[0] + lvl_paras['S']['size'][0] - lvl_paras['S']['margins'][0],
            start_point[1] + lvl_paras['S']['size'][1] - lvl_paras['S']['margins'][1]
        )
    elif lvl=='J':
        l,ll = ids
        start_point = (
            lvl_paras['J']['margins'][0]+l*2.5+ll*7,
            lvl_paras['L']['margins'][1] + (1/2. + ll)*lvl_paras['L']['size'][1]
        )
        end_point = (
            start_point[0] + lvl_paras['J']['size'][0],
            lvl_paras['L']['margins'][1] + (1/2. + l)*lvl_paras['L']['size'][1]
        )

    return start_point, end_point

def run_code(params,sim_params,n):

    options = {}
    for lvl in ['L','P','S']:
        for key in params[lvl].keys():
            options[key] = list(flatten(params[lvl][key]))

        options[lvl] = list(flatten(n[lvl]))
    options['J0_l'] = params['J']
    options['mode'] = 0
    options['mode_calc'] = 0
    options['mode_stats'] = 0

    options['simulation'] = {}
    for idx in [1,2]:
        key = sim_params[idx]['key']
        options['simulation'][key] = sim_params[idx][key]

    options['simulation']['sim_prim'] = [sim_params[1]['L'],sim_params[1]['P'],sim_params[1]['S']]
    options['simulation']['sim_sec'] = [sim_params[2]['L'],sim_params[2]['P'],sim_params[2]['S']]

    print(options)
    order = list(options['simulation'].keys())

    res = darkMatter(steps=50,options=options,rerun=True,compile=True)

    ## general plot setup
    set_plot_params()

    fig,ax = plt.subplots(1,2)
    plt_para = {
        'ax_label': [],
        'const_label': []
    }

    for p in range(2):
        plot_fins(ax[p],res[order[0]],res[order[1]],res['gamma'][p,...],res['chi'][p,...],res['regions'][p,...],plt_para)

    draw_figure_w_toolbar(window['fig_cv'].TKCanvas, fig, window['controls_cv'].TKCanvas)


class Elements():

    register = pd.DataFrame(index=['graph_id'],columns=['lvl','L','P','S','J'])

    n = {
        'L': 0,
        'P': [],
        'S': [],
    }

    params = {
        'L': {},
        'P': {},
        'S': {},
        'J': {}
    }

    graph_ids = {'L':[],'P':[], 'S':[], 'J':[]}

    sim_params = {1: {},2: {}}

    picked = {}

    def __init__(self,window,graph):
        self.graph = graph
        self.window = window

        for lvl in ['L','P','S']:
            for key in default_params[lvl].keys():
                self.params[lvl][key] = []
        self.params['J'] = np.ones((0,0))

        self.addLayer()

    def register_graph(self,f_id,ids,lvl):
        # print(f"registering id {f_id} for level {lvl}")
        self.register.loc[f_id] = [lvl, ids[0], ids[1] if len(ids)>1 else -1, ids[2] if len(ids)>2 else -1, ids[3] if len(ids)>3 else -1]

        if lvl=='L':
            self.graph_ids['L'].append(f_id)
            self.graph_ids['P'].append([])
            self.graph_ids['S'].append([])
            self.graph_ids['J'].append([])
        elif lvl=='P':
            self.graph_ids['P'][ids[0]].append(f_id)
            self.graph_ids['S'][ids[0]].append([])
        elif lvl=='S':
            self.graph_ids['S'][ids[0]][ids[1]].append(f_id)
        elif lvl=='J':
            self.graph_ids['J'][ids[3]].append(f_id)
        # print('GRAPH_ID:',self.graph_ids)

    def get_param_value(self,lvl,key,ids):
        print(ids)
        if lvl=='L':
            val = self.params[lvl][key][ids['L']]
        elif lvl=='P':
            val = self.params[lvl][key][ids['L']][ids['P']]
        elif lvl=='S':
            val = self.params[lvl][key][ids['L']][ids['P']][ids['S']]
        elif lvl=='J':
            if key=='in':
                val = self.params[lvl][ids['L']][ids['J']]
            elif key=='out':
                val = self.params[lvl][ids['J']][ids['L']]
        return val

    def pickObject(self,f_id):

        for lvl in ['L','P','S','J']:
            window[f'-{lvl}-OPTIONS-'].update(visible=False)
        if bool(self.picked):
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

            if 'J' in self.picked.keys():
                graph_id = np.where((self.register['lvl']=='J') &
                    (self.register['L']==self.picked['L']) &
                    (self.register['J']==self.picked['J']))[0]

                if len(graph_id):
                    self.graph.Widget.itemconfig(graph_id[0], fill=lvl_paras['J']['col'][0])

        break_it = False
        for lvl in ['J','S','P','L']:
            for id in f_id:
                if lvl==self.register.loc[id]['lvl']:
                    self.picked[lvl] = self.register.loc[id][lvl]

                    self.window[f'-{lvl}-OPTIONS-'].update(visible=True)

                    if lvl=='J':
                        self.picked['L'] = self.register.loc[id]['L']
                        self.window[f'-PARA-TEXT-J-IN-'].update(f'J_{self.picked["L"]},{self.picked["J"]}')
                        self.window[f'-PARA-TEXT-J-OUT-'].update(f'J_{self.picked["J"]},{self.picked["L"]}')

                        for key in ['in','out']:
                            val = self.get_param_value(lvl,key,self.register.loc[id])
                            print('val:',val)
                            self.window[f'-PARA-{lvl}-{key}-'].update(val)
                    else:
                        for key in self.params[lvl].keys():
                            val = self.get_param_value(lvl,key,self.register.loc[id])
                            self.window[f'-PARA-{lvl}-{key}-'].update(val)

                    self.graph.Widget.itemconfig(id, fill=lvl_paras[lvl]['col'][1])
                    print('picked lvl',lvl)
                    if lvl=='J':
                        break_it = True
                        break
            if break_it:
                break

    def choose_element(self,graph_id):

        choice = {}
        break_it = False
        for lvl in ['S','P','L']:
            for id in graph_id:
                if self.register.loc[id]['lvl'] == lvl:
                    # print('id:',id)
                    # print(self.register.loc[id])
                    choice['L'] = self.register.loc[id]['L']
                    choice['P'] = self.register.loc[id]['P']
                    choice['S'] = self.register.loc[id]['S']
                    break_it = True
                    break
            if break_it:
                break
        return choice, id

    def delete_element(self,graph_id):
        print('now deleting ',graph_id)

        lvl = self.register.loc[graph_id]['lvl']
        l = self.register.loc[graph_id]['L']
        p = self.register.loc[graph_id]['P']
        s = self.register.loc[graph_id]['S']
        print('ids:',l,p,s)

        # going through ...
        #   ... parameters to delete entries
        #   ... number array to delete
        #   ... children elements to delete
        if lvl=='L':
            for key in self.params['L'].keys():
                self.params['L'][key].pop(l)
            for key in self.params['P'].keys():
                self.params['P'][key].pop(l)
            for key in self.params['S'].keys():
                self.params['S'][key].pop(l)
            self.n['L'] -= 1
            self.n['P'].pop(l)
            self.n['S'].pop(l)

            for id in list(flatten(self.graph_ids['S'][l])):
                print('delete:',id)
                graph.delete_figure(id)
            for id in self.graph_ids['P'][l]:
                print('delete:',id)
                graph.delete_figure(id)
        elif lvl=='P':
            for key in self.params['P'].keys():
                self.params['P'][key][l].pop(p)
            for key in self.params['S'].keys():
                self.params['S'][key][l].pop(p)

            self.n['P'][l] -= 1
            self.n['S'][l].pop(p)

            for id in self.graph_ids['S'][l][p]:
                print('delete:',id)
                graph.delete_figure(id)
        elif lvl=='S':
            for key in self.params['S'].keys():
                self.params['S'][key][l][p].pop(s)

            self.n['S'][l][p] -= 1

        # finally, deleting element
        graph.delete_figure(graph_id)

        # relocate all other elements
        for id,row in self.register.iterrows():
            # print('------------------------')
            # print(id,row)
            if not np.isnan(row['L']):
                start_point, end_point = get_element_location(row['lvl'],[row['L'],row['P'],row['S']])
                # print(start_point, end_point)
                # graph.Relocate(id,start_point,end_point)
                # row = self.register.loc[id]



    def change_params(self,lvl,key,val):

        print('change params:',lvl,key)
        try:
            val = (int if type(default_params[lvl][key])==int else float)(val)
            if lvl=='L':
                self.params[lvl][key][self.picked['L']] = val
            elif lvl=='P':
                self.params[lvl][key][self.picked['L']][self.picked['P']] = val
            elif lvl=='S':
                self.params[lvl][key][self.picked['L']][self.picked['P']][self.picked['S']] = val
            elif lvl=='J':
                if key=='in':
                    self.params[lvl][self.picked['L']][self.picked['J']] = val
                if key=='out':
                    self.params[lvl][self.picked['J']][self.picked['L']] = val

        except:
            val = self.get_param_value(lvl,key,self.picked)
            self.window[f'-PARA-{lvl}-{key}-'].update(val)

    def change_sim_params(self,idx,vals):

        key = vals[f'-SIM-PARA-{idx}-KEY-']

        lower = vals[f'-SIM-PARA-{idx}-LOWER-']
        upper = vals[f'-SIM-PARA-{idx}-UPPER-']
        try:
            lower = 0 if lower=='' else float(lower)
            upper = 0 if upper=='' else float(upper)
        except:
            lower = 0
            upper = 0
        self.sim_params[idx] = {key: [lower,upper], 'key':key}
        # self.window[f'-SIM-PARA-{idx}-LOWER-'].update(lower)
        # self.window[f'-SIM-PARA-{idx}-UPPER-'].update(upper)

        for lvl in ['L','P','S']:
            p = vals[f'-SIM-PARA-{idx}-{lvl}-']
            p = 0 if p=='' else int(p)
            self.sim_params[idx][lvl] = p
            self.window[f'-SIM-PARA-{idx}-{lvl}-'].update(p)


    def addLayer(self):
        # print('adding a new layer')
        l = self.n['L']
        if l>=4:
            print("No adding of further layers possible")
            return

        self.n['L'] += 1
        start_point, end_point = get_element_location('L',[l,None,None])

        id = graph.draw_rectangle(start_point, end_point,fill_color=lvl_paras['L']['col'][0], line_color='red')
        self.register_graph(id,[l],'L')

        self.n['P'].append(0)
        self.n['S'].append([])

        for key in default_params['L'].keys():
            self.params['L'][key].append(default_params['L'][key])
        for key in default_params['P'].keys():
            self.params['P'][key].append([])
        for key in default_params['S'].keys():
            self.params['S'][key].append([])

        self.pickObject([id])

        self.addWeights(l)
        self.addPopulation(l)
        self.addPopulation(l)
        print(self.params)

    def addWeights(self,l):

        print(self.params)

        if self.n['L'] > 1:
            for ll in range(self.n['L']):
                if ll!=l:
                    start_point, end_point = get_element_location('J',[l,ll])

                    id = graph.draw_rectangle(start_point,end_point,fill_color=lvl_paras['J']['col'][0])
                    self.register_graph(id,[l,None,None,ll],'J')

                    # win_key = f'-PARA-J-{l}-{ll}-'
                    # if not (win_key in self.window.AllKeysDict):
                    #     self.window.extend_layout(self.window['-J-OPTIONS-'],[[
                    #         sg.Text(f'J_{ll}'),
                    #         sg.In(size=(5,1),enable_events=True, key=f'-PARA-J-{l}-{ll}-'),
                    #     ]])

        J_l = np.ones((self.n['L'],self.n['L']))*default_params['J']['in']
        np.fill_diagonal(J_l,0)
        for l,J_row in enumerate(self.params['J']):
            for ll,el in enumerate(J_row):
                J_l[l,ll] = el
        self.params['J'] = J_l


    def addPopulation(self,l):

        p = self.n['P'][l]
        if p>=2:
            print("No adding of further populations possible")
            return
        # print(f"adding population to layer {l}")
        self.n['P'][l] += 1
        start_point, end_point = get_element_location('P',[l,p,None])

        id = graph.draw_rectangle(start_point, end_point,fill_color=lvl_paras['P']['col'][0], line_color='red', line_width=5)

        self.register_graph(id,[l,p],'P')

        self.n['S'][l].append(0)
        for key in default_params['P'].keys():
            self.params['P'][key][l].append(default_params['P'][key])

        for key in default_params['S'].keys():
            self.params['S'][key][l].append([])

        self.pickObject([id])
        self.addSynapse(l,p)


    def addSynapse(self,l,p):

        s = self.n['S'][l][p]
        if s>=2:
            print("No adding of further synapses possible")
            return

        # print(f"adding synapse to population {p} in layer {l}")
        # s_idx = paras['P'][l] + p
        self.n['S'][l][p] += 1
        start_point, end_point = get_element_location('S',[l,p,s])

        id = graph.draw_rectangle(start_point, end_point,fill_color=lvl_paras['S']['col'][0], line_color='red')

        for key in default_params['S'].keys():
            self.params['S'][key][l][p].append(default_params['S'][key])

        self.register_graph(id,[l,p,s],'S')
        self.pickObject([id])

graph = window["-GRAPH-"]
els = Elements(window, graph)

# get the graph element for ease of use later

def flatten(element):
    if isinstance(element, list):
        for item in element:
            yield from flatten(item)
    else:
        yield(element)

interact = {
    'active': False,
    'action': '',
    'id': None,
}

# choosing = False
# deleting = False
while True:
    event,values = window.read()

    if event=="Exit" or event==sg.WIN_CLOSED or event is None:
        break

    if event == "-GRAPH-":
        x, y = values["-GRAPH-"]
        drag_figures = graph.get_figures_at_location((x,y))

        if len(drag_figures)>0:
            if interact['active']:
                if interact['action'] == 'choosing':
                    choice,_ = els.choose_element(drag_figures)
                    for lvl in ['L','P','S']:
                        window[f'-SIM-PARA-{interact["id"]}-{lvl}-'].update(choice[lvl])
                        els.sim_params[interact['id']][lvl] = choice[lvl]
                if interact['action'] == 'deleting':
                    choice,id = els.choose_element(drag_figures)
                    els.delete_element(id)

                window.set_cursor('arrow')
                interact = {'active': False, 'action': '', id: None}
            else:
                els.pickObject(drag_figures)
        else:
            els.picked = {}

    if event=="-ADD-L-":
        els.addLayer()
    if event=="-ADD-P-":
        els.addPopulation(els.picked['L'])
    if event=="-ADD-S-":
        els.addSynapse(els.picked['L'],els.picked['P'])
    if event=="-DELETE-ELEMENT-":
        interact = {'active': True, 'action': 'deleting', 'id': None}
        window.set_cursor('cross')

    if event.startswith('-PARA-'):
        lvl = event.split('-')[-3]
        key = event.split('-')[-2]

        els.change_params(lvl,key,values[event])

    if event=='-SIM-STEPS-':
        print('changed steps')

    if event.startswith('-SIM-PARA-'):
        idx = int(event.split('-')[-3])

        els.change_sim_params(idx,values)

    if event.startswith('-SIM-CHOOSE-'):
        # change behavior on click
        interact = {'active': True, 'action': 'choosing', 'id': int(event.split('-')[-2])}
        # change look of mouse pointer
        window.set_cursor('hand2')

    if event=="-RUN-":
        run_code(els.params,els.sim_params,els.n)

window.close()
