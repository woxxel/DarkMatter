





## I.load_data('empirical',filePath='../../data/BuscheLab/2P_data.xlsx')
data_observed = I.data[I.data_mask]
#data_observed = data_observed[data_observed>0.01]

para_steps = 51
paras = {
    'gamma': 1.5,
    'delta': 5.,
    'nu_max': 40.
}

sim = {
    'first': {
        'key':'gamma',
        'val': np.linspace(0.5,3.,para_steps)
    },
    'second': {
        'key': 'nu_max',
        'val': np.linspace(10,100.,para_steps)
    }
}

paras[sim['first']['key']] = sim['first']['val']
paras[sim['second']['key']] = sim['second']['val']

d = np.tile(data_observed,(para_steps,1,1))

for key in paras:
    paras[key] = np.tile(paras[key],(1,1,1)).T

paras[sim['first']['key']] = np.transpose(paras[sim['first']['key']],(1,0,2))

logp = logp_nu(d,paras['gamma'],paras['delta'],paras['nu_max'])
print(logp.shape)
print(d.shape)
#logp[...,data_observed<0.01] = logp[...,data_observed<0.01]#*0.1
logp_sum = logp.sum(axis=2)

X,Y = np.meshgrid(sim['first']['val'],sim['second']['val'])






import ipywidgets as widgets
from mpl_toolkits.mplot3d import Axes3D

plt.ion()
fig = plt.figure(figsize=(10,4))

ax_3d = fig.add_subplot(1, 2, 1, projection='3d')
ax_cumsum = fig.add_subplot(1, 2, 2)

ax_3d.plot_surface(X,Y,logp_sum,alpha=0.5)

global pos, dist
pos = ax_3d.scatter(sim['first']['val'][15],sim['second']['val'][15],logp_sum[15,15],color='r')
plt.setp(ax_3d,xlabel=f'$\\{sim["first"]["key"]}$',ylabel=f'$\\{sim["second"]["key"]}$',zlim=[-10000,2000])

x_arr = np.linspace(0,20,10001)
global par
par = {
    'gamma': 1.5,
    'delta': 5.,
    'nu_max': 40
}
ax_dist = ax_cumsum.twinx()
dist, = ax_dist.plot(x_arr,p_nu(x_arr,**par),'k-')
ax_dist.set_yscale('log')

#plt.scatter(np.arange(len(data_observed)),np.sort(data_observed))
idx_sort = np.argsort(data_observed)
#logp_cum, = ax_cumsum.plot(data_observed[idx_sort], np.cumsum(logp[15,15,idx_sort]),label='cumsum')
logp_cum, = ax_cumsum.plot(data_observed[idx_sort], logp[15,15,idx_sort],color='r',marker='.',linestyle='none',label='cumsum')


def update_plot(idx_first,idx_second):
    global pos, par
    plt.pause(0.1)
    gamma = sim['first']['val'][idx_first]
    nu_max = sim['second']['val'][idx_second]
    print(f"{sim['first']['key']}={gamma}, {sim['second']['key']}={nu_max}, logp={logp_sum[idx_second,idx_first]}")

    #print(sim['first']['val'][idx_first])
    #cumsum = np.cumsum(logp[idx_second,idx_first])
    cumsum = logp[idx_second,idx_first,idx_sort]
    print(f'cumulative sum: {cumsum}')
    logp_cum.set_ydata(cumsum)

    pos.remove()
    pos = ax_3d.scatter(sim['first']['val'][idx_first],sim['second']['val'][idx_second],logp_sum[idx_second,idx_first],color='r',marker='D',s=20)

    par[sim['first']['key']] = sim['first']['val'][idx_first]
    par[sim['second']['key']] = sim['second']['val'][idx_second]

    dist.set_ydata(p_nu(x_arr,**par))

    #pos._offsets3d = (gamma,
    #                  nu_max,
    #                  logp_sum[idx_first,idx_second])
    #pos.set_xdata(sim['first']['val'][idx_first])
    #pos.set_ydata(sim['second']['val'][idx_second])
    #pos.set_zdata(logp_sum[idx_first,idx_second])
    #plt.setp(ax_cumsum,ylim=[min(np.nanmin(cumsum),
    #                             ax_cumsum.get_ylim()[0]),
    #                         max(np.nanmax(cumsum),
    #                             ax_cumsum.get_ylim()[1])])
    plt.setp(ax_dist,ylim=[10**-4,10**2])
    plt.setp(ax_cumsum,ylim=[np.log(10**-4),np.log(10**2)])

    #print(np.nanmin(cumsum))
    #print(np.nanmax(cumsum))
    #print(ax_cumsum.get_ylim())
    #plt.draw()

plt.legend()
#plt.show()

first_widget = widgets.IntSlider(15,min=0,max=para_steps,step=1,orientation='horizontal',description=f"$\displaystyle \\{sim['first']['key']}$")
second_widget = widgets.IntSlider(15,min=0,max=para_steps,step=1,orientation='horizontal',description=f"$\displaystyle \\{sim['second']['key']}$")
#nu_max_widget = widgets.IntSlider(0,min=0,max=para_steps,step=0.1,orientation='horizontal',description=r'$\displaystyle \nu_{max}$')

widgets.interactive(update_plot,idx_first=first_widget,idx_second=second_widget)
#np.isnan(logp[13,9]).sum()
#plt.tight_layout()
