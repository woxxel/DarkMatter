import numpy as np
from numpy.ma import masked_array
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from scipy.io import netcdf
import os
import matplotlib.gridspec as gridspec

def haifische(steps=100,rate=[0.,20.],alpha_0=[0,0.16],tau_G=[0.005],n=[0.],eps=0):
  
  network = {}
  
  network['tau_A'] = 0.005				#the excitatory time constant for AMPA currents
  network['tau_N'] = 0.200				#the excitatory time constant for NMDA currents
  network['tau_M'] = 0.010				#the membrane time constant of all neurons
  
  # alpha vs nu (sharkfins)
  rate = [[0,20]]#,[0,15],[0,10],[0,7],[0,6.5]]
  alpha_0 = [[0,0.16],[0,0.12],[0,0.08],[0,0.05],[0,0.035]]
  tau_G = [[0.005],[0.01],[0.02],[0.04],[0.06]]
  
  ## tau vs nu
  #rate = [[0,20],[0,20],[0,20]]
  #alpha_0 = [[0.01],[0.02],[0.04]]
  #tau_G = [[0,0.1],[0,0.1],[0,0.1]]
  
  ### alpha vs tau
  #rate = [[1]]#[[1],[2],[5]]
  #alpha_0 = [[0,0.1],[0,0.1],[0,0.1]]
  #tau_G = [[0,0.1],[0,0.1],[0,0.1]]
  
  
  n = [[0],[0],[0],[0],[0]]
  
  
  
  #fig, ax = plt.subplots(len(rate)/2+1,2,figsize=(6.75,1+2.3*(len(rate)/2+1)))
  fig = plt.figure(figsize=(6,5))
  
  for i in range(len(rate)):
    #ax[0] = plt.subplot(gs[i/2,i%2])
  
    inpara = {}
    inpara['rateWnt'] = rate[i]
    inpara['alpha_0'] = alpha_0[i]
    inpara['tau_G'] = tau_G[i]				#the inhibitory time constant for GABA currents
    inpara['n'] = n[i]
    
    network['kappa'] = 1
    network['eta'] = 0.9
    network['eps'] = eps
    
    for key in inpara:
      assert type(inpara[key]) == list, 'Please specify all iteration parameters as lists!'
      #assert len(inpara[key]) in [1,2], 'Specify 1 or 2 values for iteration parameters!'
      if len(inpara[key]) == 2:
	network[key] = np.linspace(inpara[key][0] + inpara[key][-1]/float(steps),inpara[key][-1],steps).astype('float')
      else:
	network[key] = np.array(inpara[key]).astype('float')
    para_list = ['n','alpha_0','tau_G','rateWnt']
    idx = 0
    ax_list = []
    const_labels = []
    sv_str = 'steps=%d_' % steps
    sv_str += 'const'
    for key in para_list:
      if len(network[key]) == 1:
	ax_list.append(idx)
	sv_str += '_%s=%g' % (key,network[key][0])
	if key == 'n':
	  const_labels.append(r'$\displaystyle n = %g$' % network[key])
	if key == 'alpha_0':
	  const_labels.append(r'$\displaystyle \alpha_0 = %g$' % network[key])
	if key == 'tau_G':
	  const_labels.append(r'$\displaystyle \tau_I = %g\,$s' % network[key])
	if key == 'rateWnt':
	  const_labels.append(r'$\displaystyle \bar{\nu} = %g\,$Hz' % network[key])
      idx += 1
    
    idx = 0
    #ax_vals = []
    ax_labels = []
    sv_str += '_iter'
    for key in para_list:
      if len(network[key]) == steps:
	sv_str += '_%s=%g' % (key,network[key][-1])
	ax_list.append(idx)
	#ax_vals.append(key)
	if key == 'n':
	  ax_labels.append(r'$\displaystyle n$')
	if key == 'alpha_0':
	  ax_labels.append(r'$\displaystyle \alpha_0$')
	if key == 'tau_G':
	  ax_labels.append(r'$\displaystyle \tau_G$ in s')
	if key == 'rateWnt':
	  ax_labels.append(r'$\displaystyle \bar{\nu}$ in Hz')
	
      idx += 1
    ax_list.append(4)
    
    #print ax_list
    #print network
    fileNet = 'parameters.nc'
    fileResults = './data/results_%s.nc' % sv_str
    
    if not os.path.exists(fileResults):
      print 'saving file as %s ...' % fileResults

      ncid = netcdf.netcdf_file(fileNet,'w')
      
      ncid.createDimension('one',1)
      
      for key in network.keys():
	if type(network[key]) == np.ndarray:
	  var_dim = key + '_dim'
	  ncid.createDimension(var_dim,len(network[key]))
	  var_type = network[key][0].dtype.char
	else:
	  var_dim = 'one'
	  var_type = type(network[key])
	  var_type = 'd'

	Var = ncid.createVariable(key,var_type,(var_dim,))
	
	Var[:] = network[key]
      
      ncid.close()
      
        
      run_str = './Haifische ' + fileNet + ' ' + fileResults
      os.system(run_str)
      
      os.remove(fileNet)
    
    print 'reading from file %s ...' % fileResults

    ncid = netcdf.netcdf_file(fileResults,'r',mmap=None)
    
    gamma = np.zeros((steps,steps,2))
    chi = np.zeros((steps,steps,2))
    
    gamma[:] = ncid.variables['gamma'][:].transpose(ax_list)
    chi[:] = ncid.variables['chi'][:].transpose(ax_list)
    gamma = gamma[:,:,0]
    chi = chi[:,:,0]
    
    plot_matrix = np.zeros((steps,steps))
    
    mask_inconsistent = (gamma == -1)
    mask_no_peak = (gamma == -2)
    mask_dark_matter = (gamma < 1)
    plot_gamma = masked_array(1 - gamma**2,mask_inconsistent + mask_no_peak)#masked_array(gamma,gamma>0)
    plot_chi = masked_array(chi,mask_inconsistent + mask_no_peak + mask_dark_matter)
    plot_regions = masked_array(gamma,np.invert(mask_inconsistent + mask_no_peak))
    
    ncid.close()
    
    levs = range(160)
    bnw = mcolors.LinearSegmentedColormap.from_list(name='black_n_white',colors=[(0,(1,1,1)),(1,(0,0,0))],N=len(levs)-1)
    heat = mcolors.LinearSegmentedColormap.from_list(name='heat',colors=[(0,'y'),(0.5,'r'),(1,(0.1,0.1,0.5))],N=len(levs)-1)
    
    plt.rc('text', usetex=True)
    plt.rc('font', family='Times')
    
    ax = plt.axes([0.12,0.1,0.7,0.8])
    
    #ax[i/2,i%2].tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='both', which='major', labelsize=10)
    
    #plt.figure(figsize=(cm2inch(6),cm2inch(6)))
    #ax = plt.axes([0.2,0.16,0.7,0.75])
    
    #plt.figure(figsize=(cm2inch(10),cm2inch(8)))
    #ax = plt.axes([0.15,0.15,0.6,0.75])
    
    #ax[i/2,i%2].pcolormesh(network[para_list[ax_list[3]]],network[para_list[ax_list[2]]],plot_no_peak,cmap=bnw,vmin=-1,vmax=2)
    #ax[i/2,i%2].pcolormesh(network[para_list[ax_list[3]]],network[para_list[ax_list[2]]],plot_inconsistent,cmap=bnw,vmin=-1,vmax=2)
    
    bnw.set_bad('k',0.)
    #pgamma = ax[i/2,i%2].pcolormesh(network[para_list[ax_list[3]]],network[para_list[ax_list[2]]],plot_gamma,cmap=bnw,vmin=-1,vmax=1)
    pgamma = ax.pcolormesh(network[para_list[ax_list[3]]],network[para_list[ax_list[2]]],plot_gamma,cmap=bnw,vmin=-1,vmax=1)
    heat.set_bad('k',0.)
    #pchi = ax[i/2,i%2].pcolormesh(network[para_list[ax_list[3]]],network[para_list[ax_list[2]]],plot_chi,cmap=heat,vmin=0,vmax=3)
    pchi = ax.pcolormesh(network[para_list[ax_list[3]]],network[para_list[ax_list[2]]],plot_chi,cmap=heat,vmin=0,vmax=3)
    #bnw.set_bad((0.9,0.9,0.9),1.)
    #pregions = ax[i/2,i%2].pcolormesh(network[para_list[ax_list[3]]],network[para_list[ax_list[2]]],plot_regions,cmap=bnw,vmin=-2,vmax=3)
    pregions = ax.pcolormesh(network[para_list[ax_list[3]]],network[para_list[ax_list[2]]],plot_regions,cmap=bnw,vmin=-2,vmax=3)
    
    idx = 0
    
    #ax[i/2,i%2].set_xlim([0,network[para_list[ax_list[3]]][-1]])
    #ax[i/2,i%2].set_ylim([0,network[para_list[ax_list[2]]][-1]])
    ax.set_xlim([0,network[para_list[ax_list[3]]][-1]])
    ax.set_ylim([0,network[para_list[ax_list[2]]][-1]])
    for key in para_list:
      if key in [para_list[ax_list[2]],para_list[ax_list[3]]]:
	ticks = 11
	if key == 'n':
	  array = np.linspace(0,1,6)
	elif key == 'alpha_0':
	  #print network[key][-1]
	  if network[key][-1] > 0.1:
	    ticks = 6
	  array = np.linspace(0,0.2,ticks)
	  #print array
	elif key == 'tau_G':
	  if network[key][-1] > 0.05:
	    ticks = 6
	  array = np.linspace(0,0.1,ticks)
	elif key == 'rateWnt':
	  if network[key][-1] > 10:
	    ticks = 6
	  array = np.linspace(0,20,ticks)
	
	if idx:
	  try:
	    x_array = array[:np.where(array>network[key][-1])[0][0]]
	  except:
	    x_array = array[:]
	  #ax[i/2,i%2].set_xticks(x_array)
	  #ax[i/2,i%2].set_xticklabels(x_array,fontsize=10)
	  ax.set_xticks(x_array)
	  ax.set_xticklabels(x_array,fontsize=12)
	else:
	  try:
	    y_array = array[:np.where(array>network[key][-1])[0][0]]
	  except:
	    y_array = array[:]
	  #ax[i/2,i%2].set_yticks(y_array)
	  #ax[i/2,i%2].set_yticklabels(y_array,fontsize=10)
	  ax.set_yticks(y_array)
	  ax.set_yticklabels(y_array,fontsize=12)
	idx += 1
    

    #ax[i/2,i%2].set_xlabel(ax_labels[1],fontsize=12)
    #ax[i/2,i%2].set_ylabel(ax_labels[0],fontsize=12)
  
  
    #ax[i/2,i%2].set_title('%s' % (const_labels[1]),fontsize=12)
    #ax.set_title('%s' % (const_labels[1]),fontsize=12)
 
  axcb1 = plt.axes([0.85,0.35,0.05,0.55])
  axcb2 = plt.axes([0.85,0.1,0.05,0.2])  
  
  axcb1.tick_params(axis='both', which='major', labelsize=12)
  axcb2.tick_params(axis='both', which='major', labelsize=12)
  
  plt.colorbar(pchi, cax = axcb1,boundaries=np.linspace(0,3,100),ticks=np.linspace(0,3,7))
  plt.colorbar(pgamma, cax = axcb2,boundaries=np.linspace(0,1,10),ticks=np.linspace(0,1,3))
  
  #axcb2.set_yticks(np.linspace(1,0,3))
  #axcb2.set_yticklabels(np.linspace(1,0,3))
  axcb1.set_ylabel(r'$\displaystyle \chi$',fontsize=14)
  axcb2.set_ylabel(r'$\displaystyle 1 - \gamma^2$',fontsize=14)
  
    #print const_labels[1]
    #plt.show()
  #plt.subplots_adjust(left=0.125, bottom=0.075, right=0.8, top=0.95, wspace=0.3, hspace=0.3)

  #left  = 0.  # the left side of the subplots of the figure
  #right = 0.9    # the right side of the subplots of the figure
  #bottom = 0.05   # the bottom of the subplots of the figure
  #top = 0.95      # the top of the subplots of the figure
  #wspace = 0.2   # the amount of width reserved for blank space between subplots
  #hspace = 0.5   # the amount of height reserved for white space between subplots_adjust
  
  big_ax = fig.add_subplot(111)
  big_ax.set_axis_bgcolor('none')
  big_ax.tick_params(labelcolor='none',top='off',bottom='off',left='off',right='off')
  big_ax.spines['top'].set_visible(False)
  big_ax.spines['right'].set_visible(False)
  big_ax.spines['bottom'].set_visible(False)
  big_ax.spines['left'].set_visible(False)
  
  ax.set_xlabel(ax_labels[1],fontsize=16)
  ax.set_ylabel(ax_labels[0],fontsize=16)
  
  #plt.suptitle('%s, %s' % (const_labels[0],const_labels[1]),fontsize=12)
  #sv_name = '../inhib only/pics/hai_%s_%s_%s=%g_%s=%g.png' % (para_list[ax_list[2]],para_list[ax_list[3]],para_list[ax_list[0]],network[para_list[ax_list[0]]],para_list[ax_list[1]],network[para_list[ax_list[1]]])
  
  
  sv_name = '../pics/hai_steps=%d_%s_%s.png' % (steps,para_list[ax_list[2]],para_list[ax_list[3]])
  print 'Figure saved as "%s"' % sv_name
  plt.savefig(sv_name,dpi=300)
  plt.show(block=False)
  
  #plt.figure(figsize=(cm2inch(2),cm2inch(20)))
  #axcb1 = plt.axes([0.1,0.35,0.35,0.55])
  #axcb2 = plt.axes([0.1,0.1,0.35,0.2])
  
  #ax[i/2,i%2].tick_params(axis='both', which='major', labelsize=10)
  #axcb1.tick_params(axis='both', which='major', labelsize=10)
  #axcb2.tick_params(axis='both', which='major', labelsize=10)
  
  #plt.colorbar(pchi, cax = axcb1,boundaries=np.linspace(0,3,100),ticks=np.linspace(0,3,6))
  #plt.colorbar(pgamma, cax = axcb2,boundaries=np.linspace(0,1,100),ticks=np.linspace(0,1,3))
  #axcb1.set_ylabel(r'$\displaystyle \chi$',fontsize=12)
  #axcb2.set_ylabel(r'$\displaystyle \gamma^2 - 1$',fontsize=12)
  #plt.savefig('../inhib only/pics/cbar.png',format='png',dpi=300)
  #plt.show()
  
  #os.remove(fileResults)
  
def do_all():
  haifische(steps=500,rate=[0,20],alpha_0=[0,0.16],tau_G=[0.005],n=[0])
  haifische(steps=500,rate=[0,15],alpha_0=[0,0.12],tau_G=[0.01],n=[0])
  haifische(steps=500,rate=[0,10],alpha_0=[0,0.08],tau_G=[0.02],n=[0])
  haifische(steps=500,rate=[0,7],alpha_0=[0,0.05],tau_G=[0.04],n=[0])
  haifische(steps=500,rate=[0,6.5],alpha_0=[0,0.035],tau_G=[0.06],n=[0])
  
def cm2inch(value):
  return value/2.54