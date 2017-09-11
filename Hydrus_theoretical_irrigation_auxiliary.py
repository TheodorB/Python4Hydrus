#!/home/theodor/anaconda3/bin/python
"""
Created on Wed Apr 26 13:11:59 2017
@author: theodor
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sb
import datetime
import get_hydrus as gh # MY HYDRUS IMPORT FUNCTIONS
import hydrus_in_out_functions as hh # MY HYDRUS CHANGE DATA FUNCTIONS

# Define some functions

def run_hydrus_read_results(run_path, running_hydrus=True, print_running=True):
    '''
    run hydrus with current settings
    print_running - if u want to print message every run
    '''
    if running_hydrus:
        hyd = gh.hydrus_handler(run_path)
        if print_running:
            print ('running hydrus')
        hyd.run_hydrus()
    hydrus_dict = dict()
    hydrus_dict['hydrus_dat'] = hyd.get_hydrus_dat() 
    hydrus_dict['node_out']  = hyd.get_node_out()
    hydrus_dict['level_out']   = hyd.get_t_level_out()
    hydrus_dict['obs_node_out'] = hyd.get_obs_node_out()
    hydrus_dict['profile_dat'] = dict(zip(['profdatdf','textpd'], hyd.get_profile_dat()))
    return hydrus_dict, hyd

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]

def rdf(z,Zm, p, zstar):
    '''
    root distribution exponential function
    z - depth
    Zm - maxumum rooting depth
    p, zstar = empirical parameters
    
    '''
    a = 1 - (z/Zm)
    b = -(p/Zm)*np.abs(zstar-z)
    c = a*np.exp(b)
    c[c<0.0] = 0.0
    return c

def plot_node_data(run_path,labels, _mypars, hydrus_dict, par = 'Moisture'):
    '''
    PLOT NODE_OUT DATA
    '''
    mydata, hdict, obsdict = extract_data(hydrus_dict)
    # par   = _par.value
#    par = 'Moisture'
    xlabel = labels[_mypars.index(par)]
    fig1, ax1 = plt.subplots()
    fig1.set_size_inches(7,7)
    cm = sb.color_palette('viridis',len(mydata.keys()))
    for n,t in enumerate(sorted(mydata.keys())):
#        label_ = str(t)+' '+hdict['TimeUnit']
        label_ = str(t)+' '+hdict['TimeUnit']
        ax1.plot(mydata[t][par], mydata[t]['Depth'],alpha = 0.5,linewidth=3,color=cm[n],label=label_)
    ax1.set_xlabel(xlabel)
    ax1.legend(loc=7,bbox_to_anchor=(1.45,0.5), frameon=True,  fancybox=True, framealpha=0.3)
    ax1.set_ylabel('$Depth \, ('+ hdict['SpaceUnit'] +')$')
    xmin, xmax = ax1.get_xlim()
    ymin, ymax = ax1.get_ylim()
    ax1.autoscale(enable=True, axis='y', tight=True)
    if par in ['Moisture', 'Conc(1..NS)','Head','Temp']:
        # plot temporal
#        what = 'time'

        if par =='Moisture':
            npar = 'theta'
        elif par == 'Head':
            npar = 'h'
        elif par == 'Conc(1..NS)':
            npar = 'Conc'
        elif par == 'Flux':
            npar = 'flux'
        else: npar = par
    xlabel = labels[_mypars.index(par)]
    fig2, ax2 = plt.subplots()
    fig2.set_size_inches(10,4)
    cm = sb.color_palette('viridis',len(([float(i) for i in obsdict.keys()])[::-1]))
    
    for n,d in enumerate(sorted([float(i) for i in obsdict.keys()])[::-1]): # iterate sorted depths
        ax2.plot(obsdict[d]['time'], obsdict[d][npar],alpha = 0.5,linewidth=3,color=cm[n],label=str(d)+' '+hdict['SpaceUnit'])
#    ax2.set_xlabel('Time ('+ hdict['TimeUnit'] +')')
    ax2.set_xlabel('Time (h)')
    ax2.legend(loc=7,bbox_to_anchor=(1.25,0.5),frameon=True, fancybox=True, framealpha=0.3)
    ax2.set_ylabel(xlabel)
    xmin, xmax = ax2.get_xlim()
    ymin, ymax = ax2.get_ylim()
    ax2.autoscale(enable=True, axis='x', tight=True) 
    return fig1, ax1, fig2, ax2


def plot_ET(hydrus_dict):
    '''
    PLOT POTENTIAL AND ACTUAL ET
    '''
    mydata, hdict, obsdict = extract_data(hydrus_dict)
    actual    = hydrus_dict['level_out']['vRoot'].values
    potential = hydrus_dict['level_out']['rRoot'].values
    time = hydrus_dict['level_out']['Time']
    fig, ax = plt.subplots(figsize=(10,4))
    # multiply by 10 for cm to mm
    ax.plot(time, 10*potential,alpha=0.5, label='T$_p$') # linestyle='--',marker='o',markersize=5,
    ax.plot(time, 10*actual,alpha=0.5, label='T$_a$')

#    ax.legend(loc=0,bbox_to_anchor=(1.005,0.6),frameon=True, fancybox=True, framealpha=0.5)

    ax.set_xlabel('Time (h)')
    ax.set_ylabel('[mm h$^{-1}$]')
    ax.autoscale(enable=True, axis='x', tight=True)
    ax.text(0.1,0.9,'T$_a$ / T$_p$ = %.2f' %(np.nansum(actual)/np.nansum(potential)),
                             ha='left',
                             va='bottom',
                             transform = ax.transAxes)
    
    ax.legend(loc=7,bbox_to_anchor=(1.15,0.5),frameon=True, fancybox=True, framealpha=0.5)

    return fig, ax
    # print np.sum(actual)
    # fig.savefig(os.path.join(run_path,'transpiration'+ '.png'),dpi=250,bbox_inches='tight')
    
def run_n_plot(run_path, labels, _mypars, mydata):
    '''
    run hydrus and plot results
    '''
    hydrus_dict, hyd = run_hydrus_read_results(run_path, running_hydrus=True)
    fig1, ax1, fig2, ax2 = plot_node_data(run_path,labels, _mypars, hydrus_dict)
    fig, ax = plot_ET(hydrus_dict)
    return hydrus_dict, hyd, fig1, ax1, fig2, ax2, fig, ax


def write_data2hydrus(run_path, etp, irr_data, ctop='same'):
    '''
    Write new atmosphere in data
     etp, irr_dat = input in mm. function convertes to cm
     irr_data - input positive values, function convertes to negative flux downward
    '''
    atmdf, atmtext = hh.get_atmosph_in(run_path) 
    MaxAL = len(etp) #(MaxAL  =  number  of  atmospheric  data-records)

    new_atmdf = pd.DataFrame(columns = atmdf.columns,index=list(range(MaxAL)))
    new_atmdf[:] = 0.0 
    new_atmdf['tAtm'] = atmdf['tAtm'] #range(MaxAL)
    new_atmdf.loc[0, 'tAtm'] = 0.024
    # new_atmdf['RootDepth'] = ''
    new_atmdf['hCritA'] = 10000
    if ctop == 'same':
        new_atmdf['cTop'] = atmdf['cTop']
    else:
        new_atmdf['cTop'] = ctop
        
    # Interpolate NaNs
    nans, x = nan_helper(etp)
    etp[nans]= np.interp(x(nans), x(~nans), etp[~nans])
    #################
    new_atmdf['rRoot'] =  etp /10 
    new_atmdf['Prec'] =  -irr_data /10  # divide by 10 for mm--> cm
    hh.write_atmosph_in(run_path, df=new_atmdf, mytext=atmtext)
    return atmdf, atmtext, new_atmdf

def calc_extra_water(run_path, new_atmdf, irr_data):
    '''
    calc irrigation - actual ET
    '''
    new_atmdf[['rRoot','Prec']].sum()
#    atmdf, atmtext = hh.get_atmosph_in(run_path) 
    print ('Extra Water = %.2f mm ' %(np.nansum(irr_data) - 10*new_atmdf['rRoot'].sum()) )
    return 

def plot_rdf(rdf1, z_):
    fig, ax = plt.subplots()
    fig.set_size_inches(4,8)
    ax.plot(rdf1, -z_, 'g', label='$\\beta(z)_{model}$',lw=4,alpha=0.5)
    ax.autoscale(enable=True, axis='y', tight=True)

    # ax.invert_yaxis()
    ax.set_ylabel('$z \,(cm)$')
    ax.set_xlabel('$\\beta (z)$')
    ax.legend(loc='best')
    return fig, ax
##############################################################################################################
def create_recommanded_etp(one_day_data,length=14*24, daily_mean_et=4.67,kc_s=[1.25,1.3],irrigation_frq_days = 4):
    '''
    recommanded_etp - recommanded for farmer growers
    length - how many hours or whatever time unit
    mean_et: perennial mean 
    kc1:     first week of April crop coefficient
    kc2:     first week of April crop coefficient   
    one_day_data = a normilized one day etp behavior (cos like)
    '''
    num_weeks = len(kc_s)
#    irrigation_frq_days = 4

    length += (((length - 48)/24)%irrigation_frq_days)*24    # complete irrigation cycle
    one_day_data = daily_mean_et*one_day_data
    etp = np.tile(one_day_data, int(length/24))
    # etp times Kc
    len_et = len(etp)
    etp[:int(len_et*(1/num_weeks))] *= kc_s[0]
    for i in range(1,num_weeks):
        
        etp[int(len_et*((i/num_weeks))):int(len_et*((i+1)/(num_weeks)))] *= kc_s[i]
        
    etp[int(len_et*(1-((num_weeks-1)/num_weeks))):] *=  kc_s[-1] 
    nans, x= nan_helper(etp)
    etp[nans]= etp[nans]= np.interp(x(nans), x(~nans), etp[~nans])
    return etp # in mm

def create_exact_etp(data_dict, start_date, end_date, kc_s=[1.25,1.3], irrigation_frq_days = 4):
    '''
    recommanded_etp - from forecast
    mean_et: perennial mean 
    kc1:     first week of April crop coefficient
    kc2:     first week of April crop coefficient    
    '''
    num_weeks = len(kc_s)
    tdelta = end_date - start_date
    length = (tdelta.days+1) # days in total
    delta_len = 4 - ((length) -2)%4    # length for full cycles of irrigation
    end_date += datetime.timedelta(days=delta_len)
    etp = data_dict['Besor_Farm_2013'].loc[start_date:end_date, 'ETP'].values  # divide by 10 for mm--> cm
    #    etp times Kc
    len_et = len(etp)
    
    etp[:int(len_et*(1/num_weeks))] *= kc_s[0] # first week
    
    for i in range(1,num_weeks):
        etp[int(len_et*((i/num_weeks))):int(len_et*((i+1)/(num_weeks)))] *= kc_s[i]    
        
    etp[int(len_et*(1-((num_weeks-1)/num_weeks))):] *=  kc_s[-1]   # last week

    nans, x= nan_helper(etp)
    etp[nans]= np.interp(x(nans), x(~nans), etp[~nans])
    return etp # in mm


def create_forecasted_etp(forecasted_etp, kc_s=(1.4,1.38), length=28*24):
    '''
    recommanded_etp - from forecast
    mean_et: perennial mean 
    kc1:     first week of April crop coefficient
    kc2:     first week of April crop coefficient    
    '''
    num_weeks = len(kc_s)
    etp = forecasted_etp[:length+1].copy()  # divide by 10 for mm--> cm
    len_et = len(etp)
    #    etp times Kc
    etp[:int(len_et*(1/num_weeks))] *= kc_s[0] # first week
    
    for i in range(1,num_weeks):
        etp[int(len_et*((i/num_weeks))):int(len_et*((i+1)/(num_weeks)))] *= kc_s[i]    
        
    etp[int(len_et*(1-((num_weeks-1)/num_weeks))):] *=  kc_s[-1]   # last week
    nans, x= nan_helper(etp)
    etp[nans]= np.interp(x(nans), x(~nans), etp[~nans])
    return etp # in mm

def create_irrigation(etp, irr_frq=96, length=14*24):

    '''
    create irrigation data
    irr_frq - irrigation every n timesteps
    '''
    lengthfull = len(etp) # full cycles length
    # half an irrigation +  num_irrigations
    num_irrigations = int(np.ceil(np.floor(lengthfull-(irr_frq/2))/(irr_frq)))  # 
    irr_indices_ = [range(0,4)]
    irr_indices_.extend([range(int(irr_frq*(0.5+i)), int(irr_frq*(0.5+i)+6)) for i in range(num_irrigations )])
    irr_indices_ = [list(x) for x in irr_indices_] # list range for python 3 compatibillity
    
    etp_indices_ = [range(0,int(irr_frq*0.5))]
    etp_indices_.extend([range(int(irr_frq*(i+0.5)), int(irr_frq*(i+1.5))) for i in range(num_irrigations-1)])
    etp_indices_.append(range(int(irr_frq*(1.5+(num_irrigations-2))),lengthfull))
    etp_indices_ = [list(x) for x in etp_indices_]
    
    irr_data_ = np.zeros(lengthfull)
    
    for i,indx in enumerate(irr_indices_):
        irr_data_[indx] = np.nansum(etp[etp_indices_[i]])/len(indx)
    return irr_data_[:length]

def extract_data(hydrus_dict):
    mydata  = hydrus_dict['node_out'].copy()
    hdict   = hydrus_dict['hydrus_dat'].copy()
    obsdict = hydrus_dict['obs_node_out'].copy()
    return mydata, hdict, obsdict

def define_pars(run_path, mydata, hdict):
    _mypars = list(mydata.values())[0].columns[1:].tolist()
    try:
        munit = hh.get_munit_selectorin(run_path)
    except (RuntimeError, TypeError, NameError):
        munit = 'ppm'                            
    labels = [r'$\psi \, ('+ hdict['SpaceUnit'] +')$',
              r'$\theta \,('+hdict['SpaceUnit'] +'^3\,'+hdict['SpaceUnit'] +'^{-3})$',
              r'$K \, ('+hdict['SpaceUnit'] +'\,h^{-1})$',
              r'$C$',
              r'$Flux\, ('+hdict['SpaceUnit']+'\,' +hdict['TimeUnit'] +'^{-1})$',
              r'Sink',
              r'Kappa',
              r'v/KsTop',
              r'$Temp\,(^{\circ}C)$',
              r'$Conc. \, ('+munit+' \,'+hdict['SpaceUnit'] +'^{-3})$',
              r'$orb. \, ('+munit+' \,'+hdict['SpaceUnit'] +'^{-3})$']
    return _mypars, labels