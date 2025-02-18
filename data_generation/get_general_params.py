#Retrieves supporting data needed for figures (Thomson, Greenwald fraction, PSOL) and saves them as local files
#This requires the UKAEA exhaust analysis tools git repo (git.ccfe.ac.uk/jrh/mastu_exhaust_analysis), DOI: 10.1088/1361-6587/ad4058.

import numpy as np
shotnrs = [46860, 46895, 47115, 46866, 46705, 46904, 47085, 47079, 46762, 46707, 46769, 46776, 46791, 46792, 46794, 47082, 47115, 47118, 46762, 46864, 46867, 46868, 46891, 46903,48008]
twindow = [[0.4, 0.8], [0.25, 0.9], [0.3, 0.95], [0.3, 0.9], [0.4, 0.95], [0.3, 0.85], [0.3, 0.67], [0.3, 0.8], [0.4, 0.8], [0.4, 0.8], [0.4, 0.8], [0.4, 0.8], [0.4, 0.8], [0.4, 0.8], [0.4, 0.8], [0.4, 0.8], [0.4, 0.8], [0.4, 0.8], [0.4, 0.8], [0.4, 0.8], [0.4, 0.8], [0.4, 0.8], [0.4, 0.8], [0.4, 0.8],[0.4,0.8]]

t_window_f = [0.4, 0.8]

get_fgw = True
get_Psep = True
get_pn = True
get_thomson = True


if get_fgw:
    #get core Greenwald fractions (with polynomial fits)
    import mastu_exhaust_analysis as mea
    print('Retrieving core Greenwald fractions')
    fgw = dict()
    for i in range(0,len(shotnrs)):
        dum = mea.calc_ne_bar(shotnrs[i])
        p = np.polyfit(dum['t'][np.logical_and(dum['t']>t_window_f[0], dum['t']<t_window_f[1])],dum['greenwald_fraction'][np.logical_and(dum['t']>t_window_f[0], dum['t']<t_window_f[1])],2)
        fgw[str(shotnrs[i])] = dum
        fgw[str(shotnrs[i])]['p'] = p
    np.save('ngw_info',fgw)

if get_Psep:
    #get Psep estimates
    import mastu_exhaust_analysis as mea
    print('Retrieving Psol')
    fgw = dict()
    for i in range(0,len(shotnrs)):
        try:
            dum = mea.calc_psol(shotnrs[i],smooth_dt=0.075)
            p = np.polyfit(dum['t'][np.logical_and(dum['t']>t_window_f[0], dum['t']<t_window_f[1])],dum['psol'][np.logical_and(dum['t']>t_window_f[0], dum['t']<t_window_f[1])],2)
            fgw[str(shotnrs[i])] = dum
            fgw[str(shotnrs[i])]['p'] = p
        except:
            print('')
    np.save('psol_info',fgw)

if get_pn:
    #get neutral pressures
    import pyuda
    c = pyuda.Client()
    data = dict()
    for i in range(0, len(shotnrs)):
        b_midplane = True
        b_lower_divertor = True
        try:
            FIG_midplane = c.get('aga/hm12', source=shotnrs[i])  # midplane FIG
        except:
            b_midplane = False
        try:
            FIG_lower_divertor = c.get('aga/hl11', source=shotnrs[i])  # divertor FIG
        except:
            b_lower_divertor = False
        data[str(shotnrs[i])] = dict()
        if b_lower_divertor:
            data[str(shotnrs[i])]['t_l_div']=FIG_lower_divertor.time.data
            data[str(shotnrs[i])]['pn_l_div']=FIG_lower_divertor.data
        else:
            data[str(shotnrs[i])]['t_l_div']=np.nan
            data[str(shotnrs[i])]['pn_l_div']=np.nan
        if b_midplane:
            data[str(shotnrs[i])]['t_mid']=FIG_midplane.time.data
            data[str(shotnrs[i])]['pn_mid']=FIG_midplane.data
        else:
            data[str(shotnrs[i])]['t_mid']=np.nan
            data[str(shotnrs[i])]['pn_mid']=np.nan

    np.save('pn_info', data)

if get_thomson:
    from mastu_exhaust_analysis.pyThomson import Thomson

    for i in range(0,len(shotnrs)):
        ts_data = Thomson(shot = shotnrs[i])
        ts_data.save_to_matfile('TS_'+str(shotnrs[i])+'.mat',trange=twindow[i])
