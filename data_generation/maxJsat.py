#courtesy of P. Ryan
#Calculates peak particle fluxes
#This requires LP analysis tools that are embedded in the UKAEA exhaust analysis tools git repo (git.ccfe.ac.uk/jrh/mastu_exhaust_analysis), DOI: 10.1063/5.0152680.

import matplotlib.pyplot as plt
import numpy as np
import pickle
filepath='/common/uda-scratch/pryan/MU02_RT22_05_07_minTe'
bad_probes=[]
output_path = '/home/kver/PycharmProjects/mu02_exhaust_scripts/LP_analysis/'

calculate = False
plotting = True

if calculate:
    from mastu_exhaust_analysis.pyLangmuirProbe import LangmuirProbe, probe_array, compare_shots
    import pyuda

    client = pyuda.Client()
    from mastu_exhaust_analysis.calc_ne_bar import calc_ne_bar
    ##################################################
    #find the max value of jsat.
    ane_ave=[]
    ane_std=[]
    fgw_ave=[]
    fgw_std=[]
    Te=[]
    Te_error=[]
    jsat=[]
    jsat_error=[]
    for shot_index in range(len(output_jsat)):
        for time_index in range(len(output_jsat[shot_index]['x'][0])):
            ane_ave.append(output_jsat[shot_index]['ane_ave'][0][time_index])
            ane_std.append(output_jsat[shot_index]['ane_std'][0][time_index])
            fgw_ave.append(output_jsat[shot_index]['fgw_ave'][0][time_index])
            fgw_std.append(output_jsat[shot_index]['fgw_std'][0][time_index])
            max_index=np.nanargmax(output_jsat[shot_index]['y'][0][time_index])
            Te.append(output_Te[shot_index]['y'][0][time_index][max_index])
            Te_error.append(output_Te[shot_index]['y_error'][0][time_index][max_index])
            jsat.append(output_jsat[shot_index]['y'][0][time_index][max_index])
            jsat_error.append(output_jsat[shot_index]['y_error'][0][time_index][max_index])
        #    print(output_jsat[shot_index]['fgw_ave'][0][time_index])
        #    print(output_Te[shot_index]['y'][0][time_index][max_index])
            #print(shot_index)
        #    print(output_Te[shot_index]['time_ave'][0][time_index])
        #    print(output_Te[shot_index]['x'][0][time_index][max_index])
        #    print('next')

    ane_ave=np.array(ane_ave)
    ane_std=np.array(ane_std)
    fgw_ave=np.array(fgw_ave)
    fgw_std=np.array(fgw_std)
    Te=np.array(Te)
    Te_error=np.array(Te_error)
    jsat=np.array(jsat)
    jsat_error=np.array(jsat_error)

    plt.figure()
    plt.errorbar(x=fgw_ave,xerr=fgw_std,y=Te,yerr=Te_error,marker='o',linestyle='')
    plt.figure()
    plt.errorbar(x=fgw_ave,xerr=fgw_std,y=jsat,yerr=jsat_error,marker='o',linestyle='')
    plt.show()

    store={}
    store['ane_ave']=ane_ave
    store['ane_std']=ane_std
    store['fgw_ave']=fgw_ave
    store['fgw_std']=fgw_std
    store['Te']=Te
    store['Te_error']=Te_error
    store['jsat']=jsat
    store['jsat_error']=jsat_error

    #manually change the file name
    with open('SXD_lower_maxJsatTile_Te_V2.pickle','wb') as f:
        pickle.dump([store],f)

    #####################################################################
    ######################################################################
    #collect data below

    #conventional upper
    shot=[46866,46867,46868,46891,46903]
    trange=[ [0.25,0.81] ,[0.25,0.81] ,[0.25,0.81] , [0.25,0.84], [0.25,0.42]]
    output_jsat,Eich=compare_shots(psi_n_range=[0.96,1.1],output_ane=True,operation='profile',filepath=filepath,shot=shot,combine_sectors=True,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor=['upper','upper'], sectors=[4,10], quantity = 'jsat_tile', coordinate='R',tiles=['T2','T3'],time_combine=True,show=False)
    plt.close()
    output_Te,Eich=compare_shots(psi_n_range=[0.96,1.1],output_ane=False,operation='profile',filepath=filepath,shot=shot,combine_sectors=True,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor=['upper','upper'], sectors=[4,10], quantity = 'Te', coordinate='R',tiles=['T2','T3'],time_combine=True,show=False)
    plt.close()


    #conventional lower
    shot=[46866,46867,46891,46903]
    trange=[ [0.25,0.81] ,[0.25,0.81] , [0.25,0.84], [0.25,0.42]]
    output_jsat,Eich=compare_shots(psi_n_range=[0.96,1.1],output_ane=True,operation='profile',filepath=filepath,shot=shot,combine_sectors=True,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor=['lower','lower'], sectors=[4,10], quantity = 'jsat_tile', coordinate='R',tiles=['T2','T3'],time_combine=True,show=False)
    plt.close()
    output_Te,Eich=compare_shots(psi_n_range=[0.96,1.1],output_ane=False,operation='profile',filepath=filepath,shot=shot,combine_sectors=True,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor=['lower','lower'], sectors=[4,10], quantity = 'Te', coordinate='R',tiles=['T2','T3'],time_combine=True,show=False)
    plt.close()

    #ED lower
    shot=[47079,47082,47115]
    trange=[[0.425,0.52],[0.425,0.87],[0.45,0.87]]
    output_jsat,Eich=compare_shots(psi_n_range=[0.96,1.1],output_ane=True,operation='profile',filepath=filepath,shot=shot,combine_sectors=False,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor='lower', sectors=4, quantity = 'jsat_tile', coordinate='R',tiles=['T2','T3','T4','T5'],time_combine=True,show=False)
    plt.close()
    output_Te,Eich=compare_shots(psi_n_range=[0.96,1.1],output_ane=False,operation='profile',filepath=filepath,shot=shot,combine_sectors=False,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor='lower', sectors=4, quantity = 'Te', coordinate='R',tiles=['T2','T3','T4','T5'],time_combine=True,show=False)
    plt.close()


    #ED upper
    shot=[47079,47082,47115]
    trange=[[0.425,0.52],[0.425,0.87],[0.45,0.87]]
    output_jsat,Eich=compare_shots(max_y_error_remove=1000,psi_n_range=[0.96,1.1],output_ane=True,operation='profile',filepath=filepath,shot=shot,combine_sectors=True,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor=['upper','upper'], sectors=[10,4], quantity = 'jsat_tile', coordinate='R',tiles=['T2','T3','T4','T5'],time_combine=True,show=False)
    plt.close()
    output_Te,Eich=compare_shots(max_y_error_remove=15,psi_n_range=[0.96,1.1],output_ane=False,operation='profile',filepath=filepath,shot=shot,combine_sectors=True,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor=['upper','upper'], sectors=[10,4], quantity = 'Te', coordinate='R',tiles=['T2','T3','T4','T5'],time_combine=True,show=False)
    plt.close()




    #SXD lower
    shot=[46707,46776,46791,46792,46794,46860]
    trange_46707=[0.44,0.8]
    trange_46776=[0.45,0.77]
    trange_46791=[0.47,0.825]
    trange_46792=[0.47,0.825]
    trange_46794=[0.45,0.825]
    trange_46860=[0.45,0.77]
    trange=[trange_46707,trange_46776,trange_46791,trange_46792,trange_46794,trange_46860]
    output_jsat,Eich=compare_shots(psi_n_range=[0.96,1.1],output_ane=True,operation='profile',filepath=filepath,shot=shot,combine_sectors=False,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor='lower', sectors=4, quantity = 'jsat_tile', coordinate='R',tiles=['T5'],time_combine=True,show=False)
    plt.close()
    output_Te,Eich=compare_shots(psi_n_range=[0.96,1.1],output_ane=False,operation='profile',filepath=filepath,shot=shot,combine_sectors=False,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor='lower', sectors=4, quantity = 'Te', coordinate='R',tiles=['T5'],time_combine=True,show=False)
    plt.close()


    #SXD upper
    shot=[46707,46776,46791,46792,46794,46860]
    trange_46707=[0.44,0.8]
    trange_46776=[0.45,0.77]
    trange_46791=[0.47,0.825]
    trange_46792=[0.47,0.825]
    trange_46794=[0.45,0.825]
    trange_46860=[0.45,0.77]
    trange=[trange_46707,trange_46776,trange_46791,trange_46792,trange_46794,trange_46860]
    output_jsat,Eich=compare_shots(psi_n_range=[0.96,1.1],output_ane=True,operation='profile',filepath=filepath,shot=shot,combine_sectors=False,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor='upper', sectors=4, quantity = 'jsat_tile', coordinate='R',tiles=['T5'],time_combine=True,show=False)
    plt.close()
    output_Te,Eich=compare_shots(psi_n_range=[0.96,1.1],output_ane=False,operation='profile',filepath=filepath,shot=shot,combine_sectors=False,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor='upper', sectors=4, quantity = 'Te', coordinate='R',tiles=['T5'],time_combine=True,show=False)
    plt.close()

if plotting:
    #######################################################################
    #Plotting

    CD_lower  = pickle.load(open(output_path + 'CD_lower_maxJsatTile_Te_V2.pickle','rb'))[0]
    CD_upper  = pickle.load(open(output_path + 'CD_upper_maxJsatTile_Te_V2.pickle','rb'))[0]
    ED_lower = pickle.load(open(output_path + 'ED_lower_maxJsatTile_Te_V2.pickle','rb'))[0]
    ED_upper = pickle.load(open(output_path + 'ED_upper_maxJsatTile_Te_V2.pickle','rb'))[0]
    SXD_lower = pickle.load(open(output_path + 'SXD_lower_maxJsatTile_Te_V2.pickle','rb'))[0]
    SXD_upper = pickle.load(open(output_path + 'SXD_upper_maxJsatTile_Te_V2.pickle','rb'))[0]

    x = np.concatenate([SXD_lower['fgw_ave'],SXD_upper['fgw_ave']])
    y = np.concatenate([SXD_lower['jsat'],SXD_upper['jsat']])
    p_SXD = np.polyfit(x[np.logical_not(np.isnan(x*y))],y[np.logical_not(np.isnan(x*y))],2)

    x = np.concatenate([ED_lower['fgw_ave'],ED_upper['fgw_ave']])
    y = np.concatenate([ED_lower['jsat'],ED_upper['jsat']])
    p_ECD = np.polyfit(x[np.logical_not(np.isnan(x*y))],y[np.logical_not(np.isnan(x*y))],2)

    x = np.concatenate([CD_lower['fgw_ave'],CD_upper['fgw_ave']])
    y = np.concatenate([CD_lower['jsat'],CD_upper['jsat']])
    p_CD = np.polyfit(x[np.logical_not(np.isnan(x*y))],y[np.logical_not(np.isnan(x*y))],2)

    polys = dict()
    polys['SXD'] = p_SXD
    polys['ECD'] = p_ECD
    polys['CD'] = p_CD
    np.save('PartFluxPeak_poly', polys)

    fontsize=20
    plt.figure()
    #plt.rc('xtick', labelsize=fontsize) #fontsize of the x tick labels
    #plt.rc('ytick', labelsize=fontsize) #fontsize of the y tick labels
    plt.errorbar(x=CD_lower['fgw_ave'],xerr=CD_lower['fgw_std'],y=CD_lower['jsat']/1.602e-19,yerr=CD_lower['jsat_error']/1.602e-19,linestyle='',marker='+',markerfacecolor='r',markeredgecolor='r',ecolor='r',label='CD lower')
    plt.errorbar(x=ED_lower['fgw_ave'],xerr=ED_lower['fgw_std'],y=ED_lower['jsat']/1.602e-19,yerr=ED_lower['jsat_error']/1.602e-19,linestyle='',marker='o',markerfacecolor='g',markeredgecolor='g',ecolor='g',label='ED lower')
    plt.errorbar(x=SXD_lower['fgw_ave'],xerr=SXD_lower['fgw_std'],y=SXD_lower['jsat']/1.602e-19,yerr=SXD_lower['jsat_error']/1.602e-19,linestyle='',marker='x',markerfacecolor='b',ecolor='b',markeredgecolor='b',label='SXD lower')
    #plt.xlabel('Average core density (m$^{-3}$)',fontsize=fontsize)
    plt.xlabel('Greenwald Fraction',fontsize=fontsize)
    plt.ylabel('Max $j^+_{tile}$ (kAm$^{-2}$)',fontsize=fontsize)
    plt.yscale('log')
    plt.ylim(3e21,2e23)
    #plt.xlim(0,4.5e19)
    plt.xlim(0.25,0.5)
    plt.grid()
    plt.savefig('Peak_part_flux_log.eps')
    #leg=plt.legend(fontsize=fontsize)
    #leg.set_draggable(state=True)


print('done')
