#Courtesy of P. Ryan
#Calculates integrated particle fluxes

import matplotlib.pyplot as plt
import numpy as np
import pickle
filepath='/common/uda-scratch/pryan/MU02_RT22_05_07_minTe'
bad_probes=[]

calculate = False
plotting = True
smooth_fgw = True

output_path = '/home/pryan/LP_analysis/MU02/'
output_path = '/home/kver/PycharmProjects/mu02_exhaust_scripts/LP_analysis/'

twindow = [0.4,0.8]
time_filtering = True

#####################################################################
######################################################################

if calculate:
    from mastu_exhaust_analysis.pyLangmuirProbe import LangmuirProbe, probe_array, compare_shots
    import pyuda
    client=pyuda.Client()
    from mastu_exhaust_analysis.calc_ne_bar import calc_ne_bar
    #conventional upper
    shot=[46866,46867,46868,46891,46903]
    trange=[ [0.25,0.81] ,[0.25,0.81] ,[0.25,0.81] , [0.25,0.84], [0.25,0.42]]
    output_jsat,Eich=compare_shots(output_ane=True,operation='integrate',filepath=filepath,shot=shot,combine_sectors=True,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor=['upper','upper'], sectors=[4,10], quantity = 'jsat_tile', coordinate='R',tiles=['T2','T3'],time_combine=True,show=False)
    plt.close()
    with open('CD_upper_V2.pickle','wb') as f:
        pickle.dump([output_jsat],f)

    #conventional lower
    shot=[46866,46867,46891,46903]
    trange=[ [0.25,0.81] ,[0.25,0.81] , [0.25,0.84], [0.25,0.42]]
    output_jsat,Eich=compare_shots(output_ane=True,operation='integrate',filepath=filepath,shot=shot,combine_sectors=True,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor=['lower','lower'], sectors=[4,10], quantity = 'jsat_tile', coordinate='R',tiles=['T2','T3'],time_combine=True,show=False)
    plt.close()
    with open('CD_lower_V2.pickle','wb') as f:
        pickle.dump([output_jsat],f)

    #ED lower
    shot=[47079,47082,47115]
    trange=[[0.425,0.52],[0.425,0.87],[0.45,0.87]]
    output_jsat,Eich=compare_shots(output_ane=True,operation='integrate',filepath=filepath,shot=shot,combine_sectors=False,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor='lower', sectors=4, quantity = 'jsat_tile', coordinate='R',tiles=['T2','T3','T4','T5'],time_combine=True,show=False)
    plt.close()
    with open('ED_lower_V2.pickle','wb') as f:
        pickle.dump([output_jsat],f)

    #ED upper
    shot=[47079,47082,47115]
    trange=[[0.425,0.52],[0.425,0.87],[0.45,0.87]]
    output_jsat,Eich=compare_shots(output_ane=True,operation='integrate',filepath=filepath,shot=shot,combine_sectors=True,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor=['upper','upper'], sectors=[10,4], quantity = 'jsat_tile', coordinate='R',tiles=['T2','T3','T4','T5'],time_combine=True,show=False)
    with open('ED_upper_V2.pickle','wb') as f:
        pickle.dump([output_jsat],f)




    #SXD lower
    shot=[46707,46776,46791,46792,46794,46860]
    shot = [46860]
    trange_46707=[0.44,0.8]
    trange_46776=[0.45,0.77]
    trange_46791=[0.47,0.825]
    trange_46792=[0.47,0.825]
    trange_46794=[0.45,0.825]
    trange_46860=[0.45,0.77]
    trange=[trange_46707,trange_46776,trange_46791,trange_46792,trange_46794,trange_46860]
    output_jsat,Eich=compare_shots(output_ane=True,operation='integrate',filepath=filepath,shot=shot,combine_sectors=False,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor='lower', sectors=4, quantity = 'jsat_tile', coordinate='R',tiles=['T5'],time_combine=True,show=False)
    plt.close()
    with open('SXD_lower_V2.pickle','wb') as f:
        pickle.dump([output_jsat],f)

    #SXD upper
    shot=[46707,46776,46791,46792,46794,46860]
    shot = [46860]
    trange_46707=[0.44,0.8]
    trange_46776=[0.45,0.77]
    trange_46791=[0.47,0.825]
    trange_46792=[0.47,0.825]
    trange_46794=[0.45,0.825]
    trange_46860=[0.45,0.77]
    trange=[trange_46707,trange_46776,trange_46791,trange_46792,trange_46794,trange_46860]
    output_jsat,Eich=compare_shots(output_ane=True,operation='integrate',filepath=filepath,shot=shot,combine_sectors=True,dt=0.01,Eich_fit=False,bin_x_step=1e-3,bad_probes=bad_probes,trange=trange,divertor=['upper','upper'], sectors=[10,4], quantity = 'jsat_tile', coordinate='R',tiles=['T5'],time_combine=True,show=False)
    plt.close()
    with open('SXD_upper_V2.pickle','wb') as f:
        pickle.dump([output_jsat],f)


#######################################################################
#Plotting

if plotting:

    print('Generating figure 2a')
    #Lower rollover
    CD_lower  = pickle.load(open(output_path + 'CD_lower_V2.pickle','rb'))
    CD_upper  = pickle.load(open(output_path + 'CD_upper_V2.pickle','rb'))
    ECD_lower = pickle.load(open(output_path + 'ED_lower_V2.pickle','rb'))
    ECD_upper = pickle.load(open(output_path + 'ED_upper_V2.pickle','rb'))
    SXD_lower = pickle.load(open(output_path + 'SXD_lower_V2.pickle','rb'))
    SXD_upper = pickle.load(open(output_path + 'SXD_upper_V2.pickle','rb'))

    vars=['CD_lower','CD_upper','ECD_lower','ECD_upper','SXD_lower','SXD_upper']
    #performs a polyfit for fGW and overwrites fGW
    import numpy as np
    for i in range(0,len(vars)):
        data = eval(vars[i])
        for j in range(0,len(data)):
            for k in range(0,len(data[j])):
                for l in range(0,len(data[j][k]['time'])):
                    data[j][k]['time'][l] = np.array(data[j][k]['time'][l])
                    data[j][k]['fgw_ave'][l] = np.array(data[j][k]['fgw_ave'][l])
                    data[j][k]['integrate'][l] = np.array(data[j][k]['integrate'][l])
                    time = np.array([np.mean(col) for col in data[j][k]['time'][l]])
                    if time_filtering:
                        I = np.logical_and(time>twindow[0], time<twindow[1])
                    else:
                        I = np.ones(len(time),dtype=bool)
                    if smooth_fgw:
                        p = np.polyfit(time[I],data[j][k]['fgw_ave'][l][I],2)
                        data[j][k]['fgw_ave'][l][I] = np.polyval(p,time[I])
                    data[j][k]['integrate'][l][np.logical_not(I)] = np.nan
                    data[j][k]['fgw_ave'][l][np.logical_not(I)] = np.nan

    CD_y_lower=np.array([])
    CD_y_err_lower=np.array([])
    ECD_y_lower=np.array([])
    ECD_y_err_lower=np.array([])
    SXD_y_lower=np.array([])
    SXD_y_err_lower=np.array([])
    CD_x_lower=np.array([])
    CD_x_err_lower=np.array([])
    ECD_x_lower=np.array([])
    ECD_x_err_lower=np.array([])
    SXD_x_lower=np.array([])
    SXD_x_err_lower=np.array([])

    for ii in range(1,len(CD_lower[0])):
        CD_y_lower=np.concatenate([CD_y_lower,CD_lower[0][ii]['integrate'][0]])
        CD_y_err_lower=np.concatenate([CD_y_err_lower,CD_lower[0][ii]['integrate_error'][0]])
        CD_x_lower=np.concatenate([CD_x_lower,CD_lower[0][ii]['fgw_ave'][0]])
        CD_x_err_lower=np.concatenate([CD_x_err_lower,CD_lower[0][ii]['fgw_std'][0]])

    for ii in range(len(ECD_lower[0])):
        ECD_y_lower=np.concatenate([ECD_y_lower,ECD_lower[0][ii]['integrate'][0]])
        ECD_y_err_lower=np.concatenate([ECD_y_err_lower,ECD_lower[0][ii]['integrate_error'][0]])
        ECD_x_lower=np.concatenate([ECD_x_lower,ECD_lower[0][ii]['fgw_ave'][0]])
        ECD_x_err_lower=np.concatenate([ECD_x_err_lower,ECD_lower[0][ii]['fgw_std'][0]])

    for ii in range(len(SXD_lower[0])):
        SXD_y_lower=np.concatenate([SXD_y_lower,SXD_lower[0][ii]['integrate'][0]])
        SXD_y_err_lower=np.concatenate([SXD_y_err_lower,SXD_lower[0][ii]['integrate_error'][0]])
        SXD_x_lower=np.concatenate([SXD_x_lower,SXD_lower[0][ii]['fgw_ave'][0]])
        SXD_x_err_lower=np.concatenate([SXD_x_err_lower,SXD_lower[0][ii]['fgw_std'][0]])


    CD_y_upper=np.array([])
    CD_y_err_upper=np.array([])
    ECD_y_upper=np.array([])
    ECD_y_err_upper=np.array([])
    SXD_y_upper=np.array([])
    SXD_y_err_upper=np.array([])
    CD_x_upper=np.array([])
    CD_x_err_upper=np.array([])
    ECD_x_upper=np.array([])
    ECD_x_err_upper=np.array([])
    SXD_x_upper=np.array([])
    SXD_x_err_upper=np.array([])

    for ii in range(len(CD_upper[0])):
        CD_y_upper=np.concatenate([CD_y_upper,CD_upper[0][ii]['integrate'][0]])
        CD_y_err_upper=np.concatenate([CD_y_err_upper,CD_upper[0][ii]['integrate_error'][0]])
        CD_x_upper=np.concatenate([CD_x_upper,CD_upper[0][ii]['fgw_ave'][0]])
        CD_x_err_upper=np.concatenate([CD_x_err_upper,CD_upper[0][ii]['fgw_std'][0]])

    for ii in range(len(ECD_upper[0])):
        ECD_y_upper=np.concatenate([ECD_y_upper,ECD_upper[0][ii]['integrate'][0]])
        ECD_y_err_upper=np.concatenate([ECD_y_err_upper,ECD_upper[0][ii]['integrate_error'][0]])
        ECD_x_upper=np.concatenate([ECD_x_upper,ECD_upper[0][ii]['fgw_ave'][0]])
        ECD_x_err_upper=np.concatenate([ECD_x_err_upper,ECD_upper[0][ii]['fgw_std'][0]])

    for ii in range(len(SXD_upper[0])):
        SXD_y_upper=np.concatenate([SXD_y_upper,SXD_upper[0][ii]['integrate'][0]])
        SXD_y_err_upper=np.concatenate([SXD_y_err_upper,SXD_upper[0][ii]['integrate_error'][0]])
        SXD_x_upper=np.concatenate([SXD_x_upper,SXD_upper[0][ii]['fgw_ave'][0]])
        SXD_x_err_upper=np.concatenate([SXD_x_err_upper,SXD_upper[0][ii]['fgw_std'][0]])

    x = np.concatenate([SXD_x_lower,SXD_x_upper])
    y = np.concatenate([SXD_y_lower,SXD_y_upper])
    p_SXD = np.polyfit(x[np.logical_not(np.isnan(x*y))],y[np.logical_not(np.isnan(x*y))],2)

    x = np.concatenate([ECD_x_lower,ECD_x_upper])
    y = np.concatenate([ECD_y_lower,ECD_y_upper])
    p_ECD = np.polyfit(x[np.logical_not(np.isnan(x*y))],y[np.logical_not(np.isnan(x*y))],2)

    #x = np.concatenate([CD_x_lower,CD_x_upper])
    #y = np.concatenate([CD_y_lower,CD_y_upper])
    x = np.concatenate([CD_x_lower])
    y = np.concatenate([CD_y_lower])
    p_CD = np.polyfit(x[np.logical_not(np.isnan(x*y))],y[np.logical_not(np.isnan(x*y))],2)

    polys = dict()
    polys['SXD'] = p_SXD
    polys['ECD'] = p_ECD
    polys['CD'] = p_CD
    np.save('PartFlux_poly',polys)

    pax = np.arange(0.25,0.55,0.01)
    fontsize=20
    plt.figure()
    #plt.rc('xtick', labelsize=fontsize) #fontsize of the x tick labels
    #plt.rc('ytick', labelsize=fontsize) #fontsize of the y tick labels
    plt.errorbar(x=CD_x_lower,xerr=CD_x_err_lower,y=CD_y_lower/1.6e-19/1e22,yerr=CD_y_err_lower/1.6e-19/1e22,linestyle='',marker='+',markerfacecolor='r',markeredgecolor='r',ecolor='r',label='CD lower')
    plt.errorbar(x=ECD_x_lower,xerr=ECD_x_err_lower,y=ECD_y_lower/1.6e-19/1e22,yerr=ECD_y_err_lower/1.6e-19/1e22,linestyle='',marker='o',markerfacecolor='g',markeredgecolor='g',ecolor='g',label='ED lower')
    plt.errorbar(x=SXD_x_lower,xerr=SXD_x_err_lower,y=SXD_y_lower/1.6e-19/1e22,yerr=SXD_y_err_lower/1.6e-19/1e22,linestyle='',marker='x',markerfacecolor='b',ecolor='b',markeredgecolor='b',label='SXD lower')
    #plt.xlabel('Average core density (m$^{-3}$)',fontsize=fontsize)
    plt.plot(pax,np.polyval(p_SXD,pax)/1.6e-19/1e22,'b',label='SXD fit')
    plt.plot(pax,np.polyval(p_ECD,pax)/1.6e-19/1e22,'g',label='ED fit')
    plt.plot(pax[pax>0.3],np.polyval(p_CD,pax[pax>0.3])/1.6e-19/1e22,'r',label='CD fit')
    #plt.xlabel('Greenwald Fraction')
    #plt.ylabel('Ions/s')
    plt.ylim(0,3.5)
    #plt.xlim(0,4.5e19)
    plt.xlim(0.25,0.5)
    #plt.savefig('PartFlux_integr.eps')
    #plt.grid()
    #leg=plt.legend(fontsize=fontsize)
    #leg.set_draggable(state=True)

    fontsize=20
    #plt.figure()
    #plt.rc('xtick', labelsize=fontsize) #fontsize of the x tick labels
    #plt.rc('ytick', labelsize=fontsize) #fontsize of the y tick labels
    #plt.errorbar(x=CD_x_upper,xerr=CD_x_err_upper,y=CD_y_upper/1.6e-19/1e22,yerr=CD_y_err_upper/1.6e-19/1e22,linestyle='',marker='o',markerfacecolor='r',markeredgecolor='r',ecolor='r',label='CD upper')
    plt.errorbar(x=ECD_x_upper,xerr=ECD_x_err_upper,y=ECD_y_upper/1.6e-19/1e22,yerr=ECD_y_err_upper/1.6e-19/1e22,linestyle='',marker='s',markerfacecolor='g',markeredgecolor='g',ecolor='g',label='ED upper')
    plt.errorbar(x=SXD_x_upper,xerr=SXD_x_err_upper,y=SXD_y_upper/1.6e-19/1e22,yerr=SXD_y_err_upper/1.6e-19/1e22,linestyle='',marker='v',markerfacecolor='b',ecolor='b',markeredgecolor='b',label='SXD upper')
    plt.xlabel('Greenwald Fraction',fontsize=fontsize)
    #plt.xlabel('Average core density (m$^{-3}$)',fontsize=fontsize)
    plt.ylabel('Ions/s (10$^{22}$)',fontsize=fontsize)
    import dms.general_tools as gt
    plt.plot(pax,(np.polyval(p_SXD,pax)/1.6e-19/1e22), 'b', label='SXD fit')
    plt.plot(pax,np.polyval(p_ECD,pax)/1.6e-19/1e22,'g',label='ED fit')
    plt.plot(pax[pax>0.28],np.polyval(p_CD,pax[pax>0.28])/1.6e-19/1e22,'r',label='CD fit')
    #plt.ylim(0,3.5)
    plt.xlim(0,4.5e19)
    plt.xlim(0.25,0.5)
    #plt.grid()
    #leg=plt.legend(fontsize=fontsize)
    #leg.set_draggable(state=True)
    plt.savefig('LP_integrated.eps')
    plt.show()

print('Done')
