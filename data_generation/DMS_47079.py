#script for processing the spectroscopic data for #47079, using the DMS diagnostic and repository (https://git.ccfe.ac.uk/dms_development/dms.git), hosted on the UKAEA git server, DOI: 10.1088/1741-4326/aca10a, 10.1088/1741-4326/acf946, 10.1088/1741-4326/ad5851.

import numpy as np
import dms.analysis.line_integrator as line_integrator
import dms.general_tools as gt
import dms.data_read_gui.data_loading as dl
import numpy as np
import matplotlib.pyplot as plt
import dms.analysis.fit.fit_discharge as fd
from scipy.interpolate import interp1d
import mat73
from datetime import datetime

shotnrs = [47079]

#first prepare input file
shotnr = shotnrs[0]
inputfile = '/home/kver/'+str(shotnr)+'_BaSPMI.npy' #input data
fdir_temp = '/hdd/'+str(shotnr)+'_BaySPMI_low_Te_160623/'  #storage of result
fitfile = '/common/projects/diagnostics/MAST/SPEXBDMS/analysis_results/' + str(shotnr) + '/n=6 Stark mc.npz'
Gfile_calc = '/home/kver/G_46762_1.npy'

compute_output = False
prepare_input = True
calc_line_int = False
calculate_geom = False
compute_rates = True

indx_SP = 10

if prepare_input:
    # first scrape data for default analysis

    # performs default line integration
    if calc_line_int:
        line_integrator.auto_analysis_default(shotnrs, default_path='/home/kver/DMS_auto/')

    # visualise data
    o = np.load('/home/kver/DMS_auto/' + str(shotnr) + '/dms_auto_line_integral_' + str(shotnr) +'.npy', allow_pickle=True)

    n5 = o[1]['value'] * 4 * np.pi * (1 / 0.9)  # take out sterradian and include window transmission (90%)
    n6 = o[0]['value'] * 4 * np.pi * (1 / 0.9)
    Da = o[3]['value'] * 4 * np.pi * (1 / 0.9)
    Fulcher = o[2]['value'] * 4 * np.pi * (1 / 0.9)
    t = o[1]['axis_value'][0]

    # take out core line of sight
    n5 = n5[:,0:30]
    n6 = n6[:,0:30]
    Da = Da[:,0:30]
    Fulcher = Fulcher[:,0:30]

    # resolve spatial misalignment
    I = np.append(np.arange(0, 10, 0.5), np.arange(10, 20) - 3)
    Ireq = np.append(np.arange(0, 10, 0.5), np.linspace(10,16,10))

    # take out broken lines of sight
    Da[:,[18,20]] = np.nan
    Fulcher[:,[18,20]] = np.nan
    n5[:, [18,19]] = np.nan
    n6[:, [18,19]] = np.nan

    # correct misalignment
    n6 = np.transpose(gt.fix_alignment(I,Ireq,np.transpose(n6)))
    n5 = np.transpose(gt.fix_alignment(I,Ireq,np.transpose(n5)))
    Da = np.transpose(gt.fix_alignment(I,Ireq,np.transpose(Da)))
    Fulcher = np.transpose(gt.fix_alignment(I,Ireq,np.transpose(Fulcher)))

    # adjust timebase
    f = interp1d(o[2]['axis_value'][0],np.transpose(Fulcher),bounds_error=False,fill_value='extrapolate')
    Fulcher = np.transpose(f(t))
    f = interp1d(o[2]['axis_value'][0],np.transpose(Da),bounds_error=False,fill_value='extrapolate')
    Da = np.transpose(f(t))

    LRS = gt.smooth2a(n6 / n5, x_stddev=0.65, y_stddev=0.65)

    # prevent issues at boundary from smoothing
    LRS[:, -1] = n6[:, -1] / n5[:, -1]
    LRS[:, 0] = n6[:, 0] / n5[:, 0]

    # obtain electron density (45371)
    time_min = 0.4
    time_max = 0.85

    nemax = 4e19
    nemin = 1e19

    # load Stark broadening results
    fit = gt.fixload(np.load(fitfile, allow_pickle=True))

    # extract ne and filter

    neL = fit['output']['result'][1]['nestark_values'][:, :, 0]
    neH = fit['output']['result'][1]['nestark_values'][:, :, 2]
    ne = fit['output']['result'][1]['nestark_values'][:, :, 1]

    RelErr = (neH - neL) / ne

    neminR = 8e18
    RelErrMax = 1

    ne[ne > nemax] = np.nan
    ne[ne < nemin] = np.nan
    ne[RelErr > RelErrMax] = np.nan
    ne[:, [18,19]] = np.nan
    ne = ne [:,:30]

    ne = gt.inpaint_nans(ne)  # interpolate over filtered results
    ne = gt.smooth2a(ne)
    # correct misalignment
    ne = np.transpose(gt.fix_alignment(I,Ireq,np.transpose(ne)))

    ne = gt.inpaint_nans(ne)  # interpolate over filtered results
    ne = gt.smooth2a(ne)
    tstark = fit['output']['time']

    # resample density inferences to analyse grid by interpolating between them
    ne_mat = np.zeros(np.shape(LRS))
    for i in range(0, 30):
        f = interp1d(tstark, ne[:, i])
        ne_mat[:, i] = f(t)

    # smooth density profile
    ne_mat[np.logical_and(t > time_min, t < time_max), :] = gt.smooth2a(
        ne_mat[np.logical_and(t > time_min, t < time_max), :], x_stddev=0.8, y_stddev=0.4)

    # prepare output for analysis

    input = dict()
    input['AbsErr'] = 0.125  # Absolute uncertainty (68% confidence)
    input['RelErr'] = 0.075  # Relative uncertainty (line ratios - 68% confidence)
    input['Iter'] = 500  # Monte Carlo points
    input['FRec_reject'] = 1.3  # Reject values where FRec goes beyond a certain fraction of the FRec regime
    input['FRecForceTrend'] = 1  # Force FRec trend 0 - no, 1 - yes (this implicity assumes that the discharge analysed is a density ramp discharge where the electron temperature gradually decreases (and thus the EIR fraction of the emission is either constant or increasing)
    input['DoParFor'] = 1  # parallel computation
    input['DenMin'] = 7.5e18  # minimum electron density
    input['DenMax'] = 5e19  # maximum electron density
    input['TeE_filter'] = 0  # Use the excitation emission derived temperature as a filter 0 - no, 1 - yes (this implicity assumes that the discharge analysed is a density ramp discharge where the electron temperature gradually decreases (and thus the EIR fraction of the emission is either constant or increasing)
    input['TeR_filter'] = 0  # Filters outcomes where the EIR brightness exceed that which can be obtained from ADAS (e.g. TeR < 0.2 eV) in that case 0.2 eV temperature is assumed - 0 - no, 1 - yes
    input['CoeffUncertainty'] = 0.1  # uniform uncertainty in all atomic emission coefficients (the assumed uncertainty in all the molecular coefficients is double the atomic ones)
    input['Den'] = np.transpose(ne_mat)  # electron density
    input['DenErr'] = 3e19 * np.ones(np.shape(input['Den']))  # uncertainty electron density (68 % confidence interval)
    input['n1Int'] = np.transpose(n5)  # n1 brightness (lowest-n medium-n Balmer line used)
    input['n1Int'][:, t < time_min] = np.nan  # filter time
    input['n1Int'][:, t > time_max] = np.nan  # filter time
    input['n2Int'] = input['n1Int'] * np.transpose(LRS)  # n2 brightness (highest-n medium-n Balmer line used) computed by multiplying line ratio with n1 brightness
    input['none'] = np.ones(np.shape(input['n1Int']))  # neutral fraction (not used - obsolete)
    input['n1'] = 5  # n1 Balmer line index
    input['n2'] = 6  # n2 Balmer line index
    input['DaMea'] = np.transpose(Da)  # measured Da
    input['nStark'] = 6  # Balmer line index of Stark broadening (not used - for reference)
    input['none_loguniform'] = 1  # 1 - neutral fraction is sampled using a log-uniform distribution; 0 - uniform distribution is used
    input['pMolG'] = [-1.9835,17.5]  # polynomial (in log-log space) of DL*nH2 as function of Te (used to estimate the D2 contribution to the Balmer line emisison - NOT the D2+ & D- contributions)
    input['shot'] = shotnr  # discharge number (for reference)
    input['sys'] = 1  # system number (for reference)
    input['Time'] = t  # time vector

    # perform magnetic mapping of LoS with plasma geometry. Calculate line-of-sight intersections with the separatrix, leading to R, Z positions. Estimate pathlengths of emission region along line of sight
    import dms.analysis.geometry as geom

    input['R'] = np.zeros(np.shape(input['n1Int']))
    input['Z'] = np.zeros(np.shape(input['n1Int']))
    input['DL'] = np.zeros(np.shape(input['n1Int']))

    G = np.load('/home/kver/G_' + str(shotnr) + '_1.npy', allow_pickle=True)
    G=G[()]

    input['R'] = np.zeros(np.shape(input['n1Int']))
    input['Z'] = np.zeros(np.shape(input['n1Int']))
    input['DL'] = np.zeros(np.shape(input['n1Int']))

    # map timebase to same timebase as analysis through interpolation
    for i in range(0, np.shape(n5)[1]):
        f = interp1d(G['time'], G['R'][:, i])
        input['R'][i, :] = f(t)  # R position of intersection
        f = interp1d(G['time'], G['Z'][:, i])
        input['Z'][i, :] = f(t)  # Z position of intersection
        f = interp1d(G['time'], G['DL'][:, i])
        input['DL'][i, :] = f(t)  # pathlength estimated

    input['DLL'] = input['DL'] / 3  # Lower pathlength estimate (68% confidence). Pathlengths are sampled according to asymmetric Gaussian distribution with witdths DL-DLL and DLH-DL.
    input['DLH'] = input['DL'] * 1.3  # Upper pathlength estimate (68% confidence). Pathlengths are sampled according to asymmetric Gaussian distribution with witdths DL-DLL and DLH-DL.
    input['Fulcher'] = np.transpose(Fulcher) * input['n1Int'] / input['n1Int']  # Fulcher brightness (not used in matlab routine - for reference)
    input['fname_fulcher'] = '/home/kver/PycharmProjects/dms_nf_letter/Fulch_constr.npy'  # directory of Fulcher band constraints (not used in Matlab routine - for reference)

    # more rigorous none estimation
    Te_max = 10
    Te_min = 3
    noneDL_scaling = 10 ** (-1.4721 - 1.126 * np.log10([Te_min, Te_max]))  # MAST-U
    DL_charac = 0.15
    ErrFrac = 5
    none_min = (noneDL_scaling[1] / DL_charac) / ErrFrac
    none_max = (noneDL_scaling[0] / DL_charac) * ErrFrac
    input['noneH'] = none_max * np.ones(np.shape(input['n1Int']))  # maximum neutrla fraction (used)
    input['noneL'] = none_min * np.ones(np.shape(input['n1Int']))  # minimum neutral fraction (used)

    ne_charac = 1e19
    input['pMolG'] = [-2.3075, -1.176 + np.log10(ne_charac) + np.log10(0.15)]  # more accurate pmolG scaling

    input['D2pDm_model'] = 2  # Model used to detemrine D2+/D- separation of Balmer line emission (only used if no Db emission provided - 2 implies AMJUEL is used together with Janev, Reiter, 2018 mapping from H to D)

    input['AdasLowTeEIR'] = 1  # low Te ADAS switch for EIR calculation

    input['settings'] = dict()
    input['settings']['fname'] = '/home/kver/BaSPMI_'+str(shotnr)+'_lowTeEIR_v4/'  # filename of output
    input['settings']['starting_contamin'] = dict()
    input['settings']['starting_contamin']['n1'] = 0.5 * np.ones(
        np.shape(input['n1Int']))  # assumed initial condition of n1 molecular contamination (by default 0)
    LR_norm = (input['n2Int'] / input['n1Int'] - 0.16) / (np.nanmax(input['n2Int'] / input['n1Int']) - 0.16)
    LR_norm[-1, :] = np.nan
    LR_norm = gt.inpaint_nans(LR_norm)
    LR_norm_accum = np.maximum.accumulate(np.maximum.accumulate(LR_norm, 1)[::-1], 0)[::-1]
    input['settings']['starting_contamin']['n2'] = input['settings']['starting_contamin']['n1'] * (1 - (1 - (
                0.5 / 2.1)) * LR_norm_accum)  # assumed initial condition of n2 molecular contamination (by default 0). Here this is calculated such that the atomic only line ratios reach high values at th eend of the discharge (0.5)
    import scipy.io as sio

    input['FulcherMax'] = 1.7e19

    sio.savemat('/home/kver/'+str(shotnr)+'_BaSPMI.mat', {'input': input})  # save as matlab file
    np.save('/home/kver/'+str(shotnr)+'_BaSPMI.npy', input)  # save as numpy file

selection_t = np.arange(np.argwhere(np.nansum(input['n1Int'],axis=0)>0)[0], np.argwhere(np.nansum(input['n1Int'],axis=0)>0)[-1])  # time indexes to sample
selection_l = np.arange(0, 30)  # line of sight indeces to analyse

if compute_output:
    #analyse input file
    inputfile = '/home/kver/'+str(shotnr)+'_BaSPMI.npy'
    A = np.load(inputfile, allow_pickle=True)
    input = A[()]
    # #runs analysis and performs sampling
     # temporary storage directory
    import dms.general_tools as gt

    gt.mkdir(fdir_temp)
    import shutil

    import dms.analysis.emission.Balmer_analysis_bayes as bay
    import dms.analysis.emission.Balmer_analysis as bal
    input['Iter'] = 500
    bay.run_analysis(input, selection_t=selection_t, selection_l=selection_l, fdir_temp=fdir_temp, low_te_EIR=True,
                 not_only_sampling=True,data_reduction=True)
    np.save(fdir_temp + 'input.npy', input)

G = np.load(Gfile_calc, allow_pickle=True)
G=G[()]
offset = np.nanmedian(G['L'][np.logical_and(G['time']>0.4, G['time']<0.8),indx_SP])

if compute_rates:

    selection_l = np.arange(indx_SP,30)

    import dms.analysis.emission.Balmer_analysis_bayes as bay
    import dms.analysis.emission.Balmer_analysis as bal
    #get physics output
    output = bay.rates_extrap(fdir_temp, selection_t, selection_l) #rates & uncertainty propagation
    output = bal.rates_integr(output,select_l=selection_l) #integrate rates
    np.save(fdir_temp + 'output_proc.npy', output)
