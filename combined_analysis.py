#Script for performing a combined analysis given the output from the data generation routines.
#
# Generates figures 1d, 2, 3, 4, 5g-i (experimental), 6 (extended data)
#
#This uses information from:
#1) Spectroscopy analysis (performed by scrips DMS_$SHOTNUMBER$.py
#      fdir_SXD, fdir_ED, fdir_CD (analysis results), Gfile_SXD, Gfile_CD, Gfile_ED (geometrical mapping files)
#2) Particle flux information, obtained from integrate_jsat.py and maxJsat.py
#       file_partflux, file_peakpartflux
#3) General parameters obtained from get_general_params.py
#       file_coredens (Greenwald fraction), file_psol (PSOL), file_pn (Neutral pressure), file_TS_... (TS data)

#This requires the DMS data processing routines https://git.ccfe.ac.uk/dms_development/dms.git

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsqr
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter
from scipy.integrate import cumtrapz
from scipy.io import loadmat

#set up file locations
fdir_SXD = '/hdd/46860_BaySPMI_low_Te_150623/'  #storage of result
Gfile_SXD = '/home/kver/G_46860_1.npy'

fdir_CD = '/hdd/46762_BaySPMI_low_Te_150523/'  #storage of result
Gfile_CD = '/home/kver/G_46762_1.npy'

fdir_ED = '/hdd/47079_BaySPMI_low_Te_160623/'  #storage of result
Gfile_ED = '/home/kver/G_47079_1.npy'

file_partflux = '/home/kver/PycharmProjects/mu02_exhaust_scripts/PartFlux_poly.npy'
file_peakpartflux = '/home/kver/PycharmProjects/mu02_exhaust_scripts/PartFluxPeak_poly.npy'

file_coredens = '/home/kver/PycharmProjects/mu02_exhaust_scripts/ngw_info.npy'
file_psol = '/home/kver/PycharmProjects/mu02_exhaust_scripts/psol_info.npy'
file_pn = '/home/kver/PycharmProjects/mu02_exhaust_scripts/pn_info.npy'
file_TS_SXD = '/home/kver/PycharmProjects/mu02_exhaust_scripts/TS_values/TS_46707.mat'
file_TS_CD = '/home/kver/PycharmProjects/mu02_exhaust_scripts/TS_values/TS_46866.mat'
file_TS_ED = '/home/kver/PycharmProjects/mu02_exhaust_scripts/TS_values/TS_47079.mat'
file_TS_SPscan = '/home/kver/PycharmProjects/mu02_exhaust_scripts/TS_values/TS_46895.mat'



#load data - DMS & mappings
output_SXD = np.load(fdir_SXD + 'output_proc.npy',allow_pickle=True)[()]
G_SXD = np.load(Gfile_SXD, allow_pickle=True)[()]

output_ED = np.load(fdir_ED + 'output_proc.npy',allow_pickle=True)[()]
G_ED = np.load(Gfile_ED, allow_pickle=True)[()]

output_CD = np.load(fdir_CD + 'output_proc.npy',allow_pickle=True)[()]
G_CD = np.load(Gfile_CD, allow_pickle=True)[()]

#load data - core densities
p_dens = np.load(file_coredens,allow_pickle=True)[()]

#load data - particle flux fits
p_partflux = np.load(file_partflux,allow_pickle=True)[()]
p_peakpartflux = np.load(file_peakpartflux,allow_pickle=True)[()]

Figure4=True
Figure6=True
Figure2=True
Figure5=True
Figure3=True
Ploss_PMI=True

#calculations (prior)

den_SXD = np.polyval(p_dens['46860']['p'], output_SXD['input']['Time'])
den_SXD[np.logical_or(output_SXD['input']['Time']<0.4, output_SXD['input']['Time']>0.78)] = np.nan
partflux_SXD = np.polyval(p_partflux['SXD'], den_SXD)
den_CD = np.polyval(p_dens['46762']['p'], output_CD['input']['Time'])
den_CD[np.logical_or(output_CD['input']['Time']<0.3, output_CD['input']['Time']>0.8)] = np.nan
partflux_CD = np.polyval(p_partflux['CD'], den_CD)
den_ED = np.polyval(p_dens['47079']['p'], output_ED['input']['Time'])
den_ED[np.logical_or(output_ED['input']['Time']<0.35, output_ED['input']['Time']>0.8)] = np.nan
partflux_ED = np.polyval(p_partflux['ECD'], den_ED)

#estimate target power loading
gamma = 7
epsilon = 13.6+2.2
ev = 1.602e-19

#Strike point position estimates with respect to DMS view
indx_SP_SXD = 2
indx_SP_CD = 27
indx_SP_ED = 12

#Average emission weighted temperature must assume a characteristic temperature for MAR, MAD, MAI, assume that this characteristic temperature is any random temperature between the characteristic temperature of electron-impact excitation and electron-ion recombination
#calculate power loads
R = np.random.rand(output_SXD['input']['Iter'])
Te_avg_SXD = (output_SXD['ResultMC']['TeEMC']*output_SXD['ResultMC']['DaEIE'] + output_SXD['ResultMC']['TeRMC']*output_SXD['ResultMC']['DaEIR'] + output_SXD['ResultMC']['DaMOL']*(output_SXD['ResultMC']['TeRMC'] + R*(output_SXD['ResultMC']['TeEMC']-output_SXD['ResultMC']['TeRMC'])))/(output_SXD['ResultMC']['DaMOL'] + output_SXD['ResultMC']['DaEIE'] + output_SXD['ResultMC']['DaEIR'])
Ptarget_SXD = (np.polyval(p_partflux['SXD'], den_SXD)[:,None]*[0.8,1,1.2]/ev)*(epsilon + gamma*np.transpose(np.quantile(Te_avg_SXD[indx_SP_SXD,:,:],[0.16,0.5,0.84],1)))*ev
Te_avg_ED = (output_ED['ResultMC']['TeEMC']*output_ED['ResultMC']['DaEIE'] + output_ED['ResultMC']['TeRMC']*output_ED['ResultMC']['DaEIR'] + output_ED['ResultMC']['DaMOL']*(output_ED['ResultMC']['TeRMC'] + R*(output_ED['ResultMC']['TeEMC']-output_ED['ResultMC']['TeRMC'])))/(output_ED['ResultMC']['DaMOL'] + output_ED['ResultMC']['DaEIE'] + output_ED['ResultMC']['DaEIR'])
Ptarget_ED = (np.polyval(p_partflux['ECD'], den_ED)[:,None]*[0.8,1,1.2]/ev)*(epsilon + gamma*np.transpose(np.quantile(Te_avg_ED[indx_SP_ED,:,:],[0.16,0.5,0.84],1)))*ev
Te_avg_CD = (output_CD['ResultMC']['TeEMC']*output_CD['ResultMC']['DaEIE'] + output_CD['ResultMC']['TeRMC']*output_CD['ResultMC']['DaEIR'] + output_CD['ResultMC']['DaMOL']*(output_CD['ResultMC']['TeRMC'] + R*(output_CD['ResultMC']['TeEMC']-output_CD['ResultMC']['TeRMC'])))/(output_CD['ResultMC']['DaMOL'] + output_CD['ResultMC']['DaEIE'] + output_CD['ResultMC']['DaEIR'])
Ptarget_CD = (np.polyval(p_partflux['CD'], den_CD)[:,None]*[0.8,1,1.2]/ev)*(epsilon + gamma*np.transpose(np.quantile(Te_avg_CD[indx_SP_CD,:,:],[0.16,0.5,0.84],1)))*ev

#peak heat load
Pptarget_SXD =(np.polyval(p_peakpartflux['SXD'], den_SXD)[:,None]*[0.8,1,1.2]/ev)*(epsilon + gamma*np.transpose(np.quantile(Te_avg_SXD[indx_SP_SXD,:,:],[0.16,0.5,0.84],1)))*ev
Pptarget_ED =(np.polyval(p_peakpartflux['ECD'], den_ED)[:,None]*[0.8,1,1.2]/ev)*(epsilon + gamma*np.transpose(np.quantile(Te_avg_ED[indx_SP_ED,:,:],[0.16,0.5,0.84],1)))*ev
Pptarget_CD =(np.polyval(p_peakpartflux['CD'], den_CD)[:,None]*[0.8,1,1.2]/ev)*(epsilon + gamma*np.transpose(np.quantile(Te_avg_CD[indx_SP_CD,:,:],[0.16,0.5,0.84],1)))*ev

import matplotlib.pyplot as plt
import dms.general_tools as gt

#codes that help the plotting - unc_plot, inpaint_nans, smooth
def inpaint_nans(Ain):
    #original code by John D'Errico for Matlab, this is a conversion for case=4 (spring metaphor)

    n, m = np.shape(Ain)
    import copy
    A = copy.deepcopy(Ain)
    A = np.ndarray.flatten(A,order='F')
    nm = n*m
    k = np.isnan(A)

    nan_list = np.argwhere(k)
    known_list = np.argwhere(np.logical_not(k))
    nan_count = len(nan_list)

    # convert NaN indices to (r,c) form
    nr, nc = np.unravel_index(nan_list,[n,m],'F')

    # both forms of index in one array:
    # column 1 == unrolled index
    # column 2 == row index
    # column 3 == column index
    nan_list = np.concatenate((nan_list,nr,nc),axis=1)

    #spring analogy - interpolating operator

    #list of all springs between a node and a horizontal/vertical neighbour

    hv_list = np.array([[-1,-1,0],[1,1,0],[-n,0,-1],[n,0,1]])
    hv_springs = np.empty([0,2])
    import numpy.matlib
    for i in range(0,4):
        hvs = nan_list + np.matlib.repmat(hv_list[i,:],nan_count,1)
        kk = np.logical_and(np.logical_and(hvs[:,1]>=0, hvs[:,1]<=(n-1)), np.logical_and(hvs[:,2]>=0, hvs[:,2]<=(m-1)))
        hv_springs = np.concatenate((hv_springs, np.concatenate((np.expand_dims(nan_list[kk,0],axis=1), np.expand_dims(hvs[kk,0],axis=1)),axis=1)),axis=0)

    hv_springs = np.unique(np.sort(hv_springs,axis=1),axis=0)

    # build sparse matrix of connections, springs
    # connecting diagonal neighbors are weaker than
    # the horizontal and vertical springs

    nhv = len(hv_springs[:,0])

    #matlab sparse(i,j,v,m,n) -> python scipy scipy.sparse.csr_matrix((v,(i,j)), shape=[nhv,nm])

    i = np.matlib.repmat(np.expand_dims(np.arange(0,nhv,1),axis=1),1,2)
    j = hv_springs
    v = np.matlib.repmat(np.expand_dims([1,-1],axis=1).T,nhv,1)

    #make sparse matrix
    springs = csr_matrix((np.ndarray.flatten(v,'F'),(np.ndarray.flatten(i,'F').astype(int),np.ndarray.flatten(j,'F').astype(int))),shape=(nhv,nm))

    rhs=-springs[:,np.ndarray.flatten(known_list,'F')]*A[np.ndarray.flatten(known_list,'F')]

    B = copy.deepcopy(A)

    solv = lsqr(springs[:,nan_list[:,0]],rhs)

    B[nan_list[:,0]]=solv[0]

    B=np.reshape(B,(n,m),order='F')

    return B

def smooth(x,window=5,p=2):
    return savgol_filter(x,window,p)

def unc_plot(x,data,color='r',apply_smooth=True):
    #plot with uncertainty bounds
    import matplotlib.pyplot as plt
    if apply_smooth:
        datan = inpaint_nans(data)
        y = smooth(np.ravel(datan[:,1]))
        yL= smooth(np.ravel(datan[:,0]))
        yU= smooth(np.ravel(datan[:,2]))
    else:
        y = np.ravel(data[:,1])
        yL= np.ravel(data[:,0])
        yU= np.ravel(data[:,2])
    plt.plot(np.ravel(x),y,color=color,linewidth=2)
    plt.fill_between(np.ravel(x),yL,yU,alpha=0.2,antialiased=True,color=color)

#Calculate figures

if Figure4:
    print('Making particle balance figures (integrated) - Figure 4a,b,c')

    plt.figure()
    unc_plot(den_SXD, (partflux_SXD/1.602e-19)[:,None]*[0.8,1,1.2],color='k')
    unc_plot(den_SXD,output_SXD['Result']['ion_I'],color='r')
    unc_plot(den_SXD, (partflux_SXD/1.602e-19)[:,None]*[0.8,1,1.2] + output_SXD['Result']['ion_sink_tot_I'],color='m')
    unc_plot(den_SXD, output_SXD['Result']['EIR_I'],color='b')
    unc_plot(den_SXD, output_SXD['Result']['mar_h2p_I'],color='g')
    plt.xlabel('Greenwald fraction')
    plt.ylabel('Particles / s')
    plt.title('Super-X')
    plt.axis([0.25, 0.5, 0, 4e22])
    plt.savefig('SXD_partflux.eps')

    plt.figure()
    unc_plot(den_ED, (partflux_ED/1.602e-19)[:,None]*[0.8,1,1.2],color='k')
    unc_plot(den_ED,output_ED['Result']['ion_I'],color='r')
    unc_plot(den_ED, (partflux_ED/1.602e-19)[:,None]*[0.8,1,1.2] + output_ED['Result']['ion_sink_tot_I'],color='m')
    unc_plot(den_ED, output_ED['Result']['EIR_I'],color='b')
    unc_plot(den_ED, output_ED['Result']['mar_h2p_I'],color='g')
    plt.xlabel('Greenwald fraction')
    plt.ylabel('Particles / s')
    plt.title('Elongated Divertor')
    plt.axis([0.25, 0.5, 0, 4e22])
    plt.savefig('ED_partflux.eps')

    plt.figure()
    unc_plot(den_CD, (partflux_CD/1.602e-19)[:,None]*[0.8,1,1.2],color='k')
    unc_plot(den_CD,output_CD['Result']['ion_I'],color='r')
    unc_plot(den_CD, (partflux_CD/1.602e-19)[:,None]*[0.8,1,1.2] + output_CD['Result']['ion_sink_tot_I'],color='m')
    unc_plot(den_CD, output_CD['Result']['EIR_I'],color='b')
    unc_plot(den_CD, output_CD['Result']['mar_h2p_I'],color='g')
    plt.xlabel('Greenwald fraction')
    plt.ylabel('Particles / s')
    plt.title('Conventional Divertor')
    plt.axis([0.25, 0.5, 0, 4e22])
    plt.savefig('CD_partflux.eps')

    indx_SP_SXD = 0
    indx_SP_CD = 22
    indx_SP_ED = 10

    print('Generate ion source/sink profiles @ 35% Greenwald fraction - Figure 4d,e,f')
    # Profiles (ion sources/sinks)
    print('Calculate ion source / sink profiles')
    for nGW_req in [0.3,0.31,0.32,0.33,0.34,0.35]:
        # SXD
        # adjustments for spatial misalignment, based on aligning CD brightness profile SP scan (46895)
        #I = np.append(np.arange(0, 10, 0.5), np.arange(10, 20) - 3)
        #Ireq = np.linspace(0, 16, 30)
        #I = Ireq
        indx = np.nanargmin((den_SXD - nGW_req)**2)
        indx_G = np.nanargmin((G_SXD['time'] - output_SXD['input']['Time'][indx])**2)
        plt.figure()
        p_spat = np.polyfit(np.arange(0,30),G_SXD['L_m'][indx_G,:30],2)
        ioni = output_SXD['ResultMC']['Ion'] + output_SXD['ResultMC']['mai_h2'] + output_SXD['ResultMC']['mai_h2p']
        # correct for misalignments
        unc_plot(np.polyval(p_spat,np.arange(0,30)), np.transpose(np.quantile(ioni[:, indx, :], [0.16, 0.5, 0.84], axis=1)), color='r')
        mar = output_SXD['ResultMC']['mar_hm'] + output_SXD['ResultMC']['mar_h2p']
        unc_plot(np.polyval(p_spat,np.arange(0,30)), np.transpose(np.quantile(mar[:, indx, :], [0.16, 0.5, 0.84], axis=1)), color='g')
        unc_plot(np.polyval(p_spat,np.arange(0,30)), np.transpose(np.quantile(output_SXD['ResultMC']['EIR'][:, indx, :], [0.16, 0.5, 0.84], axis=1)),color='b')
        plt.axis([2, 2.8, 0, 8e21])
        plt.xlabel('Poloidal distance wrt midplane (m)')
        plt.ylabel('Particles / m^2 / s')
        plt.title('SXD @ ' + str(nGW_req) + ' n_{GW}')
        plt.savefig('SXD_nGW_'+str(nGW_req)+'.eps')
        # CD
        indx = np.nanargmin((den_CD - nGW_req)**2)
        indx_G = np.nanargmin((G_CD['time'] - output_CD['input']['Time'][indx])**2)
        plt.figure()
        p_spat = np.polyfit(np.arange(0,30),G_SXD['L_m'][indx_G,:30],2)
        ioni = output_CD['ResultMC']['Ion'] + output_CD['ResultMC']['mai_h2'] + output_CD['ResultMC']['mai_h2p']
        unc_plot(np.polyval(p_spat,np.arange(indx_SP_CD,30)), np.transpose(np.quantile(ioni[indx_SP_CD:, indx, :], [0.16, 0.5, 0.84], axis=1)), color='r')
        mar = output_CD['ResultMC']['mar_hm'] + output_CD['ResultMC']['mar_h2p']
        unc_plot(np.polyval(p_spat,np.arange(indx_SP_CD,30)), np.transpose(np.quantile(mar[indx_SP_CD:, indx, :], [0.16, 0.5, 0.84], axis=1)), color='g')
        unc_plot(np.polyval(p_spat,np.arange(indx_SP_CD,30)), np.transpose(np.quantile(output_CD['ResultMC']['EIR'][indx_SP_CD:, indx, :], [0.16, 0.5, 0.84], axis=1)),color='b')
        plt.axis([2, 2.8, 0, 8e21])
        plt.xlabel('Poloidal distance wrt midplane (m)')
        plt.ylabel('Particles / m^2 / s')
        plt.title('CD @ ' + str(nGW_req) + ' n_{GW}')
        plt.savefig('CD_nGW_'+str(nGW_req)+'.eps')
        # ED
        indx = np.nanargmin((den_ED - nGW_req)**2)
        indx_G = np.nanargmin((G_ED['time'] - output_ED['input']['Time'][indx])**2)
        plt.figure()
        p_spat = np.polyfit(np.arange(0,30),G_SXD['L_m'][indx_G,:30],2)
        ioni = output_ED['ResultMC']['Ion'] + output_ED['ResultMC']['mai_h2'] + output_ED['ResultMC']['mai_h2p']
        unc_plot(np.polyval(p_spat,np.arange(indx_SP_ED,30)), np.transpose(np.quantile(ioni[indx_SP_ED:, indx, :], [0.16, 0.5, 0.84], axis=1)), color='r')
        mar = output_ED['ResultMC']['mar_hm'] + output_ED['ResultMC']['mar_h2p']
        unc_plot(np.polyval(p_spat,np.arange(indx_SP_ED,30)), np.transpose(np.quantile(mar[indx_SP_ED:, indx, :], [0.16, 0.5, 0.84], axis=1)), color='g')
        unc_plot(np.polyval(p_spat,np.arange(indx_SP_ED,30)), np.transpose(np.quantile(output_ED['ResultMC']['EIR'][indx_SP_ED:, indx, :], [0.16, 0.5, 0.84], axis=1)),color='b')
        plt.axis([2, 2.8, 0, 8e21])
        plt.xlabel('Poloidal distance wrt midplane (m)')
        plt.ylabel('Particles / m^2 / s')
        plt.title('ED @ ' + str(nGW_req) + ' n_{GW}')
        plt.savefig('ED_nGW_'+str(nGW_req)+'.eps')

if Figure6:
    print('Generate detachment front evolution (ionisation front and EIR front) - figure 6a,b')

    plt.figure() #detachment front
    indx_SP_SXD = 0
    indx_SP_CD = 26
    indx_SP_ED = 12
    start_ED = 10
    start_CD = 22
    DetachFront_L = np.zeros(np.shape(den_SXD)[0]) + np.nan
    DetachFront_LL = np.zeros(np.shape(den_SXD)[0]) + np.nan
    DetachFront_LH = np.zeros(np.shape(den_SXD)[0]) + np.nan
    ioni = output_SXD['ResultMC']['Ion'] + output_SXD['ResultMC']['mai_h2'] + output_SXD['ResultMC']['mai_h2p']
    ioni = gaussian_filter(inpaint_nans(np.nanmedian(ioni,axis=2)),1)
    for i in range(0,len(DetachFront_L)):
        if not np.isnan(den_SXD[i]):
            indx = np.nanargmin((den_SXD - den_SXD[i])**2)
            indx_G = np.nanargmin((G_SXD['time'] - output_SXD['input']['Time'][indx])**2)
            p_spat = np.polyfit(np.arange(0,30),G_SXD['L'][indx_G,:30] - G_SXD['L'][indx_G,indx_SP_SXD],2)
            # correct for misalignments
            f = interp1d(ioni[indx_SP_SXD:,indx],np.polyval(p_spat,np.arange(0,30))[indx_SP_SXD:],'linear',bounds_error=False,fill_value=np.nan)
            DetachFront_L[i] = f(1.5e21)
            DetachFront_LL[i] = f(1.25e21)
            DetachFront_LH[i] = f(1.75e21)
    DetachFront = np.vstack((DetachFront_LL,DetachFront_L,DetachFront_LH)).transpose()
    unc_plot(den_SXD,DetachFront,color='b')
    DetachFront_L = np.zeros(np.shape(den_ED)[0]) + np.nan
    DetachFront_LL = np.zeros(np.shape(den_ED)[0]) + np.nan
    DetachFront_LH = np.zeros(np.shape(den_ED)[0]) + np.nan
    ioni = output_ED['ResultMC']['Ion'] + output_ED['ResultMC']['mai_h2'] + output_ED['ResultMC']['mai_h2p']
    ioni = gaussian_filter(inpaint_nans(np.nanmedian(ioni,axis=2)),1)
    for i in range(0,len(DetachFront_L)):
        if not np.isnan(den_ED[i]):
            indx = np.nanargmin((den_ED - den_ED[i])**2)
            indx_G = np.nanargmin((G_SXD['time'] - output_ED['input']['Time'][indx])**2)
            p_spat = np.polyfit(np.arange(0,30),G_SXD['L'][indx_G,:30]- G_SXD['L'][indx_G,indx_SP_ED],2)
            # correct for misalignments
            f = interp1d(ioni[start_ED:,indx],np.polyval(p_spat,np.arange(0,30))[start_ED:],'linear',bounds_error=False,fill_value=np.nan)
            DetachFront_L[i] = f(1.5e21)
            DetachFront_LL[i] = f(1.25e21)
            DetachFront_LH[i] = f(1.75e21)
    DetachFront = np.vstack((DetachFront_LL,DetachFront_L,DetachFront_LH)).transpose()
    unc_plot(den_ED,DetachFront,color='g')
    DetachFront_L = np.zeros(np.shape(den_CD)[0]) + np.nan
    DetachFront_LL = np.zeros(np.shape(den_CD)[0]) + np.nan
    DetachFront_LH = np.zeros(np.shape(den_CD)[0]) + np.nan
    ioni = output_CD['ResultMC']['Ion'] + output_CD['ResultMC']['mai_h2'] + output_CD['ResultMC']['mai_h2p']
    ioni = gaussian_filter(inpaint_nans(np.nanmedian(ioni,axis=2)),1)
    for i in range(0,len(DetachFront_L)):
        if not np.isnan(den_CD[i]) and i>30:
            indx = np.nanargmin((den_CD - den_CD[i])**2)
            indx_G = np.nanargmin((G_CD['time'] - output_CD['input']['Time'][indx])**2)
            p_spat = np.polyfit(np.arange(0,30),G_SXD['L'][indx_G,:30]- G_SXD['L'][indx_G,indx_SP_CD],1)
            # correct for misalignments
            f = interp1d(ioni[start_CD:, indx],np.polyval(p_spat,np.arange(0,30))[start_CD:],'linear',bounds_error=False,fill_value=(0, np.max(np.polyval(p_spat,np.arange(0,30)))))
            DetachFront_L[i] = f(1.5e21)
            DetachFront_LL[i] = f(1.25e21)
            DetachFront_LH[i] = f(1.7e21)
    DetachFront = np.sort(np.vstack((DetachFront_LL,DetachFront_L,DetachFront_LH)).transpose())
    unc_plot(den_CD,DetachFront,color='r')
    plt.axis([0.25,0.5,0, 0.7])
    plt.savefig('DetachFront_Ion.eps')

    #EIR front
    plt.figure() #detachment front
    indx_SP_SXD = 0
    indx_SP_CD = 26
    indx_SP_ED = 14
    start_ED = 10
    start_CD = 22
    DetachFront_L = np.zeros(np.shape(den_SXD)[0]) + np.nan
    DetachFront_LL = np.zeros(np.shape(den_SXD)[0]) + np.nan
    DetachFront_LH = np.zeros(np.shape(den_SXD)[0]) + np.nan
    EIR = output_SXD['ResultMC']['EIR']
    EIR = gaussian_filter(inpaint_nans(np.nanmedian(EIR,axis=2)),1)
    level = 3e20
    error = 0.5e20
    for i in range(0,len(DetachFront_L)):
        if not np.isnan(den_SXD[i]):
            indx = np.nanargmin((den_SXD - den_SXD[i])**2)
            indx_G = np.nanargmin((G_SXD['time'] - output_SXD['input']['Time'][indx])**2)
            p_spat = np.polyfit(np.arange(0,30),G_SXD['L'][indx_G,:30] - G_SXD['L'][indx_G,indx_SP_SXD],1)
            # correct for misalignments
            f = interp1d(EIR[indx_SP_SXD:,indx],np.polyval(p_spat,np.arange(0,30))[indx_SP_SXD:],'linear',bounds_error=False,fill_value=np.nan)
            DetachFront_L[i] = f(level)
            DetachFront_LL[i] = f(level-error)
            DetachFront_LH[i] = f(level+error)
    DetachFront = np.vstack((DetachFront_LL,DetachFront_L,DetachFront_LH)).transpose()
    print(den_SXD[np.argwhere(np.logical_not(np.isnan(np.nanmedian(DetachFront,axis=1))))[0]]*1.3)
    unc_plot(den_SXD,DetachFront,color='b',apply_smooth=False)
    DetachFront_L = np.zeros(np.shape(den_ED)[0]) + np.nan
    DetachFront_LL = np.zeros(np.shape(den_ED)[0]) + np.nan
    DetachFront_LH = np.zeros(np.shape(den_ED)[0]) + np.nan
    EIR = output_ED['ResultMC']['EIR']
    EIR = gaussian_filter(inpaint_nans(np.nanmedian(EIR,axis=2)),1)
    for i in range(0,len(DetachFront_L)):
        if not np.isnan(den_ED[i]):
            indx = np.nanargmin((den_ED - den_ED[i])**2)
            indx_G = np.nanargmin((G_SXD['time'] - output_ED['input']['Time'][indx])**2)
            p_spat = np.polyfit(np.arange(0,30),G_SXD['L'][indx_G,:30]- G_SXD['L'][indx_G,indx_SP_ED],1)
            # correct for misalignments
            f = interp1d(EIR[start_ED:,indx],np.polyval(p_spat,np.arange(0,30))[start_ED:],'linear',bounds_error=False,fill_value=np.nan)
            DetachFront_L[i] = f(level)
            DetachFront_LL[i] = f(level-error)
            DetachFront_LH[i] = f(level+error)
    DetachFront = np.vstack((DetachFront_LL,DetachFront_L,DetachFront_LH)).transpose()
    unc_plot(den_ED,DetachFront,color='g',apply_smooth=False)
    DetachFront_L = np.zeros(np.shape(den_CD)[0]) + np.nan
    DetachFront_LL = np.zeros(np.shape(den_CD)[0]) + np.nan
    DetachFront_LH = np.zeros(np.shape(den_CD)[0]) + np.nan
    EIR = output_CD['ResultMC']['EIR']
    EIR = gaussian_filter(inpaint_nans(np.nanmedian(EIR,axis=2)),1)
    for i in range(0,len(DetachFront_L)):
        if not np.isnan(den_CD[i]) and i>30:
            indx = np.nanargmin((den_CD - den_CD[i])**2)
            indx_G = np.nanargmin((G_CD['time'] - output_CD['input']['Time'][indx])**2)
            p_spat = np.polyfit(np.arange(0,30),G_SXD['L'][indx_G,:30]- G_SXD['L'][indx_G,indx_SP_CD],2)
            # correct for misalignments
            f = interp1d(EIR[start_CD:, indx],np.polyval(p_spat,np.arange(0,30))[start_CD:],'linear',bounds_error=False,fill_value=np.nan)
            DetachFront_L[i] = f(level)
            DetachFront_LL[i] = f(level-error)
            DetachFront_LH[i] = f(level+error)
    DetachFront = np.sort(np.vstack((DetachFront_LL,DetachFront_L,DetachFront_LH)).transpose())
    unc_plot(den_CD,DetachFront,color='r',apply_smooth=False)
    plt.axis([0.25,0.5,0, 0.7])
    plt.savefig('DetachFront_EIR.eps')

if Figure2:
    print('Generate peak heat fluxes (perpendicular) - figure 2c')

    plt.figure()
    unc_plot(den_SXD,Pptarget_SXD/1e6,color='b')
    unc_plot(den_ED,Pptarget_ED/1e6,color='g')
    unc_plot(den_CD,Pptarget_CD/1e6,color='r')
    plt.xlabel('Greenwald fraction')
    plt.ylabel('Peak target power load (MW/m2)')
    plt.yscale('log')
    plt.axis([0.25, 0.5, 1e-2, 1])
    plt.grid(True,ls='-')
    plt.savefig('Peak_target_power_load.eps')

    print('Show divertor neutral pressure - figure 2d')
    #show divertor neutral pressures
    pn = np.load(file_pn,allow_pickle=True)[()]

    den_SXD_pn = np.polyval(p_dens['46860']['p'], pn['46860']['t_l_div'])
    den_SXD_pn[np.logical_or(pn['46860']['t_l_div']<0.4, pn['46860']['t_l_div']>0.78)] = np.nan
    den_ED_pn = np.polyval(p_dens['47115']['p'], pn['47115']['t_l_div'])
    den_ED_pn[np.logical_or(pn['47115']['t_l_div']<0.3, pn['47115']['t_l_div']>0.8)] = np.nan
    den_CD_pn = np.polyval(p_dens['46866']['p'], pn['46866']['t_l_div'])
    den_CD_pn[np.logical_or(pn['46866']['t_l_div']<0.4, pn['46866']['t_l_div']>0.78)] = np.nan

    plt.figure()
    plt.plot(den_SXD_pn, pn['46860']['pn_l_div'],'b',linewidth=2)
    plt.plot(den_ED_pn, pn['47115']['pn_l_div'],'g',linewidth=2)
    plt.plot(den_CD_pn, pn['46866']['pn_l_div'],'r',linewidth=2)
    plt.axis([0.25,0.5,0,0.7])
    plt.xlabel('Greenwald fraction')
    plt.ylabel('Neutral pressure (Pa)')
    plt.savefig('Pn_div.eps')

if Figure5:

    print('Calculating power balance (Figure 5a,b,c)')
    ec = 13.6 * 1.602e-19
    # Total power loss
    plt.figure()
    import dms.general_tools as gt
    unc_plot(den_SXD, 0.002 * (
                output_SXD['Result']['Prad_tot_I'] + 2.2 * 1.602e-19 * output_SXD['Result']['tot_neutr_mol_I']),
                color='b')
    unc_plot(den_SXD, 0.002 * (
                output_SXD['Result']['Prad_tot_I'] + Ptarget_SXD + 2.2 * 1.602e-19 * output_SXD['Result'][
            'tot_neutr_mol_I']), color='k')
    unc_plot(den_SXD, 0.002 * (Ptarget_SXD), color='r')
    plt.xlabel('Greenwald fraction')
    plt.ylabel('Power loss (kW)')
    plt.title('Super-X Divertor')
    plt.axis([0.25, 0.5, 0, 500])
    plt.savefig('SXD_PowerLoss.eps')

    # Total power loss
    plt.figure()
    import dms.general_tools as gt

    unc_plot(den_ED,
                0.002 * (output_ED['Result']['Prad_tot_I'] + 2.2 * 1.602e-19 * output_ED['Result']['tot_neutr_mol_I']),
                color='b')
    unc_plot(den_ED, 0.002 * (output_ED['Result']['Prad_tot_I'] + Ptarget_ED + 2.2 * 1.602e-19 * output_ED['Result'][
        'tot_neutr_mol_I']), color='k')
    unc_plot(den_ED, 0.002 * (Ptarget_ED), color='r')
    plt.xlabel('Greenwald fraction')
    plt.ylabel('Power loss (kW)')
    plt.title('Elongated Divertor')
    plt.axis([0.25, 0.5, 0, 500])
    plt.savefig('ED_PowerLoss.eps')

    # Total power loss
    plt.figure()
    import dms.general_tools as gt

    unc_plot(den_CD,
                0.002 * (output_CD['Result']['Prad_tot_I'] + 2.2 * 1.602e-19 * output_CD['Result']['tot_neutr_mol_I']),
                color='b')
    unc_plot(den_CD, 0.002 * (output_CD['Result']['Prad_tot_I'] + Ptarget_CD + 2.2 * 1.602e-19 * output_CD['Result'][
        'tot_neutr_mol_I']), color='k')
    unc_plot(den_CD, 0.002 * (Ptarget_CD), color='r')
    plt.xlabel('Greenwald fraction')
    plt.ylabel('Power loss (kW)')
    plt.title('Conventional Divertor')
    plt.axis([0.25, 0.5, 0, 500])
    plt.savefig('CD_PowerLoss.eps')

    print('Calculating and showing power flow (Figure 5d)')


    def seg_arclength(x, y):
        # calculates segments along the curve given by x, y

        seglen = np.insert(np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2), 0, 0)

        return seglen


    def arclength(x, y):
        # calculates length along the curve given by x,y

        seglen = seg_arclength(x, y)
        cumulen = np.cumsum(seglen)

        return cumulen

    def geometric_integral(integrand, r, z, cumu=False):
        #code calculating the geometrical integral to obtain integrated parameters
        if len(np.shape(integrand)) == 1:
            integrand = np.expand_dims(integrand, axis=1)

        X = arclength(r, z)
        Y = 2 * np.pi * r[:, None] * integrand
        I = np.logical_not(np.isnan(X[:, None] * Y))

        if len(np.shape(integrand)) == 1:
            integrand = np.expand_dims(integrand, axis=1)

        if not cumu:
            R = np.zeros(np.shape(integrand)[1]) + np.nan

            for i in range(0, len(R)):
                R[i] = np.trapz(Y[I[:, i], i], X[I[:, i]])
        else:
            R = np.zeros([np.shape(integrand)[1], np.shape(integrand)[0]])

            for i in range(0, len(R)):
                R[i, 1:] = cumtrapz(Y[I[:, i], i], X[I[:, i]])

        return R

    def calc_power_flow(output_SXD, indx, partflux_SXD, SP_loc=0):
        import dms.general_tools as gt
        Prad = output_SXD['ResultMC']['rad_h2p'] + output_SXD['ResultMC']['rad_h2'] + output_SXD['ResultMC']['rad_hm'] + \
               output_SXD['ResultMC']['Pexc'] + output_SXD['ResultMC']['Prec']
        Dissoc = output_SXD['ResultMC']['mar_h2p'] + output_SXD['ResultMC']['mad_h2p'] + output_SXD['ResultMC'][
            'mai_h2p']
        Ploss_tot = Prad + 0 * 2 * Dissoc * 2.2 * 1.602e-19
        Ploss_C = geometric_integral(Ploss_tot[SP_loc:, indx, :], output_SXD['input']['R'][SP_loc:, indx],
                                        output_SXD['input']['Z'][SP_loc:, indx], cumu=True)
        Power_target = (partflux_SXD[indx] / 1.602e-19) * ((0.4 * np.random.rand(500)) + 0.8) * 1.602e-19 * (
                    epsilon + gamma * Te_avg_SXD[SP_loc, indx, :]) + geometric_integral(
            2 * Dissoc[SP_loc:, indx, :] * 2.2 * 1.602e-19, output_SXD['input']['R'][SP_loc:, indx],
            output_SXD['input']['Z'][SP_loc:, indx])
        Power_flow = Power_target[:, None] + Ploss_C
        return Power_flow


    def calc_Ploss_tot(output_SXD):
        import dms.general_tools as gt
        Prad = output_SXD['ResultMC']['rad_h2p'] + output_SXD['ResultMC']['rad_h2'] + output_SXD['ResultMC']['rad_hm'] + \
               output_SXD['ResultMC']['Pexc'] + output_SXD['ResultMC']['Prec']
        Dissoc = output_SXD['ResultMC']['mar_h2p'] + output_SXD['ResultMC']['mad_h2p'] + output_SXD['ResultMC'][
            'mai_h2p']
        Ploss_tot = Prad + 2 * Dissoc * 2.2 * 1.602e-19
        return Ploss_tot


    indx_SP_SXD = 0
    indx_SP_CD = 27
    indx_SP_ED = 14

    Dist_Midplane_Entrance = 1.5
    f_GW = 0.35
    plt.figure()
    indx = np.nanargmin(np.abs(f_GW - den_SXD))
    power_flow_SXD = calc_power_flow(output_SXD, indx, partflux_SXD)
    indx_G = np.nanargmin((G_SXD['time'] - output_SXD['input']['Time'][indx]) ** 2)
    p_spat = np.polyfit(np.arange(0, 30), G_SXD['L_m'][indx_G, :30], 2)
    unc_plot(np.polyval(p_spat, np.arange(indx_SP_SXD, 30)) - Dist_Midplane_Entrance,
                2 * np.quantile(power_flow_SXD, [0.16, 0.5, 0.84], axis=0).transpose(), color='b', apply_smooth=False)

    indx = np.nanargmin(np.abs(f_GW - den_ED))
    power_flow_ED = calc_power_flow(output_ED, indx, partflux_ED, SP_loc=indx_SP_ED)
    indx_G = np.nanargmin((G_SXD['time'] - output_SXD['input']['Time'][indx]) ** 2)
    p_spat = np.polyfit(np.arange(0, 30), G_SXD['L_m'][indx_G, :30], 2)
    unc_plot(np.polyval(p_spat, np.arange(indx_SP_ED, 30)) - Dist_Midplane_Entrance,
                2 * np.quantile(power_flow_ED, [0.16, 0.5, 0.84], axis=0).transpose(), color='g', apply_smooth=False)

    indx = np.nanargmin(np.abs(f_GW - den_CD))
    power_flow_CD = calc_power_flow(output_CD, indx, partflux_CD, SP_loc=indx_SP_CD)
    indx_G = np.nanargmin((G_SXD['time'] - output_SXD['input']['Time'][indx]) ** 2)
    p_spat = np.polyfit(np.arange(0, 30), G_SXD['L_m'][indx_G, :30], 2)
    unc_plot(np.polyval(p_spat, np.arange(indx_SP_CD, 30)) - Dist_Midplane_Entrance,
                2 * np.quantile(power_flow_CD, [0.16, 0.5, 0.84], axis=0).transpose(), color='r', apply_smooth=False)

    plt.axis([0.5, 1.3, 0, 4e5])
    plt.savefig('Powerflow_SXD_ED_CD.eps')

if Figure3:
    print('Plotting PSOL (figure 3a)')
    psol = np.load(file_psol, allow_pickle=True)[()]

    den_SXD_psol = np.polyval(p_dens['46769']['p'], psol['46769']['t'])
    den_SXD_psol[np.logical_or(psol['46769']['t'] < 0.4, psol['46769']['t'] > 0.78)] = np.nan
    den_ED_psol = np.polyval(p_dens['47079']['p'], psol['47079']['t'])
    den_ED_psol[np.logical_or(psol['47079']['t'] < 0.3, psol['47079']['t'] > 0.8)] = np.nan
    den_CD_psol = np.polyval(p_dens['46866']['p'], psol['46866']['t'])
    den_CD_psol[np.logical_or(psol['46866']['t'] < 0.4, psol['46866']['t'] > 0.78)] = np.nan

    plt.figure()
    SXD_disch = ['46707', '46769', '46791', '46792', '46792']
    plt.plot(den_SXD_psol, psol['46769']['psol'] / 1e6, label='Super-X Divertor', linewidth=2, color='b')
    plt.plot(den_ED_psol, psol['47079']['psol'] / 1e6, label='Elongated Divertor', linewidth=2, color='g')
    plt.plot(den_CD_psol, psol['46866']['psol'] / 1e6, label='Conventional Divertor', linewidth=2, color='r')
    plt.plot(den_SXD_psol, psol['46769']['pnbi'] / 1e6, label='Super-X Divertor', linestyle='--', linewidth=2,
             color='b')
    plt.plot(den_ED_psol, psol['47079']['pnbi'] / 1e6, label='Elongated Divertor', linestyle='--', linewidth=2,
             color='g')
    plt.plot(den_CD_psol, psol['46866']['pnbi'] / 1e6, label='Conventional Divertor', linestyle='--', linewidth=2,
             color='r')
    plt.plot(den_SXD_psol, psol['46769']['prad_core'] / 1e6, label='Super-X Divertor', linestyle=':', linewidth=2,
             color='b')
    plt.plot(den_ED_psol, psol['47079']['prad_core'] / 1e6, label='Elongated Divertor', linestyle=':', linewidth=2,
             color='g')
    plt.plot(den_CD_psol, psol['46866']['prad_core'] / 1e6, label='Conventional Divertor', linestyle=':', linewidth=2,
             color='r')
    plt.axis([0.25, 0.5, 0, 1.5])
    plt.xlabel('Greenwald fraction (%)')
    plt.ylabel('P_{SOL} (MW)')
    plt.savefig('PSOL.eps')

    print('Plotting Thomson information (figure 3b,c,d,e,f,g')
    # Investigte Thomson

    rshift = -0.005
    start_p = 10
    psi_interp_l = 0.95

    TS = loadmat(file_TS_SXD)
    psi_n = np.zeros(np.shape(TS['R']))
    for i in range(0, len(TS['ne'][:, 0])):
        pos = interp1d(TS['psi_n'][i, start_p:], np.arange(start_p, np.shape(TS['psi_n'])[1]))(1)
        rn = interp1d(np.arange(start_p, np.shape(TS['psi_n'])[1]), TS['R'][i, start_p:])(pos) + rshift
        posn = interp1d(TS['R'][i, start_p:], np.arange(start_p, np.shape(TS['psi_n'])[1]))(rn)
        psi_i = interp1d(np.arange(start_p, np.shape(TS['psi_n'])[1]), TS['psi_n'][i, start_p:])
        psi_n[i, :] = TS['psi_n'][i, :] - psi_i(posn) + psi_i(pos)

    TS['Te'][np.isinf(TS['Te'])] = np.nan
    TS_SXD = TS
    pe = TS['ne'] * TS['Te']

    Te_sep_SXD = np.zeros(len(TS['ne'][:, 0]))
    ne_sep_SXD = np.zeros(np.shape(Te_sep_SXD))
    pe_sep_SXD = np.zeros(np.shape(Te_sep_SXD))
    for i in range(0, len(TS['ne'][:, 0])):
        I = np.logical_and(psi_n[i, :] > psi_interp_l, np.logical_not(np.isnan(TS['Te'][i, :])))
        p = np.polyfit(psi_n[i, I], TS['Te'][i, I], 2)
        Te_sep_SXD[i] = np.polyval(p, 1)
        p = np.polyfit(psi_n[i, I], TS['ne'][i, I], 2)
        ne_sep_SXD[i] = np.polyval(p, 1)
        p = np.polyfit(psi_n[i, I], pe[i, I], 2)
        pe_sep_SXD[i] = np.polyval(p, 1)

    den_SXD_TS = np.polyval(p_dens['46707']['p'], TS['time'])
    den_SXD_TS[np.logical_or(TS['time'] < 0.4, TS['time'] > 0.78)] = np.nan

    TS = loadmat(file_TS_ED)
    psi_n = np.zeros(np.shape(TS['R']))
    for i in range(0, len(TS['ne'][:, 0])):
        pos = interp1d(TS['psi_n'][i, start_p:], np.arange(start_p, np.shape(TS['psi_n'])[1]))(1)
        rn = interp1d(np.arange(start_p, np.shape(TS['psi_n'])[1]), TS['R'][i, start_p:])(pos) + rshift
        posn = interp1d(TS['R'][i, start_p:], np.arange(start_p, np.shape(TS['psi_n'])[1]))(rn)
        psi_i = interp1d(np.arange(start_p, np.shape(TS['psi_n'])[1]), TS['psi_n'][i, start_p:])
        psi_n[i, :] = TS['psi_n'][i, :] - psi_i(posn) + psi_i(pos)

    TS['Te'][np.isinf(TS['Te'])] = np.nan
    TS_ED = TS
    pe = TS['ne'] * TS['Te']

    Te_sep_ED = np.zeros(len(TS['ne'][:, 0]))
    ne_sep_ED = np.zeros(np.shape(Te_sep_ED))
    pe_sep_ED = np.zeros(np.shape(Te_sep_ED))
    for i in range(0, len(TS['ne'][:, 0])):
        I = np.logical_and(psi_n[i, :] > psi_interp_l, np.logical_not(np.isnan(TS['Te'][i, :])))
        p = np.polyfit(psi_n[i, I], TS['Te'][i, I], 2)
        Te_sep_ED[i] = np.polyval(p, 1)
        p = np.polyfit(psi_n[i, I], TS['ne'][i, I], 2)
        ne_sep_ED[i] = np.polyval(p, 1)
        p = np.polyfit(psi_n[i, I], pe[i, I], 2)
        pe_sep_ED[i] = np.polyval(p, 1)

    den_ED_TS = np.polyval(p_dens['47079']['p'], TS['time'])
    den_ED_TS[np.logical_or(TS['time'] < 0.4, TS['time'] > 0.78)] = np.nan

    TS = loadmat(file_TS_CD)
    psi_n = np.zeros(np.shape(TS['R']))
    for i in range(0, len(TS['ne'][:, 0])):
        pos = interp1d(TS['psi_n'][i, start_p:], np.arange(start_p, np.shape(TS['psi_n'])[1]))(1)
        rn = interp1d(np.arange(start_p, np.shape(TS['psi_n'])[1]), TS['R'][i, start_p:])(pos) + rshift
        posn = interp1d(TS['R'][i, start_p:], np.arange(start_p, np.shape(TS['psi_n'])[1]))(rn)
        psi_i = interp1d(np.arange(start_p, np.shape(TS['psi_n'])[1]), TS['psi_n'][i, start_p:])
        psi_n[i, :] = TS['psi_n'][i, :] - psi_i(posn) + psi_i(pos)

    TS['Te'][np.isinf(TS['Te'])] = np.nan
    TS_CD = TS
    pe = TS['ne'] * TS['Te']

    Te_sep_CD = np.zeros(len(TS['ne'][:, 0]))
    ne_sep_CD = np.zeros(np.shape(Te_sep_CD))
    pe_sep_CD = np.zeros(np.shape(Te_sep_CD))
    for i in range(0, len(TS['ne'][:, 0])):
        I = np.logical_and(psi_n[i, :] > psi_interp_l, np.logical_not(np.isnan(TS['Te'][i, :])))
        p = np.polyfit(psi_n[i, I], TS['Te'][i, I], 2)
        Te_sep_CD[i] = np.polyval(p, 1)
        p = np.polyfit(psi_n[i, I], TS['ne'][i, I], 2)
        ne_sep_CD[i] = np.polyval(p, 1)
        p = np.polyfit(psi_n[i, I], pe[i, I], 2)
        pe_sep_CD[i] = np.polyval(p, 1)

    den_CD_TS = np.polyval(p_dens['46762']['p'], TS['time'])
    den_CD_TS[np.logical_or(TS['time'] < 0.3, TS['time'] > 0.78)] = np.nan

    p_nesep_SXD = np.polyfit(np.squeeze(den_SXD_TS)[np.logical_not(np.isnan(np.squeeze(den_SXD_TS) * ne_sep_SXD))],
                             ne_sep_SXD[np.logical_not(np.isnan(np.squeeze(den_SXD_TS) * ne_sep_SXD))], 3)
    p_nesep_CD = np.polyfit(np.squeeze(den_CD_TS)[np.logical_not(np.isnan(np.squeeze(den_CD_TS) * ne_sep_CD))],
                            ne_sep_CD[np.logical_not(np.isnan(np.squeeze(den_CD_TS) * ne_sep_CD))], 3)
    p_nesep_ED = np.polyfit(np.squeeze(den_ED_TS)[np.logical_not(np.isnan(np.squeeze(den_ED_TS) * ne_sep_ED))],
                            ne_sep_ED[np.logical_not(np.isnan(np.squeeze(den_ED_TS) * ne_sep_ED))], 3)

    plt.figure()
    plt.plot(np.squeeze(den_SXD_TS), ne_sep_SXD, 'bx')
    plt.plot(np.squeeze(den_ED_TS), ne_sep_ED, 'go')
    plt.plot(np.squeeze(den_CD_TS), ne_sep_CD, 'r+')
    plt.axis([0.25, 0.5, 0, 2e19])
    plt.savefig('ne_sep.eps')

    plt.figure()
    plt.plot(np.squeeze(den_SXD_TS), Te_sep_SXD, 'bx')
    plt.plot(np.squeeze(den_ED_TS), Te_sep_ED, 'go')
    plt.plot(np.squeeze(den_CD_TS), Te_sep_CD, 'r+')
    plt.axis([0.25, 0.5, 0, 70])
    plt.savefig('Te_sep.eps')


    print('...')

    for nGW_req in [0.35, 0.45]:
        plt.figure()
        indx = np.nanargmin((den_SXD_TS - nGW_req) ** 2)
        plt.plot(np.sqrt(TS_SXD['psi_n'][indx, :]), TS_SXD['Te'][indx, :], 'bx')
        indx = np.nanargmin((den_ED_TS - nGW_req) ** 2)
        plt.plot(np.sqrt(TS_ED['psi_n'][indx, :]), TS_ED['Te'][indx, :], 'go')
        indx = np.nanargmin((den_CD_TS - nGW_req) ** 2)
        plt.plot(np.sqrt(TS_CD['psi_n'][indx, :]), TS_CD['Te'][indx, :], 'r+')
        plt.axis([0, 1.1, 0, 1200])
        plt.xlabel('Rho normalised')
        plt.ylabel('Te (eV)')
        plt.title('Temperature @ ' + str(nGW_req) + ' n_{GW}')
        plt.savefig('Te_nGW_' + str(nGW_req) + '.eps')
        plt.figure()
        indx = np.nanargmin((den_SXD_TS - nGW_req) ** 2)
        plt.plot(np.sqrt(TS_SXD['psi_n'][indx, :]), TS_SXD['ne'][indx, :], 'bx')
        indx = np.nanargmin((den_ED_TS - nGW_req) ** 2)
        plt.plot(np.sqrt(TS_ED['psi_n'][indx, :]), TS_ED['ne'][indx, :], 'go')
        indx = np.nanargmin((den_CD_TS - nGW_req) ** 2)
        plt.plot(np.sqrt(TS_CD['psi_n'][indx, :]), TS_CD['ne'][indx, :], 'r+')
        plt.axis([0, 1.1, 0, 6e19])
        plt.xlabel('Rho normalised')
        plt.ylabel('ne (m^-3)')
        plt.title('Density @ ' + str(nGW_req) + ' n_{GW}')
        plt.savefig('ne_nGW_' + str(nGW_req) + '.eps')

if Ploss_PMI:

    print('Calculating power balance to obtain PMI power losses')
    ec = 13.6 * 1.602e-19
    # Total power loss
    plt.figure()
    PMI_col = 140*np.ones(np.shape(den_SXD))
    PMI_col_unc = PMI_col[:, None] * [0.8, 1, 1.2]
    import dms.general_tools as gt
    def gen_sample(A,N=1000):
        #generates samples from a Nx3 matrix with uncertainties
        rand_samp = np.random.rand(N)[None,:]
        return A[:,0][:,None]+(A[:,2]-A[:,1])[:,None]*rand_samp
    PMI_PEN = np.transpose(np.quantile((0.002 *gen_sample(output_SXD['Result']['P_MARMAD_net_I']) +gen_sample(PMI_col_unc)) / (0.002*gen_sample(Ptarget_SXD) + 0.002 * gen_sample(output_SXD['Result']['P_MARMAD_net_I'])+gen_sample(PMI_col_unc) + 0.002*gen_sample(output_SXD['Result']['Prec_net_I'])),[0.16,0.5,0.84],axis=1))
    PMI_PEN_col = np.transpose(np.quantile((0.002 *gen_sample(output_SXD['Result']['P_MARMAD_net_I'])) / (0.002*gen_sample(Ptarget_SXD) + 0.002 * gen_sample(output_SXD['Result']['P_MARMAD_net_I'])+gen_sample(PMI_col_unc) + 0.002*gen_sample(output_SXD['Result']['Prec_net_I'])),[0.16,0.5,0.84],axis=1))
    unc_plot(den_SXD, PMI_PEN,color='b')
    unc_plot(den_SXD,PMI_PEN_col,color='g')
    plt.xlabel('Greenwald fraction')
    plt.ylabel('Power loss fraction - PMI - downstream ion source')
    plt.title('Super-X Divertor')
    plt.axis([0.25,0.5,0,1])
    plt.savefig('PMI_ploss_fraction')


print('Done')
