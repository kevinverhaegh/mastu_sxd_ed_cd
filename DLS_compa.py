#Calculates the presented table on the DLS comparison to the experiments. Obtained from the exhaust analysis routines hosted on the UKAEA git repository (git.ccfe.ac.uk/jrh/mastu_exhaust_analysis), DOI: 10.1088/1741-4326/acea33

import matplotlib.pyplot as plt
import mastu_exhaust_analysis as mea
import numpy as np
from scipy.interpolate import interp1d

t = np.linspace(0.2,0.8,25)
dls = np.zeros(np.shape(t))
for i in range(0,len(t)):
    X=mea.divertor_geometry(shot=46895,time=t[i])
    X.get_dls_param()
    dls[i] = X.dls_shape_param

rs = X.efit_data['outer_strikepoint_lower_r']
trs= X.efit_data['t']
rsi = interp1d(trs,rs)

X = mea.divertor_geometry(shot=46860,time=0.45)
X.get_dls_param()
dls_SXD = X.dls_shape_param
X = mea.divertor_geometry(shot=47079,time=0.45)
X.get_dls_param()
dls_ED = X.dls_shape_param
X = mea.divertor_geometry(shot=46702,time=0.45)
X.get_dls_param()
dls_CD = X.dls_shape_param

plt.figure()
plt.plot(rsi(t),dls/dls_SXD, 'k',label='CD->ED->SXD scan')
plt.hlines(dls_CD/dls_SXD,0.6,1.4,colors='r',linestyles='dashed')
plt.hlines(dls_ED/dls_SXD,0.6,1.4,colors='g',linestyles='dashed')
plt.hlines(dls_SXD/dls_SXD,0.6,1.4,colors='b',linestyles='dashed')

plt.xlabel('Target radius (m)')
plt.ylabel('DLS onset norm. SXD')
plt.savefig('DLS_norm_MASTU_SXD')

plt.show()

np.savez('DLS_data.npz',dls_CD=dls_CD,dls_ED=dls_ED,dls_SXD=dls_SXD,dls_SXD_TCV=dls_SXD_TCV,dls_SP_scan=dls,r_dls_SP_scan=rsi(t),dls_TCV_L=TCV_L,r_dls_TCV_L=TCV_Rt_L,dls_TCV_Rt=TCV_Rt,r_dls_TCV_Rt=TCV_Rt_Rt, dls_TCV_Fx=TCV_Fx,r_dls_TCV_Fx=TCV_Rt_Fx,dls_TCV_XPT=TCV_XPT,r_dls_TCV_XPT=TCV_Rt_XPT,dls_TCV_SF=TCV_SF,r_dls_TCV_SF=TCV_Rt_SF)

#experimental calculations
X = np.load('DLS_data.npz',allow_pickle=True)

def Cprop(fGW):
    return ((fGW**1.46)*(((0.6/25)*fGW + 0.54)**-0.6))

def Cprop(fGW):
    PSOL = (((1.4-1)/(47.5-25))*fGW + 0.54)
    lambda_q = PSOL**0.14 * fGW**0.30
    if PSOL>1.5:
        PSOL = 1.5
    if PSOL<1:
        PSOL = 1
    qparallel = PSOL/lambda_q
    return fGW/(qparallel**(5/7))

from scipy.optimize import fsolve
fGW_CD_detach = 39
fGW_SXD_detach = fsolve(lambda tau : Cprop(tau)-(Cprop(fGW_CD_detach)/(X['dls_CD']/X['dls_SXD'])),40)
fGW_ED_detach = fsolve(lambda tau : Cprop(tau)-(Cprop(fGW_CD_detach)/(X['dls_CD']/X['dls_ED'])),40)

fGW_EIR_SXD = 33
fGW_EIR_ED_detach = fsolve(lambda tau : Cprop(tau)-(Cprop(fGW_EIR_SXD)/(X['dls_SXD']/X['dls_ED'])),40)
fGW_EIR_CD_detach = fsolve(lambda tau : Cprop(tau)-(Cprop(fGW_EIR_SXD)/(X['dls_SXD']/X['dls_CD'])),40)
