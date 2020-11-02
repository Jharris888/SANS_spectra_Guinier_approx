# Guinier Approx. 
# 5.11.20 

import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.stats
from scipy import optimize

##################################################################################

# experimental spectrum 
PI_Q = np.load('Analysis/IQ_pambou.npy')
Q = np.load('Analysis/Q_pambou.npy')
simI_Q = np.load('Analysis/all_clusters_IQ.npy')
regI_Q = np.load('Analysis/region_44_66_IQ.npy')

PI_Q = PI_Q/PI_Q[0]
simI_Q = simI_Q/simI_Q[0]

Q2 = Q**2

lnPIQ = np.log(PI_Q)
lnIQ = np.log(simI_Q)
lnREG = np.log(regI_Q)

# Q*Rgy < 3**.5
# guess that Rgy approx. 20 
# Q**2 is less than 0.0075 

r1 = 20
r2 = np.where(np.round(Q2,3) == 0.008)[0][0]

mIp = np.round(scipy.stats.linregress(Q2[r1:r2], lnPIQ[r1:r2])[0],3) # Pambou slope
mIq = np.round(scipy.stats.linregress(Q2[r1:r2], lnIQ[r1:r2])[0],3)
mreg = np.round(scipy.stats.linregress(Q2[r1:r2], lnREG[r1:r2])[0],3)

pb = np.round(scipy.stats.linregress(Q2[r1:r2], lnPIQ[r1:r2])[1],3)
sb = np.round(scipy.stats.linregress(Q2[r1:r2], lnIQ[r1:r2])[1],3)
regb = np.round(scipy.stats.linregress(Q2[r1:r2], lnREG[r1:r2])[1],3)

mstdp = np.round(scipy.stats.linregress(Q2[r1:r2], lnPIQ[r1:r2])[-1],3) # Pambou stderr
mstdIq = np.round(scipy.stats.linregress(Q2[r1:r2], lnIQ[r1:r2])[-1],3) 
mstdreg = np.round(scipy.stats.linregress(Q2[r1:r2], lnREG[r1:r2])[-1],3)

s_line = []
p_line = []
r_line = []
for i in range(r1,r2):
    s_line.append(mIq*Q2[i]+sb)
    p_line.append(mIp*Q2[i]+pb)
    r_line.append(mreg*Q2[i]+regb)

pRgy = (mIp*(-3))**0.5
errp = (mstdp*(3))**0.5
sRgy = (mIq*(-3))**0.5
errs = (mstdIq*(3))**0.05
rRgy = (mreg*(-3))**0.5
errr = (mstdreg*(3))**0.05

plt.figure(figsize=[5,3])

plt.scatter(Q2[:r2+2], lnPIQ[:r2+2], s=10, color='k',marker='o',alpha=0.6) 
plt.scatter(Q2[:r2+2], lnREG[:r2+2], s=10, color='r',marker='s',alpha=0.6) 
plt.scatter(Q2[:r2+2], lnIQ[:r2+2], s=10, color='b',marker='D',alpha=0.6)
plt.plot(Q2[r1:r2], p_line, color='k', label=r'experiment, %s $\pm$ %s'%(np.round(pRgy,1),np.round(errp,1)),linewidth=3,alpha=0.6)
plt.plot(Q2[r1:r2], r_line, color='r', label=r'fitted simulation, %s $\pm$ %s'%(np.round(rRgy,2),np.round(errr,2)),linewidth=3,alpha=0.6) 
plt.plot(Q2[r1:r2], s_line, color='b', label=r'all regions simulation, %s $\pm$ %s'%(np.round(sRgy,2),np.round(errs,2)),linewidth=3,alpha=0.6)

plt.legend(loc='lower left',frameon=False, fontsize=9)             
plt.xlabel('$Q^2$ ($\AA^{-2}$)',fontsize=15)
plt.xticks([0.002,0.004,0.006],fontsize=12) 
plt.yticks(fontsize=12)
plt.ylabel('ln[I($Q$)/I($Q_{0}$)]',fontsize=15) 
plt.tight_layout()
plt.minorticks_on()
plt.savefig('figures/guinier_rgy.pdf')
plt.show()  


