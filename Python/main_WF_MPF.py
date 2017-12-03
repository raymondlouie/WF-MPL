

########################################################################3
# Import packages

from numpy import *
from WF_sim_traj import WF_sim_traj
from estimate_MPL import estimate_MPL
import numpy as np
import scipy.linalg
import time
import matplotlib.pyplot as plt

plt.close("all")
########################################################################3
# set initial parameters

N=1000 # population size
L=int(7) # number of loci
dt_array=arange(1,10000,1) # generations sampled
mu = 1/N # mutation probability

# selection coefficients (0, -, +)
s =  np.vstack((np.zeros((int(floor(L/3)),1 )), (-10/N)*np.ones((int(floor(L/3)),1 )), (10/N)* np.ones((L-2*int(floor(L/3)), 1))))

no_runs =10 # number of WF trajectories to generate, for testing purposes

########################################################################3
# Generate WF runs and estimate the selection coefficients

K = 2**L

p_init  = np.ones(K,)
p_init = p_init/sum(p_init)

ind_pos = [index for index, elem in enumerate(s) if elem>0]
ind_neg = [index for index, elem in enumerate(s) if elem<0]
ind_neu = [index for index, elem in enumerate(s) if elem==0]

s_MPL_array  = zeros((no_runs,L))

for ind_run in range(no_runs):
    tstart =  time.clock()
    
      ########################################################################
  # Generate trajectories and estimate selection coefficients


    single_mut, double_mut = WF_sim_traj(s,mu,L,N,p_init,dt_array)

    s_MPL = estimate_MPL(mu,dt_array,single_mut,double_mut)
    s_MPL_array[ind_run,:] = s_MPL.T
    
    tend = time.clock()
    time_run = tend - tstart   
    print("Run number = ", ind_run+1, "/",no_runs,", run time=", time_run, "sec")
    
    
s_MPL_pos = s_MPL_array[:,ind_pos].flatten()
s_MPL_neg = s_MPL_array[:,ind_neg].flatten()
s_MPL_neu = s_MPL_array[:,ind_neu].flatten()



plt.hist(s_MPL_pos, bins=30,  color='b', alpha=0.5,  label='Positive')
plt.hist(s_MPL_neg, bins=30,   color='r', alpha=0.5, label='Negative')
plt.hist(s_MPL_neu, bins=30,  color='k', alpha=0.5, label='Neutral')
for svalue in s:
    plt.axvline(x=svalue)
plt.xlabel("Selection coefficient estimates")
plt.ylabel("Count of selection coefficient estimates")
plt.legend()
plt.show()    