# June 3rd 2018

'''
find the number of GW events and the coincident evenst
'''


# ================
# load data 
#=================

# coding: UTF-8
import matplotlib.pyplot as plt
import numpy as np
import sys
from math import pi, cos
from math import log10, floor, factorial
import pandas as pd
from collections import OrderedDict
from matplotlib.ticker import MultipleLocator
from pylab import *


#file_id = sys.argv[1]
sample_number = 35. #float(sys.argv[2])
#file_id_str = str(file_id)
#title = str(sys.argv[3])

# Data table
df = pd.read_table('/Users/kentaromogushi/Documents/department/research/shortgamma/chi/Obtain_chi/Allobtain/swiftdata/plusjet/typi/onaxis/Karelle_revised_brokenOUTPUTHopkins1Gyrtheta_obsGW29theta_c1sigma.txt', sep=' ', names=('label','No','N88Mpc_No','N140Mpc_No', 'N120Mpc_No', 'N190Mpc_No', 'N300Mpc_No', 'N65Mpc_No', 'N105Mpc_No', 'N40Mpc_No', 'N125Mpc_No', 'N200Mpc_No', 'N225Mpc_No','N88MpcOnaxis_No', 'N140MpcOnaxis_No', 'N120MpcOnaxis_No', 'N190MpcOnaxis_No', 'N300MpcOnaxis_No', 'N65MpcOnaxis_No', 'N105MpcOnaxis_No', 'N40MpcOnaxis_No', 'N125MpcOnaxis_No', 'N200MpcOnaxis_No', 'N225MpcOnaxis_No', 's', 'theta_c', 'theta_obsGW'))


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   astromomical parameters and the correction factors     |
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c = 3.0*10**8   # [m/s]
H_o = 67.89/(3.0857*10**19) # [km/s/Mpc]*[Mpc/km]= [1/s]
Gpc = 3.0857*10**25 # [m/Gpc]
c_H_o = c/H_o/Gpc #[m]/[m/Gpc] = ~4.4[Gpc]
f_R = sample_number/107.
f_FOV = 0.1 # field of veiw
T = 12.6






#=====================================
# Split labels to make me easy to separate data depending on RF and theta_obsGW, etc.
#=====================================

# split label 
label_st = np.array(df['label'])
label_tmin = []
label_RF = []
label_LF = []
#label_theta_obsGW = [] # I choose theta_obsGW = 29 deg
label_theta_c = []
for ind in range(len(label_st)):
    # I need the labels of RF and theta_obs
    label_tmin.append(label_st[ind].split(',')[0])
    label_RF.append(label_st[ind].split(',')[1])
    label_LF.append(label_st[ind].split(',')[2])
    label_theta_c.append(label_st[ind].split(',')[5])


# add splitted labels
df['label_tmin'] = label_tmin
df['label_RF']=label_RF
df['label_LF'] = label_LF
df['label_theta_c']=label_theta_c


# LRD is No / (f_FOV*f_R*0.8*T*Vo) where Vo = c_H_o**(3)  [Gpc^-3 yr^-1]
df['LRD'] = df['No']/(f_FOV*f_R*0.8*T)*c_H_o**(-3) 



# re-gain N(z) which is No times N120Mpc_No 
# full angle cumulative number 
df['N88Mpc'] = df['No']*df['N88Mpc_No']   # LH combined in O2 for NSNS
df['N140Mpc'] = df['No']*df['N140Mpc_No'] # LH combined in O2 for NSBH, K at design NSNS
df['N120Mpc'] = df['No']*df['N120Mpc_No'] # LH in O3 for NSNS
df['N190Mpc'] = df['No']*df['N190Mpc_No'] # LH in O3 for NSBH and LH in O4+ for NSNS
df['N300Mpc'] = df['No']*df['N300Mpc_No'] # LH in O4+ for NSBH 
df['N65Mpc'] = df['No']*df['N65Mpc_No'] # V in O3 for NSNS and V in O4+ for NSNS
df['N105Mpc'] = df['No']*df['N105Mpc_No'] # V in O3 and O4+ for NSBH
df['N40Mpc'] = df['No']*df['N40Mpc_No'] # K in O4+ for NSNS
df['N125Mpc'] = df['No']*df['N125Mpc_No'] # V in design for NSNS
df['N200Mpc'] = df['No']*df['N200Mpc_No'] # V in design for NSBH
df['N225Mpc'] = df['No']*df['N225Mpc_No'] # K in design for NSBH



# only on-axis emisson 
df['N88MpcOnaxis'] = df['No']*df['N88MpcOnaxis_No'] # LH combined in O2 for NSNS
df['N140MpcOnaxis'] = df['No']*df['N140MpcOnaxis_No'] # LH combined in O2 for NSBH, K at design NSNS
df['N120MpcOnaxis'] = df['No']*df['N120MpcOnaxis_No'] # LH in O3 for NSNS
df['N190MpcOnaxis'] = df['No']*df['N190MpcOnaxis_No'] # LH in O3 for NSBH and LH in O4+ for NSNS
df['N300MpcOnaxis'] = df['No']*df['N300MpcOnaxis_No'] # LH in O4+ for NSBH 
df['N65MpcOnaxis'] = df['No']*df['N65MpcOnaxis_No']  # V in O3 for NSNS and V in O4+ for NSNS
df['N105MpcOnaxis'] = df['No']*df['N105MpcOnaxis_No'] # V in O3 and O4+ for NSBH
df['N40MpcOnaxis'] = df['No']*df['N40MpcOnaxis_No'] # K in O4+ for NSNS
df['N125MpcOnaxis'] = df['No']*df['N125MpcOnaxis_No'] # V in design for NSNS
df['N200MpcOnaxis'] = df['No']*df['N200MpcOnaxis_No'] # V in design for NSBH
df['N225MpcOnaxis'] = df['No']*df['N225MpcOnaxis_No'] # K in design for NSBH



# only off-axis emisson 
df['N88MpcOffaxis'] = df['N88Mpc'] - df['N88MpcOnaxis'] # LH combined in O2 for NSNS
df['N140MpcOffaxis'] = df['N140Mpc'] - df['N140MpcOnaxis'] # LH combined in O2 for NSBH, K at design NSNS
df['N120MpcOffaxis'] = df['N120Mpc'] -  df['N120MpcOnaxis'] # LH in O3 for NSNS
df['N190MpcOffaxis'] = df['N190Mpc'] - df['N190MpcOnaxis'] # LH in O3 for NSBH and LH in O4+ for NSNS
df['N300MpcOffaxis'] = df['N300Mpc'] - df['N300MpcOnaxis'] # LH in O4+ for NSBH 
df['N65MpcOffaxis'] = df['N65Mpc'] - df['N65MpcOnaxis']  # V in O3 for NSNS and V in O4+ for NSNS
df['N105MpcOffaxis'] = df['N105Mpc'] - df['N105MpcOnaxis'] # V in O3 and O4+ for NSBH
df['N40MpcOffaxis'] = df['N40Mpc'] - df['N40MpcOnaxis'] # K in O4+ for NSNS
df['N125MpcOffaxis'] = df['N125Mpc'] - df['N125MpcOnaxis']  # K in desing for NSNS
df['N200MpcOffaxis'] = df['N200Mpc'] - df['N200MpcOnaxis'] # V in design for NSBH
df['N225MpcOffaxis'] = df['N225Mpc'] - df['N225MpcOnaxis'] # K in design for NSBH



#==========================
# the number of GW events |
#==========================

import numpy as np
from itertools import combinations


# list of the NSNS inspiral ranges [Gpc]
# Hanford, Livingston, Virgo, (KAGRA)
R_NSNS_list_O3 = np.array([0.120, 0.120, 0.065]) # for O3 NSNS
R_NSNS_list_O4 = np.array([0.19, 0.19, 0.065]) # with LHV for O4 NSNS
R_NSNS_list_O4_K = np.array([0.19, 0.19, 0.065, 0.04]) # with LHVK for O4 NSNS
R_NSNS_list_design = np.array([0.19, 0.19 , 0.125, 0.14]) # for design NSNS
# the relative NSBH range to that of NSNS
r_NSBH_NSNS = 1.6
# list of the NSBH ranges
R_NSBH_list_O3 = r_NSBH_NSNS*R_NSNS_list_O3 # for O3 NSBH
R_NSBH_list_O4 = r_NSBH_NSNS*R_NSNS_list_O4 # with LHV for O4 NSBH
R_NSBH_list_O4_K = r_NSBH_NSNS*R_NSNS_list_O4_K # with LHVK for O4 NSBH
R_NSBH_list_design = r_NSBH_NSNS*R_NSNS_list_design # with LHVK for design NSBH

# duty cycle factor for each interferometer
D_single = 0.8


# ratio of the NSNS populaton to the combined population NSNS + NSBH
r_nsns = 5./6.

# ratio of the NSBH populaton to the combined population NSNS + NSBH
r_nsbh = 1./6.

# local rate desity depending on theta_c and a set of paramters of LF 
#label_LF = '1' # gives smallest
label_LF = '2' # gives biggest


# LRD at theta_c = 10.7 deg
LRD_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['LRD'].values[0]
# LRD at theta_c = 7.4 degs
LRD_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['LRD'].values[0]
# LRD at theta_c = 16.8
LRD_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['LRD'].values[0]
# LRD error for the upper bound 
LRD_up = abs(LRD_biggest -  LRD_mean)
# LRD error for the lower bound 
LRD_low = abs(LRD_smallest -  LRD_mean)




#==============
# GW events
#=============

def N_i_g(LRD_mean, LRD_up, LRD_low, Q, r_i, R_i_list):
	'''
	Calculate the numbe of GW events per year either for NSNS or NSBH
	rho_g = the local rate density of NSNS + NSBH 
	Q = the total number of interferometers.
	USAGE:
		N_NSNS_g_ave_O3, N_NSNS_g_err_O3 = N_i_g(LRD_ave, LRD_std, 3, r_nsns, R_NSNS_list_O3) # for NSNS in O3
		N_NSBH_g_ave_O3, N_NSBH_g_err_O3 = N_i_g(LRD_ave, LRD_std, 3, r_nsbh, R_NSBH_list_O3) # for NSBH in O3

		N_NSNS_g_ave_O4, N_NSNS_g_err_O4 = N_i_g(LRD_ave, LRD_std, 3, r_nsns, R_NSNS_list_O4) # for NSNS in O4
		N_NSBH_g_ave_O4, N_NSBH_g_err_O4 = N_i_g(LRD_ave, LRD_std, 3, r_nsbh, R_NSBH_list_O4) # for NSBH in O4

		N_NSNS_g_ave_O4_K, N_NSNS_g_err_O4_K = N_i_g(LRD_ave, LRD_std, 4, r_nsns, R_NSNS_list_O4_K) # for NSNS in O4 including K
		N_NSBH_g_ave_O4_K, N_NSBH_g_err_O4_K = N_i_g(LRD_ave, LRD_std, 4, r_nsbh, R_NSBH_list_O4_K) # for NSBH in O4 including K 		
		N_NSNS_g_ave_design, N_NSNS_g_err_design = N_i_g(LRD_ave, LRD_std, 4, r_nsns, R_NSNS_list_design) # for NSNS in design 
		N_NSBH_g_ave_design, N_NSBH_g_err_design = N_i_g(LRD_ave, LRD_std, 4, r_nsbh, R_NSBH_list_design) # for NSBH in design including K 
	'''
	VT_com = []
	# conincident for any two interferometr and coincident for any 3 interferomter, any i, ,,, any Q
	# the k starts with 2 ends the total number of interferometers
	for k in range(2, Q+1):
		# all the combinations of any k interferometers out of Q 
		comb = combinations(R_i_list, k)
		# convert to numpy array
		comb_inter = np.array(list(comb))
		#print('The combination of any %.f inteferometers out of %.f is' % (k, Q))
		#print(comb_inter)
		# looping with all the possible combination of k interferomters out of Q
		for j in range(len(comb_inter)):
			second_large_range = np.sort(comb_inter[j, :])[-2]
			#print('The second highest range among the combination is %f' % second_large_range)
			# the combined duty cycle and the second highest range among the combination
			VT_com.append(D_single**k*(1-D_single)**(Q-k)*4.*np.pi/3.*second_large_range**3)
	#print('the all the combination of VTs')
	#print(VT_com)
	result_mean = r_i*LRD_mean*sum(VT_com) # average value 
	result_up = r_i*LRD_up* sum(VT_com) # simply the error is linear to the error of LRD, error is due to 1-sigma of theta_c
	result_low = r_i*LRD_low* sum(VT_com) # simply the error is linear to the error of LRD, error is due to 1-sigma of theta_c
	return result_mean, result_up, result_low



#================
# GW evens     
#================

# the nsns and nsby number of events in O3
# the effective combined range is 88Mpc and the LH observatioin time is 117 days 
N_g_O2_mean = (117./365.)*(r_nsns*4*np.pi/3.*(0.088)**3 + r_nsbh*4*np.pi/3.*(0.14)**3)*LRD_mean
N_g_O2_up = (117./365.)*(r_nsns*4*np.pi/3.*(0.088)**3 + r_nsbh*4*np.pi/3.*(0.14)**3)*LRD_up
N_g_O2_low = (117./365.)*(r_nsns*4*np.pi/3.*(0.088)**3 + r_nsbh*4*np.pi/3.*(0.14)**3)*LRD_low


# the nsns and nsbh number of events in O3
# the number of NSNS events in O3, seen at least two detectors in three-detector network
N_NSNS_g_mean_O3, N_NSNS_g_up_O3, N_NSNS_g_low_O3 = N_i_g(LRD_mean, LRD_up,  LRD_low,3, r_nsns, R_NSNS_list_O3) # for NSNS in O3
# the number of NSBH events in O3
N_NSBH_g_mean_O3, N_NSBH_g_up_O3, N_NSBH_g_low_O3 = N_i_g(LRD_mean, LRD_up,  LRD_low, 3, r_nsbh, R_NSBH_list_O3) # for NSBH in O3
# the number of both NSNS and NSBH events
N_nsns_nsbh_g_mean_O3  = N_NSNS_g_mean_O3 + N_NSBH_g_mean_O3 
N_nsns_nsbh_g_up_O3 = (N_NSNS_g_up_O3**2 + N_NSBH_g_up_O3**2)**0.5 
N_nsns_nsbh_g_low_O3 = (N_NSNS_g_low_O3**2 + N_NSBH_g_low_O3**2)**0.5 

# the nsns and nsbh number of events in O4
# the number of NSNS events in O4, seen at least two detectors in three-detector network
N_NSNS_g_mean_O4, N_NSNS_g_up_O4, N_NSNS_g_low_O4 = N_i_g(LRD_mean, LRD_up,  LRD_low, 3, r_nsns, R_NSNS_list_O4) # for NSNS in O4
# the number of NSBH events in O4
N_NSBH_g_mean_O4, N_NSBH_g_up_O4, N_NSBH_g_low_O4 = N_i_g(LRD_mean, LRD_up,  LRD_low, 3, r_nsbh, R_NSBH_list_O4) # for NSBH in O4
# the number of both NSNS and NSBH events
N_nsns_nsbh_g_mean_O4 = N_NSNS_g_mean_O4  + N_NSBH_g_mean_O4 
N_nsns_nsbh_g_up_O4 = (N_NSNS_g_up_O4**2 + N_NSBH_g_up_O4**2)**0.5 
N_nsns_nsbh_g_low_O4 = (N_NSNS_g_low_O4**2 + N_NSBH_g_low_O4**2)**0.5 


# the nsns and nsbh number of events in O4 incluidng K
# the number of NSNS events in O4, seen at least two detectors in three-detector network
N_NSNS_g_mean_O4_K, N_NSNS_g_up_O4_K, N_NSNS_g_low_O4_K = N_i_g(LRD_mean, LRD_up,  LRD_low, 4, r_nsns, R_NSNS_list_O4_K) # for NSNS in O4 including K
# the number of NSBH events in O4
N_NSBH_g_mean_O4_K, N_NSBH_g_up_O4_K, N_NSBH_g_low_O4_K = N_i_g(LRD_mean, LRD_up,  LRD_low, 4, r_nsbh, R_NSBH_list_O4_K) # for NSBH in O4 including K 	
N_nsns_nsbh_g_mean_O4_K = N_NSNS_g_mean_O4_K  + N_NSBH_g_mean_O4_K 
N_nsns_nsbh_g_up_O4_K = (N_NSNS_g_up_O4_K**2 + N_NSBH_g_up_O4_K**2)**0.5 
N_nsns_nsbh_g_low_O4_K = (N_NSNS_g_low_O4_K**2 + N_NSBH_g_low_O4_K**2)**0.5 


# the nsns and nsbh number of events in design incluidng K
# the number of NSNS events in design, seen at least two detectors in three-detector network
N_NSNS_g_mean_design, N_NSNS_g_up_design, N_NSNS_g_low_design = N_i_g(LRD_mean, LRD_up,  LRD_low, 4, r_nsns, R_NSNS_list_design) # for NSNS in design including K
# the number of NSBH events in design
N_NSBH_g_mean_design, N_NSBH_g_up_design, N_NSBH_g_low_design = N_i_g(LRD_mean, LRD_up,  LRD_low, 4, r_nsbh, R_NSBH_list_design) # for NSBH in design including K 
N_nsns_nsbh_g_mean_design = N_NSNS_g_mean_design  + N_NSBH_g_mean_design 
N_nsns_nsbh_g_up_design = (N_NSNS_g_up_design**2 + N_NSBH_g_up_design**2)**0.5 
N_nsns_nsbh_g_low_design = (N_NSNS_g_low_design**2 + N_NSBH_g_low_design**2)**0.5 



#+++++++++++++++++++
# print out GW events 
#++++++++++++++++++++
print('the GW events for both NSNS and NSBH in O2 = %.1f - %.1f - %.1f' %(N_g_O2_mean - N_g_O2_low, N_g_O2_mean, N_g_O2_up + N_g_O2_mean))
print('check the NSNS events in O2')
N_nsns_O2_mean = (117./365.)*(r_nsns*4*np.pi/3.*(0.088)**3)*LRD_mean
N_nsns_O2_up = (117./365.)*(r_nsns*4*np.pi/3.*(0.088)**3)*LRD_up
N_nsns_O2_low = (117./365.)*(r_nsns*4*np.pi/3.*(0.088)**3)*LRD_low
print('NSNS events in O2 = %.1f - %.1f - %.1f' % (N_nsns_O2_mean - N_nsns_O2_low, N_nsns_O2_mean, N_nsns_O2_up + N_nsns_O2_mean))



print('='*50)
print('GW events in O3 = %.1f - %.1f - %.1f per year' % (N_nsns_nsbh_g_mean_O3 - N_nsns_nsbh_g_low_O3, N_nsns_nsbh_g_mean_O3, N_nsns_nsbh_g_mean_O3 +N_nsns_nsbh_g_up_O3, ))



print('='*50)
print('GW events in O4 = %.1f - %.1f - %.1f per year' % (N_nsns_nsbh_g_mean_O4 - N_nsns_nsbh_g_low_O4, N_nsns_nsbh_g_mean_O4, N_nsns_nsbh_g_mean_O4 +N_nsns_nsbh_g_up_O4))


print('='*50)
print('GW events in O4 with KAGRA = %.1f - %.1f - %.1f per year' % (N_nsns_nsbh_g_mean_O4_K - N_nsns_nsbh_g_low_O4_K,N_nsns_nsbh_g_mean_O4_K, N_nsns_nsbh_g_mean_O4_K + N_nsns_nsbh_g_up_O4_K))



print('='*50)
print('GW events in design with KAGRA = %.1f - %.1f - %.1f per year' % (N_nsns_nsbh_g_mean_design - N_nsns_nsbh_g_low_design, N_nsns_nsbh_g_mean_design, N_nsns_nsbh_g_mean_design + N_nsns_nsbh_g_up_design))
print('='*50)




#===================
# coincident events
#===================
# N88Mpc at theta_c = 16.8 degs
N88Mpc_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N88Mpc'].values[0]
# N88Mpc at theta_c = 10.6 degs
N88Mpc_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N88Mpc'].values[0]
# N88Mpc at theta_c =  7.4 degs
N88Mpc_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N88Mpc'].values[0]
# N88Mpc error for the upper bound 
N88Mpc_up = N88Mpc_biggest - N88Mpc_mean
# N88Mpc error for the lower bound 
N88Mpc_low = N88Mpc_mean - N88Mpc_smallest 

# N120Mpc at theta_c = 16.8 degs
N120Mpc_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N120Mpc'].values[0]
# N120Mpc at theta_c = 10.6 degs
N120Mpc_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N120Mpc'].values[0]
# N120Mpc at theta_c =  7.4 degs
N120Mpc_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N120Mpc'].values[0]
# N120Mpc error for the upper bound 
N120Mpc_up = N120Mpc_biggest - N120Mpc_mean
# N120Mpc error for the lower bound 
N120Mpc_low = N120Mpc_mean - N120Mpc_smallest


# N65Mpc at theta_c = 16.8 degs
N65Mpc_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N65Mpc'].values[0]
# N65Mpc at theta_c = 10.6 degs
N65Mpc_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N65Mpc'].values[0]
# N65Mpc at theta_c = 7.4 degs
N65Mpc_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N65Mpc'].values[0]
# N65Mpc error for the upper bound 
N65Mpc_up = N65Mpc_biggest - N65Mpc_mean
# N65Mpc error for the lower bound 
N65Mpc_low = N65Mpc_mean - N65Mpc_smallest

# N190Mpc at theta_c = 16.8 degs
N190Mpc_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N190Mpc'].values[0]
# N190Mpc at theta_c = 10.6 degs
N190Mpc_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N190Mpc'].values[0]
# N190Mpc at theta_c = 7.4 degs
N190Mpc_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N190Mpc'].values[0]
# N190Mpc error for the upper bound 
N190Mpc_up = N190Mpc_biggest - N190Mpc_mean
# N190Mpc error for the lower bound 
N190Mpc_low = N190Mpc_mean - N190Mpc_smallest

# N105Mpc at theta_c = 16.8 degs
N105Mpc_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N105Mpc'].values[0]
# N105Mpc at theta_c = 10.6 degs
N105Mpc_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N105Mpc'].values[0]
# N105Mpc at theta_c = 7.4 degs
N105Mpc_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N105Mpc'].values[0]
# N105Mpc error for the upper bound 
N105Mpc_up = N105Mpc_biggest - N105Mpc_mean
# N105Mpc error for the lower bound 
N105Mpc_low = N105Mpc_mean - N105Mpc_smallest

# N40Mpc at theta_c = 16.8 degs
N40Mpc_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N40Mpc'].values[0]
# N40Mpc at theta_c = 10.6 degs
N40Mpc_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N40Mpc'].values[0]
# N40Mpc at theta_c = 7.4 degs
N40Mpc_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N40Mpc'].values[0]
# N40Mpc error for the upper bound 
N40Mpc_up = N40Mpc_biggest - N40Mpc_mean
# N40Mpc error for the lower bound 
N40Mpc_low = N40Mpc_mean - N40Mpc_smallest

# N300Mpc at theta_c = 16.8 degs
N300Mpc_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N300Mpc'].values[0]
# N300Mpc at theta_c = 10.6 degs
N300Mpc_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N300Mpc'].values[0]
# N300Mpc at theta_c = 7.4 degs
N300Mpc_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N300Mpc'].values[0]
# N300Mpc error for the upper bound 
N300Mpc_up = N300Mpc_biggest - N300Mpc_mean
# N300Mpc error for the lower bound 
N300Mpc_low = N300Mpc_mean - N300Mpc_smallest


# N125Mpc at theta_c = 16.8 degs
N125Mpc_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N125Mpc'].values[0]
# N125Mpc at theta_c = 10.6 degs
N125Mpc_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N125Mpc'].values[0]
# N125Mpc at theta_c = 7.4 degs
N125Mpc_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N125Mpc'].values[0]
# N125Mpc error for the upper bound 
N125Mpc_up = N125Mpc_biggest - N125Mpc_mean
# N125Mpc error for the lower bound 
N125Mpc_low = N125Mpc_mean - N125Mpc_smallest

# N140Mpc at theta_c = 16.8 degs
N140Mpc_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N140Mpc'].values[0]
# N140Mpc at theta_c = 10.6 degs
N140Mpc_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N140Mpc'].values[0]
# N140Mpc at theta_c = 7.4 degs
N140Mpc_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N140Mpc'].values[0]
# N140Mpc error for the upper bound 
N140Mpc_up = N140Mpc_biggest - N140Mpc_mean
# N140Mpc error for the lower bound 
N140Mpc_low = N140Mpc_mean - N140Mpc_smallest

# N225Mpc at theta_c = 16.8 degs
N225Mpc_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N225Mpc'].values[0]
# N225Mpc at theta_c = 10.6 degs
N225Mpc_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N225Mpc'].values[0]
# N225Mpc at theta_c = 7.4 degs
N225Mpc_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N225Mpc'].values[0]
# N225Mpc error for the upper bound 
N225Mpc_up = N225Mpc_biggest - N225Mpc_mean
# N225Mpc error for the lower bound 
N225Mpc_low = N225Mpc_mean - N225Mpc_smallest

# N200Mpc at theta_c = 16.8 degs
N200Mpc_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N200Mpc'].values[0]
# N200Mpc at theta_c = 10.6 degs
N200Mpc_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N200Mpc'].values[0]
# N200Mpc at theta_c = 7.4 degs
N200Mpc_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N200Mpc'].values[0]
# N200Mpc error for the upper bound 
N200Mpc_up = N200Mpc_biggest - N200Mpc_mean
# N200Mpc error for the lower bound 
N200Mpc_low = N200Mpc_mean - N200Mpc_smallest


# the list of full angle for the coincidence
C_NSNS_list_O3_mean = np.array([N120Mpc_mean, N120Mpc_mean, N65Mpc_mean]) # for coincidence NSNS 
C_NSBH_list_O3_mean = np.array([N190Mpc_mean, N190Mpc_mean, N105Mpc_mean]) # for coincidence NSBH 
C_NSNS_list_O4_mean = np.array([N190Mpc_mean, N190Mpc_mean, N65Mpc_mean]) # for coincidence NSNS 
C_NSBH_list_O4_mean = np.array([N300Mpc_mean, N300Mpc_mean, N105Mpc_mean]) # for coincidence NSBH 
C_NSNS_list_O4_K_mean = np.array([N190Mpc_mean, N190Mpc_mean, N65Mpc_mean, N40Mpc_mean]) # for coincidence NSNS 
C_NSBH_list_O4_K_mean = np.array([N300Mpc_mean, N300Mpc_mean, N105Mpc_mean, N65Mpc_mean]) # for coincidence NSBH 
C_NSNS_list_design_mean = np.array([N190Mpc_mean, N190Mpc_mean, N125Mpc_mean, N140Mpc_mean]) # for coincidence NSNS 
C_NSBH_list_design_mean = np.array([N300Mpc_mean, N300Mpc_mean, N200Mpc_mean, N225Mpc_mean]) # for coincidence NSBH 


# the list of the errors for the upper bound 
C_NSNS_list_O3_up = np.array([N120Mpc_up, N120Mpc_up, N65Mpc_up]) # for coincidence NSNS 
C_NSBH_list_O3_up = np.array([N190Mpc_up, N190Mpc_up, N105Mpc_up]) # for coincidence NSBH 
C_NSNS_list_O4_up = np.array([N190Mpc_up, N190Mpc_up, N65Mpc_up]) # for coincidence NSNS 
C_NSBH_list_O4_up = np.array([N300Mpc_up, N300Mpc_up, N105Mpc_up]) # for coincidence NSBH 
C_NSNS_list_O4_K_up = np.array([N190Mpc_up, N190Mpc_up, N65Mpc_up, N40Mpc_up]) # for coincidence NSNS 
C_NSBH_list_O4_K_up = np.array([N300Mpc_up, N300Mpc_up, N105Mpc_up, N65Mpc_up]) # for coincidence NSBH 
C_NSNS_list_design_up = np.array([N190Mpc_up, N190Mpc_up, N125Mpc_up, N140Mpc_up]) # for coincidence NSNS 
C_NSBH_list_design_up = np.array([N300Mpc_up, N300Mpc_up, N200Mpc_up, N225Mpc_up]) # for coincidence NSBH 


# the list of the errors for the lower bound 
C_NSNS_list_O3_low = np.array([N120Mpc_low, N120Mpc_low, N65Mpc_low]) # for coincidence NSNS 
C_NSBH_list_O3_low = np.array([N190Mpc_low, N190Mpc_low, N105Mpc_low]) # for coincidence NSBH 
C_NSNS_list_O4_low = np.array([N190Mpc_low, N190Mpc_low, N65Mpc_low]) # for coincidence NSNS 
C_NSBH_list_O4_low = np.array([N300Mpc_low, N300Mpc_low, N105Mpc_low]) # for coincidence NSBH 
C_NSNS_list_O4_K_low = np.array([N190Mpc_low, N190Mpc_low, N65Mpc_low, N40Mpc_low]) # for coincidence NSNS 
C_NSBH_list_O4_K_low = np.array([N300Mpc_low, N300Mpc_low, N105Mpc_low, N65Mpc_low]) # for coincidence NSBH 
C_NSNS_list_design_low = np.array([N190Mpc_low, N190Mpc_low, N125Mpc_low, N140Mpc_low]) # for coincidence NSNS 
C_NSBH_list_design_low = np.array([N300Mpc_low, N300Mpc_low, N200Mpc_low, N225Mpc_low]) # for coincidence NSBH 






def C_i(N_c_mean, N_c_up, N_c_low, Q, r_i, Dem, fem):
	'''
	The number of coincident events with a GW detector and a EM detector for either NSNS or NSBH
	per year
	USAGE: 
		For Fermi 
		C_NSNS_ave_O3, C_NSNS_err_O3 = C_i(C_NSNS_list_O3_ave, C_NSNS_list_O3_std, 3, r_nsns, 0.8, 0.75) # for NSNS in O3 with Fermi 
		C_NSBH_ave_O3, C_NSBH_err_O3 = C_i(C_NSBH_list_O3_ave, C_NSBH_list_O3_std, 3, r_nsbh, 0.8, 0.75) # for NSNS in O3 with Fermi 

		C_NSNS_ave_O4, C_NSNS_err_O4 = C_i(C_NSNS_list_O4_ave, C_NSNS_list_O4_std, 3, r_nsns, 0.8, 0.75) # for NSNS in O4 with Fermi 
		C_NSBH_ave_O4, C_NSBH_err_O4 = C_i(C_NSBH_list_O4_ave, C_NSBH_list_O4_std, 3, r_nsbh, 0.8, 0.75) # for NSNS in O4 with Fermi 

		C_NSNS_ave_O4_K, C_NSNS_err_O4_K = C_i(C_NSNS_list_O4_K_ave, C_NSNS_list_O4_K_std, 4, r_nsns, 0.8, 0.75) # for NSNS in O4 including K with Fermi 
		C_NSBH_ave_O4_K, C_NSBH_err_O4_K = C_i(C_NSBH_list_O4_K_ave, C_NSBH_list_O4_K_std, 4, r_nsbh, 0.8, 0.75) # for NSNS in O4 including K with Fermi 

		C_NSNS_ave_design, C_NSNS_err_design = C_i(C_NSNS_list_design_ave, C_NSNS_list_design_std, 4, r_nsns, 0.8, 0.75) # for NSNS in design including K with Fermi
		C_NSBH_ave_design, C_NSBH_err_design = C_i(C_NSBH_list_design_ave, C_NSBH_list_design_std, 4, r_nsbh, 0.8, 0.75) # for NSNS in design including K with Fermi 

		For Swift 
		C_NSNS_ave_O3, C_NSNS_err_O3 = C_i(C_NSNS_list_O3_ave, C_NSNS_list_O3_std, 3, r_nsns, 0.8, 0.1) # for NSNS in O3 with Swfit  
		C_NSBH_ave_O3, C_NSBH_err_O3 = C_i(C_NSBH_list_O3_ave, C_NSBH_list_O3_std, 3, r_nsbh, 0.8, 0.1) # for NSNS in O3 with Swfit

		C_NSNS_ave_O4, C_NSNS_err_O4 = C_i(C_NSNS_list_O4_ave, C_NSNS_list_O4_std, 3, r_nsns, 0.8, 0.1) # for NSNS in O4 with Swfit 
		C_NSBH_ave_O4, C_NSBH_err_O4 = C_i(C_NSBH_list_O4_ave, C_NSBH_list_O4_std, 3, r_nsbh, 0.8, 0.1) # for NSNS in O4 with Swfit 

		C_NSNS_ave_O4_K, C_NSNS_err_O4_K = C_i(C_NSNS_list_O4_K_ave, C_NSNS_list_O4_K_std, 4, r_nsns, 0.8, 0.1) # for NSNS in O4 including K with Swfit 
		C_NSBH_ave_O4_K, C_NSBH_err_O4_K = C_i(C_NSBH_list_O4_K_ave, C_NSBH_list_O4_K_std, 4, r_nsbh, 0.8, 0.1) # for NSNS in O4 including K with Swfit 

		C_NSNS_ave_design, C_NSNS_err_design = C_i(C_NSNS_list_design_ave, C_NSNS_list_design_std, 4, r_nsns, 0.8, 0.1) # for NSNS in design including K with Swift
		C_NSBH_ave_design, C_NSBH_err_design = C_i(C_NSBH_list_design_ave, C_NSBH_list_design_std, 4, r_nsbh, 0.8, 0.1) # for NSNS in design including K with Swift 

		
	'''
	duty_cycle_N_com = []
	duty_cycle_N_com_up = []
	duty_cycle_N_com_low = []

	# conincident for any two interferometr and coincident for any 3 interferomter, any i, ,,, any Q
	# the k starts with 2 ends the total number of interferometers
	for k in range(2, Q+1):
		# all the combinations of any k interferometers out of Q 
		comb = combinations(N_c_mean, k)
		# convert to numpy array
		comb_inter = np.array(list(comb))
		#print('The combination of any %.f inteferometers out of %.f is' % (k, Q))
		#print(comb_inter)
		# looping with all the possible combination of k interferomters out of Q
		for j in range(len(comb_inter)):
			second_large_range = np.sort(comb_inter[j, :])[-2]
			# find the index in N_c_ave with components second_large_range 
			# for getting the error 
			#print('The second highest range among the combination is %f' % second_large_range)
			ind = np.where(N_c_mean == second_large_range)[0][0]
			#print('the index of the list = %i with the components = the second highest ' %ind)
			second_large_range_up = N_c_up[ind]
			second_large_range_low = N_c_low[ind]
			#print('the corresponding upper error is %f' % second_large_range_up)
			#print('the corresponding lower error is %f' % second_large_range_low)
			# the combined duty cycle and the second highest range among the combination
			duty_cycle_N_com.append(D_single**k*(1-D_single)**(Q-k)*second_large_range)
			duty_cycle_N_com_up.append((D_single**k*(1-D_single)**(Q-k)*second_large_range_up)**2)
			duty_cycle_N_com_low.append((D_single**k*(1-D_single)**(Q-k)*second_large_range_low)**2)
	#print('all the combination of duty cycle time Nofz are ')
	#print(duty_cycle_N_com)
	#print('all the combinations of the errors are ')
	#print(duty_cycle_N_com_up, duty_cycle_N_com_low)
	result_ave = r_i*Dem*fem/(f_R*f_FOV*0.8*T)*sum(duty_cycle_N_com)
	result_up = r_i*Dem*fem/(f_R*f_FOV*0.8*T)*(sum(duty_cycle_N_com_up))**0.5
	result_low = r_i*Dem*fem/(f_R*f_FOV*0.8*T)*(sum(duty_cycle_N_com_low))**0.5
	return result_ave, result_up, result_low





#+++++++++++++++++++=
# For Fermi
#++++++++++++++++++++
print('+'*50)


# the nsns and nsbh coincident events with Fermi in O2
C_NSNS_NSBH_O2_mean = (117./365.)*0.8*0.75/(f_R*f_FOV*0.8*T)*(r_nsns*N88Mpc_mean + r_nsbh*N140Mpc_mean)
C_NSNS_NSBH_O2_up = (117./365.)*0.8*0.75/(f_R*f_FOV*0.8*T)*((r_nsns*N88Mpc_up)**2 + (r_nsbh*N140Mpc_up)**2)**0.5
C_NSNS_NSBH_O2_low = (117./365.)*0.8*0.75/(f_R*f_FOV*0.8*T)*((r_nsns*N88Mpc_low)**2 + (r_nsbh*N140Mpc_low)**2)**0.5
print('The coincident events with LH and Fermi in O2= %.3f - %.3f - %.3f' % (C_NSNS_NSBH_O2_mean - C_NSNS_NSBH_O2_low, C_NSNS_NSBH_O2_mean, C_NSNS_NSBH_O2_mean + C_NSNS_NSBH_O2_up))

# the nsns and nsbh number of coincident events in O3
# the number of NSNS events in O3, seen at least two detectors in three-detector network with Fermi
C_NSNS_mean_O3, C_NSNS_up_O3, C_NSNS_low_O3 = C_i(C_NSNS_list_O3_mean, C_NSNS_list_O3_up,  C_NSNS_list_O3_low, 3, r_nsns, 0.8, 0.75) # for NSNS in O3 with Fermi 
# the number of NSBH events in O3
C_NSBH_mean_O3, C_NSBH_up_O3, C_NSBH_low_O3 = C_i(C_NSBH_list_O3_mean, C_NSBH_list_O3_up, C_NSBH_list_O3_low, 3, r_nsbh, 0.8, 0.75) # for NSNS in O3 with Fermi  
# the number of both NSNS and NSBH events
C_nsns_nsbh_g_mean_O3  = C_NSNS_mean_O3  + C_NSBH_mean_O3 
C_nsns_nsbh_g_up_O3 = (C_NSNS_up_O3**2 + C_NSBH_up_O3**2)**0.5 
C_nsns_nsbh_g_low_O3 = (C_NSNS_low_O3**2 + C_NSBH_low_O3**2)**0.5 

print('='*50)
print('The coincident events with Fermi in O3 = %f - %f - %f' % (C_nsns_nsbh_g_mean_O3 - C_nsns_nsbh_g_low_O3, C_nsns_nsbh_g_mean_O3, C_nsns_nsbh_g_mean_O3 + C_nsns_nsbh_g_up_O3))



# O4 
# the number of NSNS events in O4, seen at least two detectors in three-detector network with Fermi
C_NSNS_mean_O4, C_NSNS_up_O4, C_NSNS_low_O4 = C_i(C_NSNS_list_O4_mean, C_NSNS_list_O4_up, C_NSNS_list_O4_low,3, r_nsns, 0.8, 0.75) # for NSNS in O4 with Fermi 
# the number of NSBH events in O4
C_NSBH_mean_O4, C_NSBH_up_O4, C_NSBH_low_O4 = C_i(C_NSBH_list_O4_mean, C_NSBH_list_O4_up, C_NSBH_list_O4_low, 3, r_nsbh, 0.8, 0.75) # for NSNS in O4 with Fermi 
# the number of both NSNS and NSBH events
C_nsns_nsbh_g_mean_O4  = C_NSNS_mean_O4  + C_NSBH_mean_O4 
C_nsns_nsbh_g_up_O4 = (C_NSNS_up_O4**2 + C_NSBH_up_O4**2)**0.5 
C_nsns_nsbh_g_low_O4 = (C_NSNS_low_O4**2 + C_NSBH_low_O4**2)**0.5 

print('='*50)
print('The coincident events with Fermi in O4 = %f - %f - %f' % (C_nsns_nsbh_g_mean_O4 - C_nsns_nsbh_g_low_O4, C_nsns_nsbh_g_mean_O4, C_nsns_nsbh_g_mean_O4 + C_nsns_nsbh_g_up_O4))


# O4 icluding KAGRA
# the number of NSNS events in O4 including K , seen at least two detectors in three-detector network with Fermi
C_NSNS_mean_O4_K, C_NSNS_up_O4_K, C_NSNS_low_O4_K = C_i(C_NSNS_list_O4_K_mean, C_NSNS_list_O4_K_up, C_NSNS_list_O4_K_low, 4, r_nsns, 0.8, 0.75) # for NSNS in O4 including K with Fermi
 # the number of NSBH events in O4 including K
C_NSBH_mean_O4_K, C_NSBH_up_O4_K, C_NSBH_low_O4_K = C_i(C_NSBH_list_O4_K_mean, C_NSBH_list_O4_K_up, C_NSBH_list_O4_K_low, 4, r_nsbh, 0.8, 0.75) # for NSNS in O4 including K with Fermi 
# the number of both NSNS and NSBH events
C_nsns_nsbh_g_mean_O4_K  = C_NSNS_mean_O4_K  + C_NSBH_mean_O4_K 
C_nsns_nsbh_g_up_O4_K = (C_NSNS_up_O4_K**2 + C_NSBH_up_O4_K**2)**0.5 
C_nsns_nsbh_g_low_O4_K = (C_NSNS_low_O4_K**2 + C_NSBH_low_O4_K**2)**0.5 

print('='*50)
print('The coincident events with Fermi in O4 including K = %f - %f - %f' % (C_nsns_nsbh_g_mean_O4_K - C_nsns_nsbh_g_low_O4_K, C_nsns_nsbh_g_mean_O4_K, C_nsns_nsbh_g_mean_O4_K + C_nsns_nsbh_g_up_O4_K))




# design icluding KAGRA
# the number of NSNS events in O4 including K , seen at least two detectors in three-detector network with Fermi
C_NSNS_mean_design, C_NSNS_up_design, C_NSNS_low_design = C_i(C_NSNS_list_design_mean, C_NSNS_list_design_up, C_NSNS_list_design_low, 4, r_nsns, 0.8, 0.75) # for NSNS in design including K with Fermi
 # the number of NSBH events in O4 including K
C_NSBH_mean_design, C_NSBH_up_design, C_NSBH_low_design = C_i(C_NSBH_list_design_mean, C_NSBH_list_design_up, C_NSBH_list_design_low, 4, r_nsbh, 0.8, 0.75) # for NSNS in design including K with Fermi 
# the number of both NSNS and NSBH events
C_nsns_nsbh_g_mean_design = C_NSNS_mean_design  + C_NSBH_mean_design 
C_nsns_nsbh_g_up_design = (C_NSNS_up_design**2 + C_NSBH_up_design**2)**0.5 
C_nsns_nsbh_g_low_design = (C_NSNS_low_design**2 + C_NSBH_low_design**2)**0.5 

print('='*50)
print('The coincident events with Fermi in design = %f - %f - %f' % (C_nsns_nsbh_g_mean_design - C_nsns_nsbh_g_low_design, C_nsns_nsbh_g_mean_design, C_nsns_nsbh_g_mean_design + C_nsns_nsbh_g_up_design))





# #++++++++++++++++++++
# # For Swift
# #++++++++++++++++++++

print('+'*50)


# the nsns and nsbh coincident events with Swift in O2
C_NSNS_NSBH_O2_mean = (117./365.)*0.8*0.1/(f_R*f_FOV*0.8*T)*(r_nsns*N88Mpc_mean + r_nsbh*N140Mpc_mean)
C_NSNS_NSBH_O2_up = (117./365.)*0.8*0.1/(f_R*f_FOV*0.8*T)*((r_nsns*N88Mpc_up)**2 + (r_nsbh*N140Mpc_up)**2)**0.5
C_NSNS_NSBH_O2_low = (117./365.)*0.8*0.1/(f_R*f_FOV*0.8*T)*((r_nsns*N88Mpc_low)**2 + (r_nsbh*N140Mpc_low)**2)**0.5
print('The coincident events with LH and Swift in O2= %.3f - %.3f - %.3f' % (C_NSNS_NSBH_O2_mean - C_NSNS_NSBH_O2_low, C_NSNS_NSBH_O2_mean, C_NSNS_NSBH_O2_mean + C_NSNS_NSBH_O2_up))
print('I think it is ok')

# the nsns and nsbh number of coincident events in O3
# the number of NSNS events in O3, seen at least two detectors in three-detector network with Swift
C_NSNS_mean_O3, C_NSNS_up_O3, C_NSNS_low_O3 = C_i(C_NSNS_list_O3_mean, C_NSNS_list_O3_up,  C_NSNS_list_O3_low, 3, r_nsns, 0.8, 0.1) # for NSNS in O3 with Swift 
# the number of NSBH events in O3
C_NSBH_mean_O3, C_NSBH_up_O3, C_NSBH_low_O3 = C_i(C_NSBH_list_O3_mean, C_NSBH_list_O3_up, C_NSBH_list_O3_low, 3, r_nsbh, 0.8, 0.1) # for NSNS in O3 with Swift  
# the number of both NSNS and NSBH events
C_nsns_nsbh_g_mean_O3  = C_NSNS_mean_O3  + C_NSBH_mean_O3 
C_nsns_nsbh_g_up_O3 = (C_NSNS_up_O3**2 + C_NSBH_up_O3**2)**0.5 
C_nsns_nsbh_g_low_O3 = (C_NSNS_low_O3**2 + C_NSBH_low_O3**2)**0.5 

print('='*50)
print('The coincident events with Swift in O3 = %f - %f - %f' % (C_nsns_nsbh_g_mean_O3 - C_nsns_nsbh_g_low_O3, C_nsns_nsbh_g_mean_O3, C_nsns_nsbh_g_mean_O3 + C_nsns_nsbh_g_up_O3))



# O4 
# the number of NSNS events in O4, seen at least two detectors in three-detector network with Swift
C_NSNS_mean_O4, C_NSNS_up_O4, C_NSNS_low_O4 = C_i(C_NSNS_list_O4_mean, C_NSNS_list_O4_up, C_NSNS_list_O4_low,3, r_nsns, 0.8, 0.1) # for NSNS in O4 with Swift 
# the number of NSBH events in O4
C_NSBH_mean_O4, C_NSBH_up_O4, C_NSBH_low_O4 = C_i(C_NSBH_list_O4_mean, C_NSBH_list_O4_up, C_NSBH_list_O4_low, 3, r_nsbh, 0.8, 0.1) # for NSNS in O4 with Swift 
# the number of both NSNS and NSBH events
C_nsns_nsbh_g_mean_O4  = C_NSNS_mean_O4  + C_NSBH_mean_O4 
C_nsns_nsbh_g_up_O4 = (C_NSNS_up_O4**2 + C_NSBH_up_O4**2)**0.5 
C_nsns_nsbh_g_low_O4 = (C_NSNS_low_O4**2 + C_NSBH_low_O4**2)**0.5 

print('='*50)
print('The coincident events with Swift in O4 = %f - %f - %f' % (C_nsns_nsbh_g_mean_O4 - C_nsns_nsbh_g_low_O4, C_nsns_nsbh_g_mean_O4, C_nsns_nsbh_g_mean_O4 + C_nsns_nsbh_g_up_O4))


# O4 icluding KAGRA
# the number of NSNS events in O4 including K , seen at least two detectors in three-detector network with Swift
C_NSNS_mean_O4_K, C_NSNS_up_O4_K, C_NSNS_low_O4_K = C_i(C_NSNS_list_O4_K_mean, C_NSNS_list_O4_K_up, C_NSNS_list_O4_K_low, 4, r_nsns, 0.8, 0.1) # for NSNS in O4 including K with Swift
 # the number of NSBH events in O4 including K
C_NSBH_mean_O4_K, C_NSBH_up_O4_K, C_NSBH_low_O4_K = C_i(C_NSBH_list_O4_K_mean, C_NSBH_list_O4_K_up, C_NSBH_list_O4_K_low, 4, r_nsbh, 0.8, 0.1) # for NSNS in O4 including K with Swift 
# the number of both NSNS and NSBH events
C_nsns_nsbh_g_mean_O4_K  = C_NSNS_mean_O4_K  + C_NSBH_mean_O4_K 
C_nsns_nsbh_g_up_O4_K = (C_NSNS_up_O4_K**2 + C_NSBH_up_O4_K**2)**0.5 
C_nsns_nsbh_g_low_O4_K = (C_NSNS_low_O4_K**2 + C_NSBH_low_O4_K**2)**0.5 

print('='*50)
print('The coincident events with Swift in O4 including K = %f - %f - %f' % (C_nsns_nsbh_g_mean_O4_K - C_nsns_nsbh_g_low_O4_K, C_nsns_nsbh_g_mean_O4_K, C_nsns_nsbh_g_mean_O4_K + C_nsns_nsbh_g_up_O4_K))




# design icluding KAGRA
# the number of NSNS events in O4 including K , seen at least two detectors in three-detector network with Swift
C_NSNS_mean_design, C_NSNS_up_design, C_NSNS_low_design = C_i(C_NSNS_list_design_mean, C_NSNS_list_design_up, C_NSNS_list_design_low, 4, r_nsns, 0.8, 0.1) # for NSNS in design including K with Swift
 # the number of NSBH events in O4 including K
C_NSBH_mean_design, C_NSBH_up_design, C_NSBH_low_design = C_i(C_NSBH_list_design_mean, C_NSBH_list_design_up, C_NSBH_list_design_low, 4, r_nsbh, 0.8, 0.1) # for NSNS in design including K with Swift 
# the number of both NSNS and NSBH events
C_nsns_nsbh_g_mean_design = C_NSNS_mean_design  + C_NSBH_mean_design 
C_nsns_nsbh_g_up_design = (C_NSNS_up_design**2 + C_NSBH_up_design**2)**0.5 
C_nsns_nsbh_g_low_design = (C_NSNS_low_design**2 + C_NSBH_low_design**2)**0.5 

print('='*50)
print('The coincident events with Swift in design = %f - %f - %f' % (C_nsns_nsbh_g_mean_design - C_nsns_nsbh_g_low_design, C_nsns_nsbh_g_mean_design, C_nsns_nsbh_g_mean_design + C_nsns_nsbh_g_up_design))





#=====================================
# Off-axis events 
#====================================



# N88Mpc at theta_c = 16.8 degs
N88MpcOffaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N88MpcOffaxis'].values[0]
# N88Mpc at theta_c = 10.6 degs
N88MpcOffaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N88MpcOffaxis'].values[0]
# N88Mpc at theta_c =  7.4 degs
N88MpcOffaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N88MpcOffaxis'].values[0]
# N88Mpc error for the upper bound 
N88MpcOffaxis_up = N88MpcOffaxis_biggest - N88MpcOffaxis_mean
# N88Mpc error for the lower bound 
N88MpcOffaxis_low = N88MpcOffaxis_mean - N88MpcOffaxis_smallest 

# N120Mpc at theta_c = 16.8 degs
N120MpcOffaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N120MpcOffaxis'].values[0]
# N120Mpc at theta_c = 10.6 degs
N120MpcOffaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N120MpcOffaxis'].values[0]
# N120Mpc at theta_c =  7.4 degs
N120MpcOffaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N120MpcOffaxis'].values[0]
# N120Mpc error for the upper bound 
N120MpcOffaxis_up = N120MpcOffaxis_biggest - N120MpcOffaxis_mean
# N120Mpc error for the lower bound 
N120MpcOffaxis_low = N120MpcOffaxis_mean - N120MpcOffaxis_smallest


# N65Mpc at theta_c = 16.8 degs
N65MpcOffaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N65MpcOffaxis'].values[0]
# N65Mpc at theta_c = 10.6 degs
N65MpcOffaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N65MpcOffaxis'].values[0]
# N65Mpc at theta_c = 7.4 degs
N65MpcOffaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N65MpcOffaxis'].values[0]
# N65Mpc error for the upper bound 
N65MpcOffaxis_up = N65MpcOffaxis_biggest - N65MpcOffaxis_mean
# N65Mpc error for the lower bound 
N65MpcOffaxis_low = N65MpcOffaxis_mean - N65MpcOffaxis_smallest

# N190Mpc at theta_c = 16.8 degs
N190MpcOffaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N190MpcOffaxis'].values[0]
# N190Mpc at theta_c = 10.6 degs
N190MpcOffaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N190MpcOffaxis'].values[0]
# N190Mpc at theta_c = 7.4 degs
N190MpcOffaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N190MpcOffaxis'].values[0]
# N190Mpc error for the upper bound 
N190MpcOffaxis_up = N190MpcOffaxis_biggest - N190MpcOffaxis_mean
# N190Mpc error for the lower bound 
N190MpcOffaxis_low = N190MpcOffaxis_mean - N190MpcOffaxis_smallest

# N105Mpc at theta_c = 16.8 degs
N105MpcOffaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N105MpcOffaxis'].values[0]
# N105Mpc at theta_c = 10.6 degs
N105MpcOffaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N105MpcOffaxis'].values[0]
# N105Mpc at theta_c = 7.4 degs
N105MpcOffaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N105MpcOffaxis'].values[0]
# N105Mpc error for the upper bound 
N105MpcOffaxis_up = N105MpcOffaxis_biggest - N105MpcOffaxis_mean
# N105Mpc error for the lower bound 
N105MpcOffaxis_low = N105MpcOffaxis_mean - N105MpcOffaxis_smallest

# N40Mpc at theta_c = 16.8 degs
N40MpcOffaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N40MpcOffaxis'].values[0]
# N40Mpc at theta_c = 10.6 degs
N40MpcOffaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N40MpcOffaxis'].values[0]
# N40Mpc at theta_c = 7.4 degs
N40MpcOffaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N40MpcOffaxis'].values[0]
# N40Mpc error for the upper bound 
N40MpcOffaxis_up = N40MpcOffaxis_biggest - N40MpcOffaxis_mean
# N40Mpc error for the lower bound 
N40MpcOffaxis_low = N40MpcOffaxis_mean - N40MpcOffaxis_smallest

# N300Mpc at theta_c = 16.8 degs
N300MpcOffaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N300MpcOffaxis'].values[0]
# N300Mpc at theta_c = 10.6 degs
N300MpcOffaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N300MpcOffaxis'].values[0]
# N300Mpc at theta_c = 7.4 degs
N300MpcOffaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N300MpcOffaxis'].values[0]
# N300Mpc error for the upper bound 
N300MpcOffaxis_up = N300MpcOffaxis_biggest - N300MpcOffaxis_mean
# N300Mpc error for the lower bound 
N300MpcOffaxis_low = N300MpcOffaxis_mean - N300MpcOffaxis_smallest


# N125Mpc at theta_c = 16.8 degs
N125MpcOffaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N125MpcOffaxis'].values[0]
# N125Mpc at theta_c = 10.6 degs
N125MpcOffaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N125MpcOffaxis'].values[0]
# N125Mpc at theta_c = 7.4 degs
N125MpcOffaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N125MpcOffaxis'].values[0]
# N125Mpc error for the upper bound 
N125MpcOffaxis_up = N125MpcOffaxis_biggest - N125MpcOffaxis_mean
# N125Mpc error for the lower bound 
N125MpcOffaxis_low = N125MpcOffaxis_mean - N125MpcOffaxis_smallest

# N140Mpc at theta_c = 16.8 degs
N140MpcOffaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N140MpcOffaxis'].values[0]
# N140Mpc at theta_c = 10.6 degs
N140MpcOffaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N140MpcOffaxis'].values[0]
# N140Mpc at theta_c = 7.4 degs
N140MpcOffaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N140MpcOffaxis'].values[0]
# N140Mpc error for the upper bound 
N140MpcOffaxis_up = N140MpcOffaxis_biggest - N140MpcOffaxis_mean
# N140Mpc error for the lower bound 
N140MpcOffaxis_low = N140MpcOffaxis_mean - N140MpcOffaxis_smallest

# N225Mpc at theta_c = 16.8 degs
N225MpcOffaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N225MpcOffaxis'].values[0]
# N225Mpc at theta_c = 10.6 degs
N225MpcOffaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N225MpcOffaxis'].values[0]
# N225Mpc at theta_c = 7.4 degs
N225MpcOffaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N225MpcOffaxis'].values[0]
# N225Mpc error for the upper bound 
N225MpcOffaxis_up = N225MpcOffaxis_biggest - N225MpcOffaxis_mean
# N225Mpc error for the lower bound 
N225MpcOffaxis_low = N225MpcOffaxis_mean - N225MpcOffaxis_smallest

# N200Mpc at theta_c = 16.8 degs
N200MpcOffaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N200MpcOffaxis'].values[0]
# N200Mpc at theta_c = 10.6 degs
N200MpcOffaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N200MpcOffaxis'].values[0]
# N200Mpc at theta_c = 7.4 degs
N200MpcOffaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N200MpcOffaxis'].values[0]
# N200Mpc error for the upper bound 
N200MpcOffaxis_up = N200MpcOffaxis_biggest - N200MpcOffaxis_mean
# N200Mpc error for the lower bound 
N200MpcOffaxis_low = N200MpcOffaxis_mean - N200MpcOffaxis_smallest


# the list of full angle for the coincidence
COffaxis_NSNS_list_O3_mean = np.array([N120MpcOffaxis_mean, N120MpcOffaxis_mean, N65MpcOffaxis_mean]) # for coincidence NSNS 
COffaxis_NSBH_list_O3_mean = np.array([N190MpcOffaxis_mean, N190MpcOffaxis_mean, N105MpcOffaxis_mean]) # for coincidence NSBH 
COffaxis_NSNS_list_O4_mean = np.array([N190MpcOffaxis_mean, N190MpcOffaxis_mean, N65MpcOffaxis_mean]) # for coincidence NSNS 
COffaxis_NSBH_list_O4_mean = np.array([N300MpcOffaxis_mean, N300MpcOffaxis_mean, N105MpcOffaxis_mean]) # for coincidence NSBH 
COffaxis_NSNS_list_O4_K_mean = np.array([N190MpcOffaxis_mean, N190MpcOffaxis_mean, N65MpcOffaxis_mean, N40MpcOffaxis_mean]) # for coincidence NSNS 
COffaxis_NSBH_list_O4_K_mean = np.array([N300MpcOffaxis_mean, N300MpcOffaxis_mean, N105MpcOffaxis_mean, N65MpcOffaxis_mean]) # for coincidence NSBH 
COffaxis_NSNS_list_design_mean = np.array([N190MpcOffaxis_mean, N190MpcOffaxis_mean, N125MpcOffaxis_mean, N140MpcOffaxis_mean]) # for coincidence NSNS 
COffaxis_NSBH_list_design_mean = np.array([N300MpcOffaxis_mean, N300MpcOffaxis_mean, N200MpcOffaxis_mean, N225MpcOffaxis_mean]) # for coincidence NSBH 


# the list of the errors for the upper bound 
COffaxis_NSNS_list_O3_up = np.array([N120MpcOffaxis_up, N120MpcOffaxis_up, N65MpcOffaxis_up]) # for coincidence NSNS 
COffaxis_NSBH_list_O3_up = np.array([N190MpcOffaxis_up, N190MpcOffaxis_up, N105MpcOffaxis_up]) # for coincidence NSBH 
COffaxis_NSNS_list_O4_up = np.array([N190MpcOffaxis_up, N190MpcOffaxis_up, N65MpcOffaxis_up]) # for coincidence NSNS 
COffaxis_NSBH_list_O4_up = np.array([N300MpcOffaxis_up, N300MpcOffaxis_up, N105MpcOffaxis_up]) # for coincidence NSBH 
COffaxis_NSNS_list_O4_K_up = np.array([N190MpcOffaxis_up, N190MpcOffaxis_up, N65MpcOffaxis_up, N40MpcOffaxis_up]) # for coincidence NSNS 
COffaxis_NSBH_list_O4_K_up = np.array([N300MpcOffaxis_up, N300MpcOffaxis_up, N105MpcOffaxis_up, N65MpcOffaxis_up]) # for coincidence NSBH 
COffaxis_NSNS_list_design_up = np.array([N190MpcOffaxis_up, N190MpcOffaxis_up, N125MpcOffaxis_up, N140MpcOffaxis_up]) # for coincidence NSNS 
COffaxis_NSBH_list_design_up = np.array([N300MpcOffaxis_up, N300MpcOffaxis_up, N200MpcOffaxis_up, N225MpcOffaxis_up]) # for coincidence NSBH 


# the list of the errors for the lower bound 
COffaxis_NSNS_list_O3_low = np.array([N120MpcOffaxis_low, N120MpcOffaxis_low, N65MpcOffaxis_low]) # for coincidence NSNS 
COffaxis_NSBH_list_O3_low = np.array([N190MpcOffaxis_low, N190MpcOffaxis_low, N105MpcOffaxis_low]) # for coincidence NSBH 
COffaxis_NSNS_list_O4_low = np.array([N190MpcOffaxis_low, N190MpcOffaxis_low, N65MpcOffaxis_low]) # for coincidence NSNS 
COffaxis_NSBH_list_O4_low = np.array([N300MpcOffaxis_low, N300MpcOffaxis_low, N105MpcOffaxis_low]) # for coincidence NSBH 
COffaxis_NSNS_list_O4_K_low = np.array([N190MpcOffaxis_low, N190MpcOffaxis_low, N65MpcOffaxis_low, N40MpcOffaxis_low]) # for coincidence NSNS 
COffaxis_NSBH_list_O4_K_low = np.array([N300MpcOffaxis_low, N300MpcOffaxis_low, N105MpcOffaxis_low, N65MpcOffaxis_low]) # for coincidence NSBH 
COffaxis_NSNS_list_design_low = np.array([N190MpcOffaxis_low, N190MpcOffaxis_low, N125MpcOffaxis_low, N140MpcOffaxis_low]) # for coincidence NSNS 
COffaxis_NSBH_list_design_low = np.array([N300MpcOffaxis_low, N300MpcOffaxis_low, N200MpcOffaxis_low, N225MpcOffaxis_low]) # for coincidence NSBH 



#++++++++++++++++++++++++++++++
# For Fermi only Off axis
#++++++++++++++++++++++++++++++++

print('+'*50)


# the nsns and nsbh coincident (only off-axis) events with Fermi in O2
COffaxis_NSNS_NSBH_O2_mean = (117./365.)*0.8*0.75/(f_R*f_FOV*0.8*T)*(r_nsns*N88MpcOffaxis_mean + r_nsbh*N140MpcOffaxis_mean)
COffaxis_NSNS_NSBH_O2_up = (117./365.)*0.8*0.75/(f_R*f_FOV*0.8*T)*((r_nsns*N88MpcOffaxis_up)**2 + (r_nsbh*N140MpcOffaxis_up)**2)**0.5
COffaxis_NSNS_NSBH_O2_low = (117./365.)*0.8*0.75/(f_R*f_FOV*0.8*T)*((r_nsns*N88MpcOffaxis_low)**2 + (r_nsbh*N140MpcOffaxis_low)**2)**0.5
print('The coincident (only off-axis) events with LH and Fermi in O2 = %.3f - %.3f - %.3f' % (COffaxis_NSNS_NSBH_O2_mean - COffaxis_NSNS_NSBH_O2_low, COffaxis_NSNS_NSBH_O2_mean, COffaxis_NSNS_NSBH_O2_mean + COffaxis_NSNS_NSBH_O2_up))

# the nsns and nsbh number of coincident (only off-axis) events in O3
# the number of NSNS events in O3, seen at least two detectors in three-detector network with Fermi
COffaxis_NSNS_mean_O3, COffaxis_NSNS_up_O3, COffaxis_NSNS_low_O3 = C_i(COffaxis_NSNS_list_O3_mean, COffaxis_NSNS_list_O3_up,  COffaxis_NSNS_list_O3_low, 3, r_nsns, 0.8, 0.75) # for NSNS in O3 with Fermi 
# the number of NSBH events in O3
COffaxis_NSBH_mean_O3, COffaxis_NSBH_up_O3, COffaxis_NSBH_low_O3 = C_i(COffaxis_NSBH_list_O3_mean, COffaxis_NSBH_list_O3_up, COffaxis_NSBH_list_O3_low, 3, r_nsbh, 0.8, 0.75) # for NSNS in O3 with Fermi  
# the number of both NSNS and NSBH events
COffaxis_nsns_nsbh_g_mean_O3  = COffaxis_NSNS_mean_O3  + COffaxis_NSBH_mean_O3 
COffaxis_nsns_nsbh_g_up_O3 = (COffaxis_NSNS_up_O3**2 + COffaxis_NSBH_up_O3**2)**0.5 
COffaxis_nsns_nsbh_g_low_O3 = (COffaxis_NSNS_low_O3**2 + COffaxis_NSBH_low_O3**2)**0.5 

print('='*50)
print('The coincident (only off-axis) events with Fermi in O3 = %f - %f - %f' % (COffaxis_nsns_nsbh_g_mean_O3 - COffaxis_nsns_nsbh_g_low_O3, COffaxis_nsns_nsbh_g_mean_O3, COffaxis_nsns_nsbh_g_mean_O3 + COffaxis_nsns_nsbh_g_up_O3))



# O4 
# the number of NSNS events in O4, seen at least two detectors in three-detector network with Fermi
COffaxis_NSNS_mean_O4, COffaxis_NSNS_up_O4, COffaxis_NSNS_low_O4 = C_i(COffaxis_NSNS_list_O4_mean, COffaxis_NSNS_list_O4_up, COffaxis_NSNS_list_O4_low,3, r_nsns, 0.8, 0.75) # for NSNS in O4 with Fermi 
# the number of NSBH events in O4
COffaxis_NSBH_mean_O4, COffaxis_NSBH_up_O4, COffaxis_NSBH_low_O4 = C_i(COffaxis_NSBH_list_O4_mean, COffaxis_NSBH_list_O4_up, COffaxis_NSBH_list_O4_low, 3, r_nsbh, 0.8, 0.75) # for NSNS in O4 with Fermi 
# the number of both NSNS and NSBH events
COffaxis_nsns_nsbh_g_mean_O4  = COffaxis_NSNS_mean_O4  + COffaxis_NSBH_mean_O4 
COffaxis_nsns_nsbh_g_up_O4 = (COffaxis_NSNS_up_O4**2 + COffaxis_NSBH_up_O4**2)**0.5 
COffaxis_nsns_nsbh_g_low_O4 = (COffaxis_NSNS_low_O4**2 + COffaxis_NSBH_low_O4**2)**0.5 

print('='*50)
print('The coincident (only off-axis) events with Fermi in O4 = %f - %f - %f' % (COffaxis_nsns_nsbh_g_mean_O4 - COffaxis_nsns_nsbh_g_low_O4, COffaxis_nsns_nsbh_g_mean_O4, COffaxis_nsns_nsbh_g_mean_O4 + COffaxis_nsns_nsbh_g_up_O4))


# O4 icluding KAGRA
# the number of NSNS events in O4 including K , seen at least two detectors in three-detector network with Fermi
COffaxis_NSNS_mean_O4_K, COffaxis_NSNS_up_O4_K, COffaxis_NSNS_low_O4_K = C_i(COffaxis_NSNS_list_O4_K_mean, COffaxis_NSNS_list_O4_K_up, COffaxis_NSNS_list_O4_K_low, 4, r_nsns, 0.8, 0.75) # for NSBH in O4 including K with Fermi
 # the number of NSBH events in O4 including K
COffaxis_NSBH_mean_O4_K, COffaxis_NSBH_up_O4_K, COffaxis_NSBH_low_O4_K = C_i(COffaxis_NSBH_list_O4_K_mean, COffaxis_NSBH_list_O4_K_up, COffaxis_NSBH_list_O4_K_low, 4, r_nsbh, 0.8, 0.75) # for NSNS in O4 including K with Fermi 
# the number of both NSNS and NSBH events
COffaxis_nsns_nsbh_g_mean_O4_K  = COffaxis_NSNS_mean_O4_K  + COffaxis_NSBH_mean_O4_K 
COffaxis_nsns_nsbh_g_up_O4_K = (COffaxis_NSNS_up_O4_K**2 + COffaxis_NSBH_up_O4_K**2)**0.5 
COffaxis_nsns_nsbh_g_low_O4_K = (COffaxis_NSNS_low_O4_K**2 + COffaxis_NSBH_low_O4_K**2)**0.5 

print('='*50)
print('The coincident (only off-axis) events with Fermi in O4 including K = %f - %f - %f' % (COffaxis_nsns_nsbh_g_mean_O4_K - COffaxis_nsns_nsbh_g_low_O4_K, COffaxis_nsns_nsbh_g_mean_O4_K, COffaxis_nsns_nsbh_g_mean_O4_K + COffaxis_nsns_nsbh_g_up_O4_K))




# design icluding KAGRA
# the number of NSNS events in O4 including K , seen at least two detectors in three-detector network with Fermi
COffaxis_NSNS_mean_design, COffaxis_NSNS_up_design, COffaxis_NSNS_low_design = C_i(COffaxis_NSNS_list_design_mean, COffaxis_NSNS_list_design_up, COffaxis_NSNS_list_design_low, 4, r_nsns, 0.8, 0.75) # for NSNS in design including K with Fermi
 # the number of NSBH events in O4 including K
COffaxis_NSBH_mean_design, COffaxis_NSBH_up_design, COffaxis_NSBH_low_design = C_i(COffaxis_NSBH_list_design_mean, COffaxis_NSBH_list_design_up, COffaxis_NSBH_list_design_low, 4, r_nsbh, 0.8, 0.75) # for NSNS in design including K with Fermi 
# the number of both NSNS and NSBH events
COffaxis_nsns_nsbh_g_mean_design = COffaxis_NSNS_mean_design  + COffaxis_NSBH_mean_design 
COffaxis_nsns_nsbh_g_up_design = (COffaxis_NSNS_up_design**2 + COffaxis_NSBH_up_design**2)**0.5 
COffaxis_nsns_nsbh_g_low_design = (COffaxis_NSNS_low_design**2 + COffaxis_NSBH_low_design**2)**0.5 

print('='*50)
print('The coincident (only off-axis) events with Fermi in design = %f - %f - %f' % (COffaxis_nsns_nsbh_g_mean_design - COffaxis_nsns_nsbh_g_low_design, COffaxis_nsns_nsbh_g_mean_design, COffaxis_nsns_nsbh_g_mean_design + COffaxis_nsns_nsbh_g_up_design))


#++++++++++++++++++++++++++++++
# For Swift only Off-axis
#++++++++++++++++++++++++++++++++

print('+'*50)


# the nsns and nsbh coincident (only off-axis) events with Swift in O2
COffaxis_NSNS_NSBH_O2_mean = (117./365.)*0.8*0.1/(f_R*f_FOV*0.8*T)*(r_nsns*N88MpcOffaxis_mean + r_nsbh*N140MpcOffaxis_mean)
COffaxis_NSNS_NSBH_O2_up = (117./365.)*0.8*0.1/(f_R*f_FOV*0.8*T)*((r_nsns*N88MpcOffaxis_up)**2 + (r_nsbh*N140MpcOffaxis_up)**2)**0.5
COffaxis_NSNS_NSBH_O2_low = (117./365.)*0.8*0.1/(f_R*f_FOV*0.8*T)*((r_nsns*N88MpcOffaxis_low)**2 + (r_nsbh*N140MpcOffaxis_low)**2)**0.5
print('The coincident (only off-axis) events with LH and Swift in O2 = %.3f - %.3f - %.3f' % (COffaxis_NSNS_NSBH_O2_mean - COffaxis_NSNS_NSBH_O2_low, COffaxis_NSNS_NSBH_O2_mean, COffaxis_NSNS_NSBH_O2_mean + COffaxis_NSNS_NSBH_O2_up))

# the nsns and nsbh number of coincident (only off-axis) events in O3
# the number of NSNS events in O3, seen at least two detectors in three-detector network with Fermi
COffaxis_NSNS_mean_O3, COffaxis_NSNS_up_O3, COffaxis_NSNS_low_O3 = C_i(COffaxis_NSNS_list_O3_mean, COffaxis_NSNS_list_O3_up,  COffaxis_NSNS_list_O3_low, 3, r_nsns, 0.8, 0.1) # for NSNS in O3 with Swift 
# the number of NSBH events in O3
COffaxis_NSBH_mean_O3, COffaxis_NSBH_up_O3, COffaxis_NSBH_low_O3 = C_i(COffaxis_NSBH_list_O3_mean, COffaxis_NSBH_list_O3_up, COffaxis_NSBH_list_O3_low, 3, r_nsbh, 0.8, 0.1) # for NSNS in O3 with Swift  
# the number of both NSNS and NSBH events
COffaxis_nsns_nsbh_g_mean_O3  = COffaxis_NSNS_mean_O3  + COffaxis_NSBH_mean_O3 
COffaxis_nsns_nsbh_g_up_O3 = (COffaxis_NSNS_up_O3**2 + COffaxis_NSBH_up_O3**2)**0.5 
COffaxis_nsns_nsbh_g_low_O3 = (COffaxis_NSNS_low_O3**2 + COffaxis_NSBH_low_O3**2)**0.5 

print('='*50)
print('The coincident (only off-axis) events with Swift in O3 = %f - %f - %f' % (COffaxis_nsns_nsbh_g_mean_O3 - COffaxis_nsns_nsbh_g_low_O3, COffaxis_nsns_nsbh_g_mean_O3, COffaxis_nsns_nsbh_g_mean_O3 + COffaxis_nsns_nsbh_g_up_O3))



# O4 
# the number of NSNS events in O4, seen at least two detectors in three-detector network with Swift
COffaxis_NSNS_mean_O4, COffaxis_NSNS_up_O4, COffaxis_NSNS_low_O4 = C_i(COffaxis_NSNS_list_O4_mean, COffaxis_NSNS_list_O4_up, COffaxis_NSNS_list_O4_low,3, r_nsns, 0.8, 0.1) # for NSNS in O4 with Swift 
# the number of NSBH events in O4
COffaxis_NSBH_mean_O4, COffaxis_NSBH_up_O4, COffaxis_NSBH_low_O4 = C_i(COffaxis_NSBH_list_O4_mean, COffaxis_NSBH_list_O4_up, COffaxis_NSBH_list_O4_low, 3, r_nsbh, 0.8, 0.1) # for NSNS in O4 with Swift 
# the number of both NSNS and NSBH events
COffaxis_nsns_nsbh_g_mean_O4  = COffaxis_NSNS_mean_O4  + COffaxis_NSBH_mean_O4 
COffaxis_nsns_nsbh_g_up_O4 = (COffaxis_NSNS_up_O4**2 + COffaxis_NSBH_up_O4**2)**0.5 
COffaxis_nsns_nsbh_g_low_O4 = (COffaxis_NSNS_low_O4**2 + COffaxis_NSBH_low_O4**2)**0.5 

print('='*50)
print('The coincident (only off-axis) events with Swift in O4 = %f - %f - %f' % (COffaxis_nsns_nsbh_g_mean_O4 - COffaxis_nsns_nsbh_g_low_O4, COffaxis_nsns_nsbh_g_mean_O4, COffaxis_nsns_nsbh_g_mean_O4 + COffaxis_nsns_nsbh_g_up_O4))


# O4 icluding KAGRA
# the number of NSNS events in O4 including K , seen at least two detectors in three-detector network with Fermi
COffaxis_NSNS_mean_O4_K, COffaxis_NSNS_up_O4_K, COffaxis_NSNS_low_O4_K = C_i(COffaxis_NSNS_list_O4_K_mean, COffaxis_NSNS_list_O4_K_up, COffaxis_NSNS_list_O4_K_low, 4, r_nsns, 0.8, 0.1) # for NSBH in O4 including K with Swift
 # the number of NSBH events in O4 including K
COffaxis_NSBH_mean_O4_K, COffaxis_NSBH_up_O4_K, COffaxis_NSBH_low_O4_K = C_i(COffaxis_NSBH_list_O4_K_mean, COffaxis_NSBH_list_O4_K_up, COffaxis_NSBH_list_O4_K_low, 4, r_nsbh, 0.8, 0.1) # for NSNS in O4 including K with Swift 
# the number of both NSNS and NSBH events
COffaxis_nsns_nsbh_g_mean_O4_K  = COffaxis_NSNS_mean_O4_K  + COffaxis_NSBH_mean_O4_K 
COffaxis_nsns_nsbh_g_up_O4_K = (COffaxis_NSNS_up_O4_K**2 + COffaxis_NSBH_up_O4_K**2)**0.5 
COffaxis_nsns_nsbh_g_low_O4_K = (COffaxis_NSNS_low_O4_K**2 + COffaxis_NSBH_low_O4_K**2)**0.5 

print('='*50)
print('The coincident (only off-axis) events with Swift in O4 including K = %f - %f - %f' % (COffaxis_nsns_nsbh_g_mean_O4_K - COffaxis_nsns_nsbh_g_low_O4_K, COffaxis_nsns_nsbh_g_mean_O4_K, COffaxis_nsns_nsbh_g_mean_O4_K + COffaxis_nsns_nsbh_g_up_O4_K))




# design icluding KAGRA
# the number of NSNS events in O4 including K , seen at least two detectors in three-detector network with Fermi
COffaxis_NSNS_mean_design, COffaxis_NSNS_up_design, COffaxis_NSNS_low_design = C_i(COffaxis_NSNS_list_design_mean, COffaxis_NSNS_list_design_up, COffaxis_NSNS_list_design_low, 4, r_nsns, 0.8, 0.1) # for NSNS in design including K with Swift
 # the number of NSBH events in O4 including K
COffaxis_NSBH_mean_design, COffaxis_NSBH_up_design, COffaxis_NSBH_low_design = C_i(COffaxis_NSBH_list_design_mean, COffaxis_NSBH_list_design_up, COffaxis_NSBH_list_design_low, 4, r_nsbh, 0.8, 0.1) # for NSNS in design including K with Swift
# the number of both NSNS and NSBH events
COffaxis_nsns_nsbh_g_mean_design = COffaxis_NSNS_mean_design  + COffaxis_NSBH_mean_design 
COffaxis_nsns_nsbh_g_up_design = (COffaxis_NSNS_up_design**2 + COffaxis_NSBH_up_design**2)**0.5 
COffaxis_nsns_nsbh_g_low_design = (COffaxis_NSNS_low_design**2 + COffaxis_NSBH_low_design**2)**0.5 

print('='*50)
print('The coincident (only off-axis) events with Swift in design = %f - %f - %f' % (COffaxis_nsns_nsbh_g_mean_design - COffaxis_nsns_nsbh_g_low_design, COffaxis_nsns_nsbh_g_mean_design, COffaxis_nsns_nsbh_g_mean_design + COffaxis_nsns_nsbh_g_up_design))




#=====================================
# On-axis events  list
#====================================



# N88Mpc at theta_c = 16.8 degs
N88MpcOnaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N88MpcOnaxis'].values[0]
# N88Mpc at theta_c = 10.6 degs
N88MpcOnaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N88MpcOnaxis'].values[0]
# N88Mpc at theta_c =  7.4 degs
N88MpcOnaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N88MpcOnaxis'].values[0]
# N88Mpc error for the upper bound 
N88MpcOnaxis_up = N88MpcOnaxis_biggest - N88MpcOnaxis_mean
# N88Mpc error for the lower bound 
N88MpcOnaxis_low = N88MpcOnaxis_mean - N88MpcOnaxis_smallest 

# N120Mpc at theta_c = 16.8 degs
N120MpcOnaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N120MpcOnaxis'].values[0]
# N120Mpc at theta_c = 10.6 degs
N120MpcOnaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N120MpcOnaxis'].values[0]
# N120Mpc at theta_c =  7.4 degs
N120MpcOnaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N120MpcOnaxis'].values[0]
# N120Mpc error for the upper bound 
N120MpcOnaxis_up = N120MpcOnaxis_biggest - N120MpcOnaxis_mean
# N120Mpc error for the lower bound 
N120MpcOnaxis_low = N120MpcOnaxis_mean - N120MpcOnaxis_smallest


# N65Mpc at theta_c = 16.8 degs
N65MpcOnaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N65MpcOnaxis'].values[0]
# N65Mpc at theta_c = 10.6 degs
N65MpcOnaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N65MpcOnaxis'].values[0]
# N65Mpc at theta_c = 7.4 degs
N65MpcOnaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N65MpcOnaxis'].values[0]
# N65Mpc error for the upper bound 
N65MpcOnaxis_up = N65MpcOnaxis_biggest - N65MpcOnaxis_mean
# N65Mpc error for the lower bound 
N65MpcOnaxis_low = N65MpcOnaxis_mean - N65MpcOnaxis_smallest

# N190Mpc at theta_c = 16.8 degs
N190MpcOnaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N190MpcOnaxis'].values[0]
# N190Mpc at theta_c = 10.6 degs
N190MpcOnaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N190MpcOnaxis'].values[0]
# N190Mpc at theta_c = 7.4 degs
N190MpcOnaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N190MpcOnaxis'].values[0]
# N190Mpc error for the upper bound 
N190MpcOnaxis_up = N190MpcOnaxis_biggest - N190MpcOnaxis_mean
# N190Mpc error for the lower bound 
N190MpcOnaxis_low = N190MpcOnaxis_mean - N190MpcOnaxis_smallest

# N105Mpc at theta_c = 16.8 degs
N105MpcOnaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N105MpcOnaxis'].values[0]
# N105Mpc at theta_c = 10.6 degs
N105MpcOnaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N105MpcOnaxis'].values[0]
# N105Mpc at theta_c = 7.4 degs
N105MpcOnaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N105MpcOnaxis'].values[0]
# N105Mpc error for the upper bound 
N105MpcOnaxis_up = N105MpcOnaxis_biggest - N105MpcOnaxis_mean
# N105Mpc error for the lower bound 
N105MpcOnaxis_low = N105MpcOnaxis_mean - N105MpcOnaxis_smallest

# N40Mpc at theta_c = 16.8 degs
N40MpcOnaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N40MpcOnaxis'].values[0]
# N40Mpc at theta_c = 10.6 degs
N40MpcOnaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N40MpcOnaxis'].values[0]
# N40Mpc at theta_c = 7.4 degs
N40MpcOnaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N40MpcOnaxis'].values[0]
# N40Mpc error for the upper bound 
N40MpcOnaxis_up = N40MpcOnaxis_biggest - N40MpcOnaxis_mean
# N40Mpc error for the lower bound 
N40MpcOnaxis_low = N40MpcOnaxis_mean - N40MpcOnaxis_smallest

# N300Mpc at theta_c = 16.8 degs
N300MpcOnaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N300MpcOnaxis'].values[0]
# N300Mpc at theta_c = 10.6 degs
N300MpcOnaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N300MpcOnaxis'].values[0]
# N300Mpc at theta_c = 7.4 degs
N300MpcOnaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N300MpcOnaxis'].values[0]
# N300Mpc error for the upper bound 
N300MpcOnaxis_up = N300MpcOnaxis_biggest - N300MpcOnaxis_mean
# N300Mpc error for the lower bound 
N300MpcOnaxis_low = N300MpcOnaxis_mean - N300MpcOnaxis_smallest


# N125Mpc at theta_c = 16.8 degs
N125MpcOnaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N125MpcOnaxis'].values[0]
# N125Mpc at theta_c = 10.6 degs
N125MpcOnaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N125MpcOnaxis'].values[0]
# N125Mpc at theta_c = 7.4 degs
N125MpcOnaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N125MpcOnaxis'].values[0]
# N125Mpc error for the upper bound 
N125MpcOnaxis_up = N125MpcOnaxis_biggest - N125MpcOnaxis_mean
# N125Mpc error for the lower bound 
N125MpcOnaxis_low = N125MpcOnaxis_mean - N125MpcOnaxis_smallest

# N140Mpc at theta_c = 16.8 degs
N140MpcOnaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N140MpcOnaxis'].values[0]
# N140Mpc at theta_c = 10.6 degs
N140MpcOnaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N140MpcOnaxis'].values[0]
# N140Mpc at theta_c = 7.4 degs
N140MpcOnaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N140MpcOnaxis'].values[0]
# N140Mpc error for the upper bound 
N140MpcOnaxis_up = N140MpcOnaxis_biggest - N140MpcOnaxis_mean
# N140Mpc error for the lower bound 
N140MpcOnaxis_low = N140MpcOnaxis_mean - N140MpcOnaxis_smallest

# N225Mpc at theta_c = 16.8 degs
N225MpcOnaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N225MpcOnaxis'].values[0]
# N225Mpc at theta_c = 10.6 degs
N225MpcOnaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N225MpcOnaxis'].values[0]
# N225Mpc at theta_c = 7.4 degs
N225MpcOnaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N225MpcOnaxis'].values[0]
# N225Mpc error for the upper bound 
N225MpcOnaxis_up = N225MpcOnaxis_biggest - N225MpcOnaxis_mean
# N225Mpc error for the lower bound 
N225MpcOnaxis_low = N225MpcOnaxis_mean - N225MpcOnaxis_smallest

# N200Mpc at theta_c = 16.8 degs
N200MpcOnaxis_smallest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='3')]['N200MpcOnaxis'].values[0]
# N200Mpc at theta_c = 10.6 degs
N200MpcOnaxis_mean = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='2')]['N200MpcOnaxis'].values[0]
# N200Mpc at theta_c = 7.4 degs
N200MpcOnaxis_biggest = df[(df['label_LF']==label_LF) & (df['label_theta_c']=='1')]['N200MpcOnaxis'].values[0]
# N200Mpc error for the upper bound 
N200MpcOnaxis_up = N200MpcOnaxis_biggest - N200MpcOnaxis_mean
# N200Mpc error for the lower bound 
N200MpcOnaxis_low = N200MpcOnaxis_mean - N200MpcOnaxis_smallest


# the list of full angle for the coincidence
COnaxis_NSNS_list_O3_mean = np.array([N120MpcOnaxis_mean, N120MpcOnaxis_mean, N65MpcOnaxis_mean]) # for coincidence NSNS 
COnaxis_NSBH_list_O3_mean = np.array([N190MpcOnaxis_mean, N190MpcOnaxis_mean, N105MpcOnaxis_mean]) # for coincidence NSBH 
COnaxis_NSNS_list_O4_mean = np.array([N190MpcOnaxis_mean, N190MpcOnaxis_mean, N65MpcOnaxis_mean]) # for coincidence NSNS 
COnaxis_NSBH_list_O4_mean = np.array([N300MpcOnaxis_mean, N300MpcOnaxis_mean, N105MpcOnaxis_mean]) # for coincidence NSBH 
COnaxis_NSNS_list_O4_K_mean = np.array([N190MpcOnaxis_mean, N190MpcOnaxis_mean, N65MpcOnaxis_mean, N40MpcOnaxis_mean]) # for coincidence NSNS 
COnaxis_NSBH_list_O4_K_mean = np.array([N300MpcOnaxis_mean, N300MpcOnaxis_mean, N105MpcOnaxis_mean, N65MpcOnaxis_mean]) # for coincidence NSBH 
COnaxis_NSNS_list_design_mean = np.array([N190MpcOnaxis_mean, N190MpcOnaxis_mean, N125MpcOnaxis_mean, N140MpcOnaxis_mean]) # for coincidence NSNS 
COnaxis_NSBH_list_design_mean = np.array([N300MpcOnaxis_mean, N300MpcOnaxis_mean, N200MpcOnaxis_mean, N225MpcOnaxis_mean]) # for coincidence NSBH 


# the list of the errors for the upper bound 
COnaxis_NSNS_list_O3_up = np.array([N120MpcOnaxis_up, N120MpcOnaxis_up, N65MpcOnaxis_up]) # for coincidence NSNS 
COnaxis_NSBH_list_O3_up = np.array([N190MpcOnaxis_up, N190MpcOnaxis_up, N105MpcOnaxis_up]) # for coincidence NSBH 
COnaxis_NSNS_list_O4_up = np.array([N190MpcOnaxis_up, N190MpcOnaxis_up, N65MpcOnaxis_up]) # for coincidence NSNS 
COnaxis_NSBH_list_O4_up = np.array([N300MpcOnaxis_up, N300MpcOnaxis_up, N105MpcOnaxis_up]) # for coincidence NSBH 
COnaxis_NSNS_list_O4_K_up = np.array([N190MpcOnaxis_up, N190MpcOnaxis_up, N65MpcOnaxis_up, N40MpcOnaxis_up]) # for coincidence NSNS 
COnaxis_NSBH_list_O4_K_up = np.array([N300MpcOnaxis_up, N300MpcOnaxis_up, N105MpcOnaxis_up, N65MpcOnaxis_up]) # for coincidence NSBH 
COnaxis_NSNS_list_design_up = np.array([N190MpcOnaxis_up, N190MpcOnaxis_up, N125MpcOnaxis_up, N140MpcOnaxis_up]) # for coincidence NSNS 
COnaxis_NSBH_list_design_up = np.array([N300MpcOnaxis_up, N300MpcOnaxis_up, N200MpcOnaxis_up, N225MpcOnaxis_up]) # for coincidence NSBH 


# the list of the errors for the lower bound 
COnaxis_NSNS_list_O3_low = np.array([N120MpcOnaxis_low, N120MpcOnaxis_low, N65MpcOnaxis_low]) # for coincidence NSNS 
COnaxis_NSBH_list_O3_low = np.array([N190MpcOnaxis_low, N190MpcOnaxis_low, N105MpcOnaxis_low]) # for coincidence NSBH 
COnaxis_NSNS_list_O4_low = np.array([N190MpcOnaxis_low, N190MpcOnaxis_low, N65MpcOnaxis_low]) # for coincidence NSNS 
COnaxis_NSBH_list_O4_low = np.array([N300MpcOnaxis_low, N300MpcOnaxis_low, N105MpcOnaxis_low]) # for coincidence NSBH 
COnaxis_NSNS_list_O4_K_low = np.array([N190MpcOnaxis_low, N190MpcOnaxis_low, N65MpcOnaxis_low, N40MpcOnaxis_low]) # for coincidence NSNS 
COnaxis_NSBH_list_O4_K_low = np.array([N300MpcOnaxis_low, N300MpcOnaxis_low, N105MpcOnaxis_low, N65MpcOnaxis_low]) # for coincidence NSBH 
COnaxis_NSNS_list_design_low = np.array([N190MpcOnaxis_low, N190MpcOnaxis_low, N125MpcOnaxis_low, N140MpcOnaxis_low]) # for coincidence NSNS 
COnaxis_NSBH_list_design_low = np.array([N300MpcOnaxis_low, N300MpcOnaxis_low, N200MpcOnaxis_low, N225MpcOnaxis_low]) # for coincidence NSBH 




#++++++++++++++++++++++++++++++
# For Fermi only On-axis
#++++++++++++++++++++++++++++++++


print('+'*50)

# the nsns and nsbh coincident (only on-axis) events with Fermi in O2
COnaxis_NSNS_NSBH_O2_mean = (117./365.)*0.8*0.75/(f_R*f_FOV*0.8*T)*(r_nsns*N88MpcOnaxis_mean + r_nsbh*N140MpcOnaxis_mean)
COnaxis_NSNS_NSBH_O2_up = (117./365.)*0.8*0.75/(f_R*f_FOV*0.8*T)*((r_nsns*N88MpcOnaxis_up)**2 + (r_nsbh*N140MpcOnaxis_up)**2)**0.5
COnaxis_NSNS_NSBH_O2_low = (117./365.)*0.8*0.75/(f_R*f_FOV*0.8*T)*((r_nsns*N88MpcOnaxis_low)**2 + (r_nsbh*N140MpcOnaxis_low)**2)**0.5
print('The coincident (only on-axis) events with LH and Fermi in O2 = %.3f - %.3f - %.3f' % (COnaxis_NSNS_NSBH_O2_mean - COnaxis_NSNS_NSBH_O2_low, COnaxis_NSNS_NSBH_O2_mean, COnaxis_NSNS_NSBH_O2_mean + COnaxis_NSNS_NSBH_O2_up))

# the nsns and nsbh number of coincident (only on-axis) events in O3
# the number of NSNS events in O3, seen at least two detectors in three-detector network with Fermi
COnaxis_NSNS_mean_O3, COnaxis_NSNS_up_O3, COnaxis_NSNS_low_O3 = C_i(COnaxis_NSNS_list_O3_mean, COnaxis_NSNS_list_O3_up,  COnaxis_NSNS_list_O3_low, 3, r_nsns, 0.8, 0.75) # for NSNS in O3 with Fermi 
# the number of NSBH events in O3
COnaxis_NSBH_mean_O3, COnaxis_NSBH_up_O3, COnaxis_NSBH_low_O3 = C_i(COnaxis_NSBH_list_O3_mean, COnaxis_NSBH_list_O3_up, COnaxis_NSBH_list_O3_low, 3, r_nsbh, 0.8, 0.75) # for NSNS in O3 with Fermi  
# the number of both NSNS and NSBH events
COnaxis_nsns_nsbh_g_mean_O3  = COnaxis_NSNS_mean_O3  + COnaxis_NSBH_mean_O3 
COnaxis_nsns_nsbh_g_up_O3 = (COnaxis_NSNS_up_O3**2 + COnaxis_NSBH_up_O3**2)**0.5 
COnaxis_nsns_nsbh_g_low_O3 = (COnaxis_NSNS_low_O3**2 + COnaxis_NSBH_low_O3**2)**0.5 

print('='*50)
print('The coincident (only on-axis) events with Fermi in O3 = %f - %f - %f' % (COnaxis_nsns_nsbh_g_mean_O3 - COnaxis_nsns_nsbh_g_low_O3, COnaxis_nsns_nsbh_g_mean_O3, COnaxis_nsns_nsbh_g_mean_O3 + COnaxis_nsns_nsbh_g_up_O3))



# O4 
# the number of NSNS events in O4, seen at least two detectors in three-detector network with Fermi
COnaxis_NSNS_mean_O4, COnaxis_NSNS_up_O4, COnaxis_NSNS_low_O4 = C_i(COnaxis_NSNS_list_O4_mean, COnaxis_NSNS_list_O4_up, COnaxis_NSNS_list_O4_low,3, r_nsns, 0.8, 0.75) # for NSNS in O4 with Fermi 
# the number of NSBH events in O4
COnaxis_NSBH_mean_O4, COnaxis_NSBH_up_O4, COnaxis_NSBH_low_O4 = C_i(COnaxis_NSBH_list_O4_mean, COnaxis_NSBH_list_O4_up, COnaxis_NSBH_list_O4_low, 3, r_nsbh, 0.8, 0.75) # for NSNS in O4 with Fermi 
# the number of both NSNS and NSBH events
COnaxis_nsns_nsbh_g_mean_O4  = COnaxis_NSNS_mean_O4  + COnaxis_NSBH_mean_O4 
COnaxis_nsns_nsbh_g_up_O4 = (COnaxis_NSNS_up_O4**2 + COnaxis_NSBH_up_O4**2)**0.5 
COnaxis_nsns_nsbh_g_low_O4 = (COnaxis_NSNS_low_O4**2 + COnaxis_NSBH_low_O4**2)**0.5 

print('='*50)
print('The coincident (only on-axis) events with Fermi in O4 = %f - %f - %f' % (COnaxis_nsns_nsbh_g_mean_O4 - COnaxis_nsns_nsbh_g_low_O4, COnaxis_nsns_nsbh_g_mean_O4, COnaxis_nsns_nsbh_g_mean_O4 + COnaxis_nsns_nsbh_g_up_O4))


# O4 icluding KAGRA
# the number of NSNS events in O4 including K , seen at least two detectors in three-detector network with Fermi
COnaxis_NSNS_mean_O4_K, COnaxis_NSNS_up_O4_K, COnaxis_NSNS_low_O4_K = C_i(COnaxis_NSNS_list_O4_K_mean, COnaxis_NSNS_list_O4_K_up, COnaxis_NSNS_list_O4_K_low, 4, r_nsns, 0.8, 0.75) # for NSBH in O4 including K with Fermi
 # the number of NSBH events in O4 including K
COnaxis_NSBH_mean_O4_K, COnaxis_NSBH_up_O4_K, COnaxis_NSBH_low_O4_K = C_i(COnaxis_NSBH_list_O4_K_mean, COnaxis_NSBH_list_O4_K_up, COnaxis_NSBH_list_O4_K_low, 4, r_nsbh, 0.8, 0.75) # for NSNS in O4 including K with Fermi 
# the number of both NSNS and NSBH events
COnaxis_nsns_nsbh_g_mean_O4_K  = COnaxis_NSNS_mean_O4_K  + COnaxis_NSBH_mean_O4_K 
COnaxis_nsns_nsbh_g_up_O4_K = (COnaxis_NSNS_up_O4_K**2 + COnaxis_NSBH_up_O4_K**2)**0.5 
COnaxis_nsns_nsbh_g_low_O4_K = (COnaxis_NSNS_low_O4_K**2 + COnaxis_NSBH_low_O4_K**2)**0.5 

print('='*50)
print('The coincident (only on-axis) events with Fermi in O4 including K = %f - %f - %f' % (COnaxis_nsns_nsbh_g_mean_O4_K - COnaxis_nsns_nsbh_g_low_O4_K, COnaxis_nsns_nsbh_g_mean_O4_K, COnaxis_nsns_nsbh_g_mean_O4_K + COnaxis_nsns_nsbh_g_up_O4_K))




# design icluding KAGRA
# the number of NSNS events in O4 including K , seen at least two detectors in three-detector network with Fermi
COnaxis_NSNS_mean_design, COnaxis_NSNS_up_design, COnaxis_NSNS_low_design = C_i(COnaxis_NSNS_list_design_mean, COnaxis_NSNS_list_design_up, COnaxis_NSNS_list_design_low, 4, r_nsns, 0.8, 0.75) # for NSNS in design including K with Fermi
 # the number of NSBH events in O4 including K
COnaxis_NSBH_mean_design, COnaxis_NSBH_up_design, COnaxis_NSBH_low_design = C_i(COnaxis_NSBH_list_design_mean, COnaxis_NSBH_list_design_up, COnaxis_NSBH_list_design_low, 4, r_nsbh, 0.8, 0.75) # for NSNS in design including K with Fermi 
# the number of both NSNS and NSBH events
COnaxis_nsns_nsbh_g_mean_design = COnaxis_NSNS_mean_design  + COnaxis_NSBH_mean_design 
COnaxis_nsns_nsbh_g_up_design = (COnaxis_NSNS_up_design**2 + COnaxis_NSBH_up_design**2)**0.5 
COnaxis_nsns_nsbh_g_low_design = (COnaxis_NSNS_low_design**2 + COnaxis_NSBH_low_design**2)**0.5 

print('='*50)
print('The coincident (only on-axis) events with Fermi in design = %f - %f - %f' % (COnaxis_nsns_nsbh_g_mean_design - COnaxis_nsns_nsbh_g_low_design, COnaxis_nsns_nsbh_g_mean_design, COnaxis_nsns_nsbh_g_mean_design + COnaxis_nsns_nsbh_g_up_design))


