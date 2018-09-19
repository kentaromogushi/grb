import numpy as np
from scipy import integrate
from math import exp, pi, log, log10, pow, sqrt, gamma
import sys
import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
from scipy.special import gammainc, gdtrc
import scipy.optimize

#+++++++++++++++++++++++++++++
# input parameters           |
#+++++++++++++++++++++++++++++
result_file_name = sys.argv[1]
label = sys.argv[2]
file_id = sys.argv[3]
Rate_function_name = sys.argv[4] # unity ~ wilkins, delayed_unity ~  delayed_wilkins
log10L_o_GP = float(sys.argv[5])
a_GP = float(sys.argv[6])
b_GP = float(sys.argv[7])
log10delta1 = float(sys.argv[8])
log10delta2 = float(sys.argv[9])
omega_M = float(sys.argv[10]) # omega matter
omega_A = float(sys.argv[11]) # omega dark energy
t_min = float(sys.argv[12])  # input as Gyr
power_of_timeprb = float(sys.argv[13])
#bin_num_red_withinMinMax = int(sys.argv[14])
binwidth = round(float(sys.argv[14]), 1)
theta_obsGW = float(sys.argv[15])           # import the observation angle for GW170817
Fchosen_str = str(sys.argv[16])             # I expect an imput as 'F_pow', 'F_gau' or 'F_cos'
if Fchosen_str == 'F_pow': # choose a beam profile
    F_c = 1. # This parameter does not affect a rate at all
    theta_c = float(sys.argv[17]) # cut-off angle [radian]
elif Fchosen_str == 'F_gau':
    F_c = 1. # does not affect a rate at all
elif Fchosen_str == 'F_cos':
    F_c = 1.
log10LisoGW = 46.4615 # Isotropic equivalent luminosity for GRB170817
log10typiL_iso = 50.36 # On-axis isotroic equivalent luminosity for a typical GRB
#+++++++++++++++++++++++++++++
#  astronomical parameters   |
#+++++++++++++++++++++++++++++
c = 3.0*10**8   # [m/s]
H_o = 67.89/(3.0857*10**19) # [km/s/Mpc]*[Mpc/km]= [1/s]

#F_swift = 5*pow(10, -9)  # [erg/s/cm^2] Cao et al. 2011
#F_swift = 5.0*pow(10.0, -5.0)  # [erg/s/m^2]
# Swfit minimum flux threshold reffered from Swfit technical notebook 2017
# for 5 sigma and 1 second 
F_swift = 2.8*10.0**(-4.0)  # [erg/s/m^2]

#omega_M = 0.27
#omega_A = 0.73
#omega_M = 0.308
#omega_A = 0.692
#
##lognormal
#L_o_lm = 1.9*10**50 #[erg/s]
#sigma = 2.46153846154
##L_o_lm = 2.15315315315*10**50 #[erg/s]
##sigma = 2.70570570571
#
##shecter
#L_o_sht = 1.01*10**51
#1.06306306306*10**51
#a_sht = 0.65758
#0.544444444444
#
##GP
#L_o_GP = 4.09090909091*10**48
#a_GP = -0.686868686869
#b_GP = 0.909090909091
#
##cutoff time
#t_min = 1.0*pow(10, 9)*365*24*60*60 #[s]
#power_of_timeprb = 1.0 # greater than 1
#

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# making cumulative hist data about logarism of luminosity and redshift  |
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def make_luminohist_zhist():
    
    
    file_id_str = str(file_id)
    
    
    raw_data_matrix = np.loadtxt(file_id_str, dtype = float)
    
    
    duration_vector = raw_data_matrix[:,0]  #[s]
    fluence_vector = raw_data_matrix[:,1] # [10^(-7)erg/cm^2]
    redshift_vector = raw_data_matrix[:,2]
    
    
    
    log10_luminosity_vector = [0 for row in range(len(redshift_vector))]
    
    
    
    #+++++++++++++++++++++++++++++++++
    #                 I(z)           |
    #+++++++++++++++++++++++++++++++++
    
    
    def I(z):
        value, a = integrate.quad(lambda t: 1.0/pow(omega_M*pow(1.0+t, 3.0)+omega_A, 0.5), 0.0, z)
        return value
    
    
    #+++++++++++++++++++++++++++++++
    #              d_L             |
    #+++++++++++++++++++++++++++++++
    
    def d_L(z):
        #c = 3.0*10**8   # [m/s]
        #H_o = 23*pow(10, -19)   # [1/s]
        d_L_value = c*(z+1)*I(z)/H_o    #[m]
        return d_L_value #[m]
    
    
    
    
    for index in range(len(redshift_vector)):
        log10_luminosity_vector[index] += log10(4*pi*pow(d_L(redshift_vector[index]), 2)*fluence_vector[index]*pow(10, -3)/duration_vector[index])
    
    
    #++++++++++++++++++++++++++++++++
    #   redshift                    |
    #++++++++++++++++++++++++++++++++
    startpoint = np.amin(redshift_vector)
    endpoint = 5.0
    
    #    binwidth = (sorted(redshift_vector, reverse=True)[1] - np.amin(redshift_vector))/float(bin_num_red_withinMinMax)
    #    bin_num_red = int(round((endpoint - startpoint)/binwidth))
    
    
    (hist_redshift_vector, binedge_redshift_vector) = np.histogram(redshift_vector, bins=np.arange(0., 5., binwidth), density=False)
    
    binedge_redshift_vector = [round(i, 2) for i in binedge_redshift_vector]
    
    center_of_bin_redshift_vector = [0 for row in range(len(binedge_redshift_vector)-1)]
    for index in range(len(binedge_redshift_vector)-1):
        center_of_bin_redshift_vector[index] += (binedge_redshift_vector[index] + binedge_redshift_vector[index+1])/2
    
    center_of_bin_redshift_vector = [round(i, 2) for i in center_of_bin_redshift_vector]
    cum_hist_redshift_vector = [0 for row in range(len(hist_redshift_vector))]
    for index in range(len(hist_redshift_vector)):
        for sum_num in range(0, index+1):
            cum_hist_redshift_vector[index] += hist_redshift_vector[sum_num]



    return center_of_bin_redshift_vector, cum_hist_redshift_vector, redshift_vector

#++++++++++++++++++++++++++++++++++++++++++
#  launching the function to get values   |
#++++++++++++++++++++++++++++++++++++++++++
(center_of_bin_redshift_vector, cum_hist_redshift_vector, redshift_vector) = make_luminohist_zhist()


#+++++++++++++++++++++++++++++++++++++++
# making cumulative hist for plot      |
#+++++++++++++++++++++++++++++++++++++++
def make_lumino_cumhist_value():
    loglumino = []
    for index, ceter_of_bin in enumerate(center_of_bin_1og10lumino_vector, 0):
        for number in range(0, cum_hist_log10lumino_vector[index]):
            loglumino.append(ceter_of_bin)
    return loglumino



def make_redshift_cumhist_value():
    redshift = []
    for index, ceter_of_bin in enumerate(center_of_bin_redshift_vector, 0):
        for number in range(0, cum_hist_redshift_vector[index]):
            redshift.append(ceter_of_bin)
    return redshift






#+++++++++++++++++++++++++++++++++
#                 I(z)           |
#+++++++++++++++++++++++++++++++++
def I(z):
    #omega_M = 0.27
    #omega_A = 0.73
    value, a = integrate.quad(lambda t: 1.0/pow(omega_M*pow(1.0+t, 3.0)+omega_A, 0.5), 0.0, z)
    return value


#+++++++++++++++++++++++++++++++
#              d_L             |
#+++++++++++++++++++++++++++++++

def d_L(z):
    #c = 3.0*10**8   # [m/s]
    #H_o = 23*pow(10, -19)   # [1/s]
    d_L_value = c*(z+1)*I(z)/H_o    #[m]
    return d_L_value #[m]


#+++++++++++++++++++++++++++++++
#              L_min           |
#+++++++++++++++++++++++++++++++
def L_min(z):
    #c = 3.0*10**8   # [m/s]
    #H_o = 71/(3.0857*10**19)   # [1/s]
    d_L = c*(z+1)*I(z)/H_o    #[m]
    return F_swift*4*pi*pow(d_L, 2.0)   #[erg/s]




#++++++++++++++++++++++++++++++++
#        Rate functions         |
#++++++++++++++++++++++++++++++++
def Rate_unity(z):
    return 1.0


def Rate_porciani(z):
    if z < 100.:
        value = 23.0*exp(3.4*z)/(exp(3.4*z)+22.0) # [dimensionless]
    else:
        value = 23.
    return value




def Rate_Hernquist(z):
    if z < 300.:
        def H(z):
            return H_o*pow(pow(1+z, 3)*omega_M+omega_A, 0.5)
        def chi(z):
            return pow(H(z)/H_o, 2.0/3.0)
        alpha = 0.012
        beta = 0.041
        value = pow(chi(z), 2.0)/(1+alpha*pow(chi(z)-1.0, 3.0)*exp(beta*pow(chi(z), 7.0/4.0)))
    else:
        value = 0.
    return value


def Rate_fardal(z):
    p_1 = 0.075
    p_2 = 3.7
    p_3 = 0.84
    a = pow(1+z, -1.0)
    def H_DL(z):
        return pow(pow(1+z, 3)*omega_M+omega_A, 0.5)
    return pow(a, -p_2)*H_DL(z)/(pow(1.0+p_1*a**(-p_2), p_3+1.0))/0.875403796547




def Rate_cole(z):
    a = 0.0166
    b = 0.1848
    c = 1.9474
    d = 2.6316
    def H_DL(z):
        return pow(pow(1+z, 3.0)*omega_M+omega_A, 0.50)
    return (a+b*z)/(1+pow(z/c, d))*H_DL(z)/0.0166




def Rate_hopkins(z):
    a = 0.017
    b = 0.13
    c = 3.3
    d = 5.3
    def H_DL(z):
        return pow(pow(1+z, 3)*omega_M+omega_A, 0.50)
    return (a+b*z)/(1.0+pow(z/c, d))*H_DL(z)/0.017




def Rate_wilkins(z):
    a = 0.014
    b = 0.11
    c = 1.4
    d = 2.2
    def H_DL(z):
        return pow(pow(1+z, 3.0)*omega_M+omega_A, 0.50)
    return (a+b*z)/(1.0+pow(z/c, d))*H_DL(z)/0.014



#++++++++++++++++++++++++++++++++++++++++++++
#    delayed rate functions                  |
#++++++++++++++++++++++++++++++++++++++++++++
def T_DL(z):  # lookback time
    value, a = integrate.quad(lambda s: 1.0/(1.0+s)*pow(omega_M*pow(1.0+s, 3.0)+omega_A, -0.5), 0.0, z)  # s is reshift, output is dimensionless time
    return value




def Z(s): # analitically
    T_o = 1/H_o  #[s]
    def E(s):
        return exp(log((1.0+sqrt(omega_A))/(1.0-sqrt(omega_A)))-3.0*sqrt(omega_A)*s)
    return pow(omega_A/omega_M*4.0*E(s)/(pow(1.0-E(s), 2.0)), 1/3.0)-1.0  # s is dimensionless time


#def delayed_ratefunction_norm(rate_functoin):   #dimensionless gving the normalization factor
#    T_o = 1.0/H_o  #[s]
#    t_max = T_o - T(0.0)
#    #t_min = 1*pow(10, 9)*365*24*60*60
#    t_maxDL = t_max/T_o
#    t_minDL = t_min/T_o
#    #power_of_timeprb = 1 # greater than 1
#    if  t_maxDL >= t_minDL:
#        norm, a = integrate.quad(lambda s: 1.0/(1+Z(T(0.0)+s))*rate_function(Z(T(0)+s))*s**(-power_of_timeprb), t_minDL, t_maxDL)  #s is time
#    else:
#        norm = 0.0
#    return norm


#############  delayed unity ##################
def delayed_Rate_unity(z):#dimensionless
    T_o = 1.0/(H_o*31536000000000000)  #[Gyr]
    def T_DL(z):  # lookback time
        value, a = integrate.quad(lambda s: 1.0/(1.0+s)*pow(omega_M*pow(1.0+s, 3.0)+omega_A, -0.5), 0.0, z)  # s is reshift, output is dimensionless time
        return value
    
    def Z(s): # analitically
        T_o = 1/H_o  #[s]
        def E(s):
            return exp(log((1.0+sqrt(omega_A))/(1.0-sqrt(omega_A)))-3.0*sqrt(omega_A)*s)
        return pow(omega_A/omega_M*4.0*E(s)/(pow(1.0-E(s), 2.0)), 1/3.0)-1.0  # s is dimensionless time
    
    def delayed_ratefunction_norm():   #dimensionless gving the normalization factor
        t_maxDL = T_DL(np.inf) - T_DL(0.0)
        t_minDL = t_min/T_o   #t_min[Gyr]/T_o[Gyr]
        norm, a = integrate.quad(lambda s: 1.0/(1+Z(T_DL(0.0)+s))*Rate_unity(Z(T_DL(0.0)+s))*s**(-power_of_timeprb), t_minDL, t_maxDL)  #s is time
        return norm
    
    norm = delayed_ratefunction_norm()
    t_maxDL = T_DL(np.inf) - T_DL(z)
    t_minDL = t_min/T_o   #t_min[Gyr]/T_o[Gyr]
    if  t_maxDL >= t_minDL:
        value, a = integrate.quad(lambda s: 1.0/norm*1.0/(1+Z(T_DL(z)+s))*Rate_unity(Z(T_DL(z)+s))*s**(-power_of_timeprb), t_minDL, t_maxDL)  #s is time
    else:
        value = 0.0
    return value




#################  deleyed porciani########################
def delayed_Rate_porciani(z):#dimensionless
    T_o = 1.0/(H_o*31536000000000000)  #[Gyr]
    def T_DL(z):  # lookback time
        value, a = integrate.quad(lambda s: 1.0/(1.0+s)*pow(omega_M*pow(1.0+s, 3.0)+omega_A, -0.5), 0.0, z)  # s is reshift, output is dimensionless time
        return value
    
    def Z(s): # analitically
        T_o = 1/H_o  #[s]
        def E(s):
            return exp(log((1.0+sqrt(omega_A))/(1.0-sqrt(omega_A)))-3.0*sqrt(omega_A)*s)
        return pow(omega_A/omega_M*4.0*E(s)/(pow(1.0-E(s), 2.0)), 1/3.0)-1.0  # s is dimensionless time
    
    
    
    def delayed_ratefunction_norm():   #dimensionless gving the normalization factor
        t_maxDL = T_DL(np.inf) - T_DL(0.0)
        t_minDL = t_min/T_o   #t_min[Gyr]/T_o[Gyr]
        norm, a = integrate.quad(lambda s: 1.0/(1+Z(T_DL(0.0)+s))*Rate_porciani(Z(T_DL(0.0)+s))*s**(-power_of_timeprb), t_minDL, t_maxDL)  #s is time
        return norm
    
    norm = delayed_ratefunction_norm()
    
    
    t_maxDL = T_DL(np.inf) - T_DL(z)
    t_minDL = t_min/T_o   #t_min[Gyr]/T_o[Gyr]
    
    if  t_maxDL >= t_minDL:
        value, a = integrate.quad(lambda s: 1.0/norm*1.0/(1+Z(T_DL(z)+s))*Rate_porciani(Z(T_DL(z)+s))*s**(-power_of_timeprb), t_minDL, t_maxDL)  #s is time
    else:
        value = 0.0
    return value




###########   delayed Hernquist  ###############
def delayed_Rate_Hernquist(z):#dimensionless
    T_o = 1.0/(H_o*31536000000000000)  #[Gyr]
    def T_DL(z):  # lookback time
        value, a = integrate.quad(lambda s: 1.0/(1.0+s)*pow(omega_M*pow(1.0+s, 3.0)+omega_A, -0.5), 0.0, z)  # s is reshift, output is dimensionless time
        return value
    
    def Z(s): # analitically
        T_o = 1/H_o  #[s]
        def E(s):
            return exp(log((1.0+sqrt(omega_A))/(1.0-sqrt(omega_A)))-3.0*sqrt(omega_A)*s)
        return pow(omega_A/omega_M*4.0*E(s)/(pow(1.0-E(s), 2.0)), 1/3.0)-1.0  # s is dimensionless time
    
    def delayed_ratefunction_norm():   #dimensionless gving the normalization factor
        t_maxDL = T_DL(np.inf) - T_DL(0.0)
        t_minDL = t_min/T_o   #t_min[Gyr]/T_o[Gyr]
        norm, a = integrate.quad(lambda s: 1.0/(1+Z(T_DL(0.0)+s))*Rate_Hernquist(Z(T_DL(0.0)+s))*s**(-power_of_timeprb), t_minDL, t_maxDL)  #s is time
        return norm
    
    norm = delayed_ratefunction_norm()
    t_maxDL = T_DL(np.inf) - T_DL(z)
    t_minDL = t_min/T_o   #t_min[Gyr]/T_o[Gyr]
    
    if  t_maxDL >= t_minDL:
        value, a = integrate.quad(lambda s: 1.0/norm*1.0/(1+Z(T_DL(z)+s))*Rate_Hernquist(Z(T_DL(z)+s))*s**(-power_of_timeprb), t_minDL, t_maxDL)  #s is time
    else:
        value = 0.0
    return value




#########  delayed fardal #########################
def delayed_Rate_fardal(z):#dimensionless
    T_o = 1.0/(H_o*31536000000000000)  #[Gyr]
    def T_DL(z):  # lookback time
        value, a = integrate.quad(lambda s: 1.0/(1.0+s)*pow(omega_M*pow(1.0+s, 3.0)+omega_A, -0.5), 0.0, z)  # s is reshift, output is dimensionless time
        return value
    
    def Z(s): # analitically
        T_o = 1/H_o  #[s]
        def E(s):
            return exp(log((1.0+sqrt(omega_A))/(1.0-sqrt(omega_A)))-3.0*sqrt(omega_A)*s)
        return pow(omega_A/omega_M*4.0*E(s)/(pow(1.0-E(s), 2.0)), 1/3.0)-1.0  # s is dimensionless time
    
    def delayed_ratefunction_norm():   #dimensionless gving the normalization factor
        t_maxDL = T_DL(np.inf) - T_DL(0.0)
        t_minDL = t_min/T_o   #t_min[Gyr]/T_o[Gyr]
        norm, a = integrate.quad(lambda s: 1.0/(1+Z(T_DL(0.0)+s))*Rate_fardal(Z(T_DL(0.0)+s))*s**(-power_of_timeprb), t_minDL, t_maxDL)  #s is time
        return norm
    
    norm = delayed_ratefunction_norm()
    t_maxDL = T_DL(np.inf) - T_DL(z)
    t_minDL = t_min/T_o   #t_min[Gyr]/T_o[Gyr]
    
    if  t_maxDL >= t_minDL:
        value, a = integrate.quad(lambda s: 1.0/norm*1.0/(1+Z(T_DL(z)+s))*Rate_fardal(Z(T_DL(z)+s))*s**(-power_of_timeprb), t_minDL, t_maxDL)  #s is time
    else:
        value = 0.0
    return value




################     delayed cole ####################
def delayed_Rate_cole(z):#dimensionless
    T_o = 1.0/(H_o*31536000000000000)  #[Gyr]
    def T_DL(z):  # lookback time
        value, a = integrate.quad(lambda s: 1.0/(1.0+s)*pow(omega_M*pow(1.0+s, 3.0)+omega_A, -0.5), 0.0, z)  # s is reshift, output is dimensionless time
        return value
    
    def Z(s): # analitically
        T_o = 1/H_o  #[s]
        def E(s):
            return exp(log((1.0+sqrt(omega_A))/(1.0-sqrt(omega_A)))-3.0*sqrt(omega_A)*s)
        return pow(omega_A/omega_M*4.0*E(s)/(pow(1.0-E(s), 2.0)), 1/3.0)-1.0  # s is dimensionless time
    
    def delayed_ratefunction_norm():   #dimensionless gving the normalization factor
        t_maxDL = T_DL(np.inf) - T_DL(0.0)
        t_minDL = t_min/T_o   #t_min[Gyr]/T_o[Gyr]
        norm, a = integrate.quad(lambda s: 1.0/(1+Z(T_DL(0.0)+s))*Rate_cole(Z(T_DL(0.0)+s))*s**(-power_of_timeprb), t_minDL, t_maxDL)  #s is time
        return norm
    
    norm = delayed_ratefunction_norm()
    t_maxDL = T_DL(np.inf) - T_DL(z)
    t_minDL = t_min/T_o   #t_min[Gyr]/T_o[Gyr]
    
    if  t_maxDL >= t_minDL:
        value, a = integrate.quad(lambda s: 1.0/norm*1.0/(1+Z(T_DL(z)+s))*Rate_cole(Z(T_DL(z)+s))*s**(-power_of_timeprb), t_minDL, t_maxDL)  #s is time
    else:
        value = 0.0
    return value



################     delayed hopkins ####################
def delayed_Rate_hopkins(z):#dimensionless
    T_o = 1.0/(H_o*31536000000000000)  #[Gyr]
    def T_DL(z):  # lookback time
        value, a = integrate.quad(lambda s: 1.0/(1.0+s)*pow(omega_M*pow(1.0+s, 3.0)+omega_A, -0.5), 0.0, z)  # s is reshift, output is dimensionless time
        return value
    
    def Z(s): # analitically
        T_o = 1/H_o  #[s]
        def E(s):
            return exp(log((1.0+sqrt(omega_A))/(1.0-sqrt(omega_A)))-3.0*sqrt(omega_A)*s)
        return pow(omega_A/omega_M*4.0*E(s)/(pow(1.0-E(s), 2.0)), 1/3.0)-1.0  # s is dimensionless time
    
    def delayed_ratefunction_norm():   #dimensionless gving the normalization factor
        t_maxDL = T_DL(np.inf) - T_DL(0.0)
        t_minDL = t_min/T_o   #t_min[Gyr]/T_o[Gyr]
        norm, a = integrate.quad(lambda s: 1.0/(1+Z(T_DL(0.0)+s))*Rate_hopkins(Z(T_DL(0.0)+s))*s**(-power_of_timeprb), t_minDL, t_maxDL)  #s is time
        return norm
    
    norm = delayed_ratefunction_norm()
    t_maxDL = T_DL(np.inf) - T_DL(z)
    t_minDL = t_min/T_o   #t_min[Gyr]/T_o[Gyr]
    
    if  t_maxDL >= t_minDL:
        value, a = integrate.quad(lambda s: 1.0/norm*1.0/(1+Z(T_DL(z)+s))*Rate_hopkins(Z(T_DL(z)+s))*s**(-power_of_timeprb), t_minDL, t_maxDL)  #s is time
    else:
        value = 0.0
    return value


################     delayed wilkins ####################
def delayed_Rate_wilkins(z):#dimensionless
    T_o = 1.0/(H_o*31536000000000000)  #[Gyr]
    def T_DL(z):  # lookback time
        value, a = integrate.quad(lambda s: 1.0/(1.0+s)*pow(omega_M*pow(1.0+s, 3.0)+omega_A, -0.5), 0.0, z)  # s is reshift, output is dimensionless time
        return value
    
    def Z(s): # analitically
        T_o = 1/H_o  #[s]
        def E(s):
            return exp(log((1.0+sqrt(omega_A))/(1.0-sqrt(omega_A)))-3.0*sqrt(omega_A)*s)
        return pow(omega_A/omega_M*4.0*E(s)/(pow(1.0-E(s), 2.0)), 1/3.0)-1.0  # s is dimensionless time
    
    def delayed_ratefunction_norm():   #dimensionless gving the normalization factor
        t_maxDL = T_DL(np.inf) - T_DL(0.0)
        t_minDL = t_min/T_o   #t_min[Gyr]/T_o[Gyr]
        norm, a = integrate.quad(lambda s: 1.0/(1+Z(T_DL(0.0)+s))*Rate_wilkins(Z(T_DL(0.0)+s))*s**(-power_of_timeprb), t_minDL, t_maxDL)  #s is time
        return norm
    
    norm = delayed_ratefunction_norm()
    t_maxDL = T_DL(np.inf) - T_DL(z)
    t_minDL = t_min/T_o   #t_min[Gyr]/T_o[Gyr]
    
    if  t_maxDL >= t_minDL:
        value, a = integrate.quad(lambda s: 1.0/norm*1.0/(1+Z(T_DL(z)+s))*Rate_wilkins(Z(T_DL(z)+s))*s**(-power_of_timeprb), t_minDL, t_maxDL)  #s is time
    else:
        value = 0.0
    return value






#++++++++++++++++++++++++++++++++++++++++++++
#         co-moving volumm / dz             |
#++++++++++++++++++++++++++++++++++++++++++++


def dV_dz(z):    #[dimensionless]
    dI_dz = 1/pow(omega_M*pow(1.0+z, 3.0)+omega_A, 0.5)  #[1]
    return 4*pi*pow(I(z), 2.0)*dI_dz   # [dimensionless]




#++++++++++++++++++
# Beam profile    |
#++++++++++++++++++


def F_pow(theta_obs, *para):  # uniform and power law decay with 3 parameters
    if 0.<= theta_obs <= np.pi/2.:   #take the symmetricy of the jet into account
        theta_obs = theta_obs
    elif np.pi/2. < theta_obs <= np.pi: #covert the south jet to the north jet
        theta_obs = -theta_obs + np.pi
    F_c = para[0]         # the first paramter about the intensity
    theta_c = para[1]     # the second parameter about the cut-off angle
    s = para[2]           # the third parameter about the index of decay of the off-axis part
    if 0 <= theta_obs <= theta_c:   # descrite the on-axis part (uniform upto some angle)
        output = F_c
    elif theta_c <= theta_obs <= np.pi/2.:    # the off-axis part with decay with a power law
        output = F_c*(theta_obs/theta_c)**(-s)
    return output         # output is the flux for a given observation angle theta

def F_gau(theta_obs, *para):    # Gaussin structure jet with 2 parameters
    if 0.<= theta_obs <= np.pi/2.:   #take the symmetricy of the jet into account
        theta_obs = theta_obs
    elif np.pi/2. < theta_obs <= np.pi: #covert the south jet to the north jet
        theta_obs = -theta_obs + np.pi
    F_c = para[0]           # 1st parameter about the intensity of flux
    sigma = para[1]         # 2nd paramter about how it decays
    return F_c*np.exp(-theta_obs**2/(2*sigma))      # the output is the flux for a given observation angle, theta

def F_cos(theta_obs, *para):    # Flux proportional cos(theta) to the power of 's' as a just preliminary model
    if 0.<= theta_obs <= np.pi/2.:   #take the symmetricy of the jet into account
        theta_obs = theta_obs
    elif np.pi/2. < theta_obs <= np.pi: #covert the south jet to the north jet
        theta_obs = -theta_obs + np.pi
    F_c = para[0]           # 1st parameter about the intensity
    s = para[1]             # 2nd parameter the decay
    return F_c*(np.cos(theta_obs))**(s) # the output is the flux for a given observation angle, theta

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#      Define the detectability for a give z, with Schechter function       |
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculate the constant for a given parameters of luminosity function   |
# Calculate the constant for the total flux                              |
# to improve the computational cost                                      |
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def integrate_GP_rescaled_norm():
    return 1/(-0.5-a_GP)*(1.0-10.**((0.5+a_GP)*log10delta1)) + 1/(-0.5-b_GP)*(10.**((-0.5-b_GP)*log10delta2)-1.0)
def TotFlux(F, *para):  #Find the total flux for a give beam profile with parameters
    totalFlux, er = integrate.quad(lambda theta: np.sin(theta)*F(theta, *para), 0., np.pi/2.) #find the total flux at the distance d_L for the given jet profile
    return totalFlux


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# find the normalization factors for luminosity function|
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if Fchosen_str == 'F_pow': # the uniform and power law decay profile
    Fchosen = F_pow # chose the uniform and power law decay profile function
    def s(log10LisoGW, theta_obsGW, log10Lisoave, theta_c):
        '''
            find the power index with given values
            log10LisoGW = 'Isotropic equivalent luminosity for GRB170817'
            theta_obsGW = 'Observation angle for GRB170817'
            log10Lisoave = 'Average isotropic equivalent luminosity for other GRBs'
            theta_c = 'cut-off angle for this profile'
            This needs theta_obsGW as well as user's chosen theta_c
            '''
        #theta_obsGW = np.pi/180.*theta_obsGW
        #theta_c = np.pi/180.*theta_c
        return (log10Lisoave - log10LisoGW)/np.log10(theta_obsGW/theta_c)
    s = s(log10LisoGW, theta_obsGW, log10typiL_iso, theta_c) # get a power index for structure jet
    totalFlux=TotFlux(Fchosen, F_c, theta_c, s) # total flux over all angle
elif Fchosen_str == 'F_gau': # Gaussian profile
    Fchosen = F_gau # choose Gaussina profile function
    def sigma(log10LisoGW,theta_obsGW,log10Lisoave):
        '''
            find the width of Gaussian profile with given values
            log10LisoGW = 'Isotropic equivalent luminosity for GRB170817'
            theta_obsGW = 'Observation angle for GRB170817'
            log10Lisoave = 'Average isotropic equivalent luminosity for other GRBs'
            This profides an unique value for a given theta_obsGW
            '''
        return theta_obsGW*(np.log10(np.e)/(2*(log10Lisoave-log10LisoGW)))**(0.5) # sigma for Gaussian profile
    sigma = sigma(log10LisoGW, theta_obsGW, log10typiL_iso) # get a sigma for Gaussian profile
    totalFlux=TotFlux(Fchosen, F_c, sigma) # total flux over all angle
elif Fchosen_str == 'F_cos': # cossain profile
    Fchosen = F_cos # choose the cos profile function
    def s_cos(log10LisoGW, theta_obsGW, log10Lisoave):
        '''
            find the width of Gaussian profile with given values
            log10LisoGW = 'Isotropic equivalent luminosity for GRB170817'
            theta_obsGW = 'Observation angle for GRB170817'
            log10Lisoave = 'Average isotropic equivalent luminosity for other GRBs'
            This profides an unique value for a given theta_obsGW
            '''
        return (log10LisoGW - log10Lisoave)/np.log10(np.cos(theta_obsGW)) # power index for cos profile
    s = s_cos(log10LisoGW, theta_obsGW, log10typiL_iso) # get a power index for cos profile, assuming log10Lisoave = log10L_o_shi
    totalFlux=TotFlux(Fchosen, F_c, s) # total flux over all angle


Normforluminofunc = integrate_GP_rescaled_norm()    #Obtain the normalization factor




#+++++++++++++++++++++++++++++++++++++++++++++++++++
#              main fuction                        |
#+++++++++++++++++++++++++++++++++++++++++++++++++++


def integral_GP_logscale_Lsmin(z, F, norm, TF, *para):   #Detectability for a give z
    '''
        local scopes = z, a give profile model, profile parameters, normalization facotor, total flux for a profile model
        variables are z, F and *para
        First, define the minimum source luminosity threshold
        for a given z and observation angle
        Second, find the number of source per unit z and unit observation angle
        at a given z and observaton angle
        Third, find a fraction of detectable sources at a given z
        '''
    def L_s_min(z, theta_obs, *para):
        '''
            First, define the mimimal flux threshold for Swift
            Second, compute the total flux over all the angle
            The third, the detectable minimam source luminoisty
            for a given redshift and observation angle
            '''
        return 4*np.pi*d_L(z)**2*F_swift/F(theta_obs, *para)*TF #the minimum source luminoisty for a given z and the observation angle
    def FofZandtheta_obs(theta_obs, z): #calculate the number of source per unit z per unit observation angle
        '''
            find the number of source per unit z and unit observation angle
            at a given z and observaton angle
            '''
        if z == 0.: #the integratin is 1 as it is normalized
            output = 1.
        else: #calculate the value by integrating with respect of L_s
            def log10x_min(z, theta_obs, *para):
                return log10(L_s_min(z, theta_obs, *para))-(log10L_o_GP+log10(TF/F(0., *para)))
            if log10x_min(z, theta_obs, *para) < -log10delta1:
                output = 1.0
            elif -log10delta1 <= log10x_min(z, theta_obs, *para) < 0.:
                output = 1.0-1.0/norm/(-0.5-a_GP)*(10.**((log10x_min(z, theta_obs, *para)*(-0.5-a_GP)))-10.**(log10delta1*(0.5+a_GP)))
            elif 0. <= log10x_min(z, theta_obs, *para) <  log10delta2:
                output = 1.0 - 1/norm*(1.0/(-0.5-a_GP)*(1.0-10.**(log10delta1*(0.5+a_GP)))+1/(-0.5-b_GP)*(10.**(log10x_min(z, theta_obs, *para)*(-0.5-b_GP))-1.0))
            else:
                output = 0.0
        return output
    #calculate a fraction of the number of source upto z
    # symmetry about the equator gives the same value in theta_obs = [0, pi/2] and [pi/2, pi]
    # so I integrate only the north part, and compensate for a factor of 2 at the end
    Fofz, a = integrate.quad(lambda theta_obs : np.sin(theta_obs)*FofZandtheta_obs(theta_obs, z), 0., np.pi/2.)
    return 2*Fofz/2. # a factor of 2 in the denominator is the normalization factor for the integral of theta_obs




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  the number of detectable sources upto z with a given profile, Fchosen     |
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def N_th(z, N_o, Rate_function, F, *para): # the theoritical number
    value_N, b = integrate.quad(lambda t: N_o*pow(1.0+t, -1.0)*Rate_function(t)*dV_dz(t)*integral_GP_logscale_Lsmin(t, F, Normforluminofunc, totalFlux, *para), 0.0, z)
    return value_N





#+++++++++++++++++++++++++++++++++++
#        plot  rate functions      |
#+++++++++++++++++++++++++++++++++++
#def plot_ratefunctions():
#    x_vector = np.arange(0.0, 16.0, 0.01)
#    y_vector = [0 for row in range(len(x_vector))]
#    Rate_unity_vector = [0 for row in range(len(x_vector))]
#    Rate_porciani_vector = [0 for row in range(len(x_vector))]
#    Rate_Hernquist_vector = [0 for row in range(len(x_vector))]
#    Rate_fardal_vector = [0 for row in range(len(x_vector))]
#    Rate_cole_vector = [0 for row in range(len(x_vector))]
#    Rate_hopkins_vector = [0 for row in range(len(x_vector))]
#    Rate_wilkins_vector = [0 for row in range(len(x_vector))]
#    for index, x in enumerate(x_vector, 0):
#        Rate_unity_vector[index] += Rate_unity(x)
#        Rate_porciani_vector[index] += Rate_porciani(x)
#        Rate_Hernquist_vector[index] += Rate_Hernquist(x)
#        Rate_fardal_vector[index] += Rate_fardal(x)
#        Rate_cole_vector[index] += Rate_cole(x)
#        Rate_hopkins_vector[index] += Rate_hopkins(x)
#        Rate_wilkins_vector[index] += Rate_wilkins(x)
#    
#    
#    plt.semilogy(x_vector, Rate_unity_vector, label='unity')
#    plt.semilogy(x_vector, Rate_porciani_vector, label='porciani')
#    plt.semilogy(x_vector, Rate_Hernquist_vector, label='hernquist')
#    plt.semilogy(x_vector, Rate_fardal_vector, label='fardal')
#    plt.semilogy(x_vector, Rate_cole_vector, label='cole')
#    plt.semilogy(x_vector, Rate_hopkins_vector, label='hopkins')
#    plt.semilogy(x_vector, Rate_wilkins_vector, label='wilkins')
#
#
#
#    plt.ylim([0.5,200])
#    plt.legend(loc = 'lower right')#
#    plt.xlabel('redshift')
#    plt.ylabel('rate [demensionless]')
#    #plt.ylim(0, 100)
#    plt.title('rate functions')
#    plt.grid(True)
#    #plt.savefig("test.png")
#    plt.show()
#    return
#
##plot_ratefunctions()
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#        plot  comoving volume & integral lognormal from L_min     |
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#def plot_comoving():
#    
#    z_vector = np.arange(0.0, 5.0, 0.01)
#    y_vector = [0 for row in range(len(z_vector))]
#    k_vector = [0 for row in range(len(z_vector))]
#    l_vector = [0 for row in range(len(z_vector))]
#    for index, z in enumerate(z_vector, 0):
#        #y_vector[index] += dV_dz(z)
#        #        k_vector[index] += integral_lognormal_from_min(z)
#        l_vector[index] += integrate_GP_rescaled_fromLmin(z)
#    #plt.plot(z_vector, y_vector)
#    #    plt.semilogy(z_vector, k_vector)
#    #plt.plot(z_vector, y_vector)
#    plt.plot(z_vector, l_vector)
#    plt.xlabel('redshift')
#    plt.ylabel('dv/dz (co-moving shell) [demensionless]')
#    plt.title('co-moving shell versus redshift')
#    #    plt.ylabel('number')
#    #    plt.title('cumulative number of lognormal from L_min')
#    plt.grid(True)
#    #plt.savefig("test.png")
#    plt.show()
#    return 0
#
#plot_comoving()
#
#++++++++++++++++++++++++++++++++++++++++++++
#        plot  normalized lognormal         |
#++++++++++++++++++++++++++++++++++++++++++++
#def plot1():
#    L_vector = np.arange(40, 55, 0.1)
#    x_vector = [0 for row in range(len(L_vector))]
#    a_vector = [0 for row in range(len(L_vector))]
#    b_vector = [0 for row in range(len(L_vector))]
#    L_o = 1.9*10**50
#    for index, L in enumerate(L_vector, 0):
#        
#        x_vector[index] += pow(10, L)/L_o
#        #a_vector[index] += lognormal_normalized(x_vector[index])
#        #a_vector[index] += normalized_lognormal_testDL(x_vector[index])
#        #b_vector[index] += lognormal_normalized_rescaled(x_vector[index])
#        b_vector[index] += normalized_lognormal_rescaled_testDL(x_vector[index])
#    integral_lognormal_from_minDL
#    
#    #y_vector[index] += schechter_normalized(L_vector[index])
#    y_vector[index] += GP_normalized(L_vector[index])
#    
#    
#    plt.plot(L_vector, a_vector)
#    plt.plot(L_vector, b_vector)
#    
#    plt.xlabel('logL [log(erg)]')
#    plt.ylabel('number [1/erg]')
#    plt.title('luminosity distribution')
#    plt.grid(True)
#    plt.savefig("test.png")
#    plt.show()
#    return 0
#

#
##+++++++++++++++++++++++++++++++++++++++++++
##        plot N_th versus N_ex             |
##+++++++++++++++++++++++++++++++++++++++++++
#def plot_N(N_o, Rate_function):
#    x_vector = np.arange(0.0, 5, 0.5)
#    y_vector = [0 for row in range(len(x_vector))]
#    l_vector = [0 for row in range(len(x_vector))]
#    for index, x in enumerate(x_vector, 0):
#        y_vector[index] += N_th(x, N_o, Rate_function)
#    #y_vector[index] += L_min(x)
#    #l_vector[index] += L_min(x)/(6*pow(10, 49))
#    plt.plot(x_vector, y_vector)
#    #plt.semilogy(x_vector, y_vector)
#    #plt.semilogy(x_vector, l_vector)
#
#    #file_name = sys.argv[2]
#    #exp_data_vector = np.loadtxt(file_name)
#
#    z_vector = center_of_bin_redshift_vector
#    N_ex_vector = cum_hist_redshift_vector
#    #plt.plot(z_vector, N_ex_vector)
#    plt.hist(make_redshift_cumhist_value(), bins=np.arange(0., 5., binwidth), color='red', alpha=0.5)
#    plt.xlabel('z')
#    plt.ylabel('number')
#    plt.title('N_th versus N_ex ')
#    plt.grid(True)
#    #plt.savefig("test.png")
#    plt.show()
#    return 0




##++++++++++++++++++++++++++++++++++++++
##  plot d_L (luminosity distance      |
##++++++++++++++++++++++++++++++++++++++
#
#def plot_d_L():
#    z_vector = np.arange(0.0, 0.3, 0.01)
#    y_vector = [0 for row in range(len(z_vector))]
#    for index, z in enumerate(z_vector, 0):
#        y_vector[index] += c/(3.086*10**22)
#    plt.plot(z_vector, y_vector)
#    plt.ylim([600,1000])
#    plt.xlabel('z')
#    plt.ylabel('d_L [Mpc]')
#    plt.title('luminosity distance vs redshift')
#    plt.grid(True)
#    #plt.savefig("test.png")
#    plt.show()
#    return 0


#####################################################################
#print (integrate_lognormal())
#print (integrate_lognormal_rescaled())
#print (integral_lognormal_from_minDL(0))
#print (d_L(0.018)/(3.086*10**22))



#corect = correction_N_o(1)
#print (corect)

#plot_ratefunctions()
#plot_comoving()
#plot_N()
#plot1()
#print (N_th(1, ))
#print (dN_dt_dV(0, 139.13963963944462))
#plot_d_L()
#++++++++++++++++++++++++++++++++++++++++++++



if Rate_function_name == "unity":
    Rate_function = (lambda z: Rate_unity(z))
elif Rate_function_name == "hernquist":
    Rate_function = (lambda z: Rate_Hernquist(z))
elif Rate_function_name == "porciani":
    Rate_function = (lambda z: Rate_porciani(z))
elif Rate_function_name == "fardal":
    Rate_function = (lambda z: Rate_fardal(z))
elif Rate_function_name == "cole":
    Rate_function = (lambda z: Rate_cole(z))
elif Rate_function_name == "hopkins":
    Rate_function = (lambda z: Rate_hopkins(z))
elif Rate_function_name == "wilkins":
    Rate_function = (lambda z: Rate_wilkins(z))
elif Rate_function_name == "delayed_unity":
    Rate_function = (lambda z: delayed_Rate_unity(z))
elif Rate_function_name == "delayed_hernquist":
    Rate_function = (lambda z: delayed_Rate_Hernquist(z))
elif Rate_function_name == "delayed_porciani":
    Rate_function = (lambda z: delayed_Rate_porciani(z))
elif Rate_function_name == "delayed_fardal":
    Rate_function = (lambda z: delayed_Rate_fardal(z))
elif Rate_function_name == "delayed_cole":
    Rate_function = (lambda z: delayed_Rate_cole(z))
elif Rate_function_name == "delayed_hopkins":
    Rate_function = (lambda z: delayed_Rate_hopkins(z))
elif Rate_function_name == "delayed_wilkins":
    Rate_function = (lambda z: delayed_Rate_wilkins(z))
else:
    print ("please the correct rate fucntion name")




#+++++++++++++++++++++++++++++++++
#    write out the results       |
#+++++++++++++++++++++++++++++++++
"""
def plot_N_all():
    x_vector = np.arange(0.0, 5, 0.1)
    unity_vector = [0 for row in range(len(x_vector))]
    hernquist_vector = [0 for row in range(len(x_vector))]
    porciani_vector = [0 for row in range(len(x_vector))]
    fardal_vector = [0 for row in range(len(x_vector))]
    cole_vector = [0 for row in range(len(x_vector))]
    hopkins_vector = [0 for row in range(len(x_vector))]
    wilkins_vector = [0 for row in range(len(x_vector))]
    for index, x in enumerate(x_vector, 0):
        unity_vector[index] += N_th(x, 110, (lambda z: Rate_unity(z)))
        hernquist_vector[index] += N_th(x,  76, (lambda z: Rate_Hernquist(z)))
        porciani_vector[index] += N_th(x, 22, (lambda z: Rate_porciani(z)))
        fardal_vector[index] += N_th(x,  30, (lambda z: Rate_fardal(z)))
        cole_vector[index] += N_th(x, 14, (lambda z: Rate_cole(z)))
        hopkins_vector[index] += N_th(x, 16, (lambda z: Rate_hopkins(z)))
        wilkins_vector[index] += N_th(x, 22, (lambda z: Rate_wilkins(z)))
    
    #######delayed
    #        unity_vector[index] += N_th(x, 93, (lambda z: delayed_Rate_unity(z)))
    #        hernquist_vector[index] += N_th(x,  64, (lambda z: delayed_Rate_Hernquist(z)))
    #        porciani_vector[index] += N_th(x, 19, (lambda z: delayed_Rate_porciani(z)))
    #        fardal_vector[index] += N_th(x,  25, (lambda z: delayed_Rate_fardal(z)))
    #        cole_vector[index] += N_th(x, 12, (lambda z: delayed_Rate_cole(z)))
    #        hopkins_vector[index] += N_th(x, 14, (lambda z: delayed_Rate_hopkins(z)))
    #        wilkins_vector[index] += N_th(x, 17, (lambda z: delayed_Rate_wilkins(z)))
    
    plt.plot(x_vector, unity_vector, label = "unity")
    plt.plot(x_vector, hernquist_vector, label = "hernquist")
    plt.plot(x_vector, porciani_vector, label = "porciani")
    plt.plot(x_vector, fardal_vector, label = "fardal")
    plt.plot(x_vector, cole_vector, label = "cole")
    plt.plot(x_vector, hopkins_vector, label = "hopkins")
    plt.plot(x_vector, wilkins_vector, label = "wilkins")



    z_vector = center_of_bin_redshift_vector
    N_ex_vector = cum_hist_redshift_vector
    #plt.plot(z_vector, N_ex_vector, label = "observed")
    plt.hist(make_redshift_cumhist_value(), bins=np.arange(0., 5., binwidth), color='red', alpha=0.5)
    plt.legend(loc='lower right')
    plt.xlabel('z')
    plt.ylabel('number')
    plt.title('N_th versus N_ex (schechter delay)')
    plt.grid(True)
    plt.savefig("test.png")
    plt.show()
    return 0

plot_N_all()  # plot all the parturns

"""

#+++++++++++++++++++++++++++++++++++++++++++++
#    find N_o by least-sqaure method fitting |
#+++++++++++++++++++++++++++++++++++++++++++++

def find_fitted_parameters_by_curve_fit():
    
    
    def N_th1(z, N_o):
        return N_th(z, N_o, Rate_function)
    
    
    parameter_initial = np.array([1.])
    
    parabounds = (400., 800.)
    
    
    
    paramater_optimal, covariance = scipy.optimize.curve_fit(np.vectorize(N_th1), center_of_bin_redshift_vector, cum_hist_redshift_vector, p0=parameter_initial)
    
    std_dv = np.sqrt(np.diag(covariance))
    return paramater_optimal, std_dv



if Fchosen_str == 'F_pow': # Fchosen will be F_pow
    N_o = len(redshift_vector)/N_th(5.0, 1.0, Rate_function, Fchosen, F_c, theta_c, s) # find N_o
    N80Mpc = N_th(0.018, 1., Rate_function, Fchosen, F_c, theta_c, s) # N(z=0.022)/N_o
    N120Mpc = N_th(0.027, 1., Rate_function, Fchosen, F_c, theta_c, s) # N(z=0.031)/N_o
elif Fchosen_str == 'F_gau':
    N_o = len(redshift_vector)/N_th(5.0, 1.0, Rate_function, Fchosen, F_c, sigma) # Fchosen is F_gau
    N80Mpc = N_th(0.018, 1., Rate_function, Fchosen, F_c, sigma) #N(z=0.022)/N_o
    N120Mpc = N_th(0.027, 1., Rate_function, Fchosen, F_c, sigma) #N(z=0.031)/N_o
elif Fchosen_str == 'F_cos':
    N_o = len(redshift_vector)/N_th(5.0, 1.0, Rate_function, Fchosen, F_c, s) # Fchosen is F_cos
    N80Mpc = N_th(0.018, 1., Rate_function, Fchosen, F_c, s) ##N(z=0.022)/N_o
    N120Mpc = N_th(0.027, 1., Rate_function, Fchosen, F_c, s) ##N(z=0.031)/N_o




#+++++++++++++++++++++++++++++++++++++++++++++
#  chi-square test                           |
#+++++++++++++++++++++++++++++++++++++++++++++






def chi_sqr_test(fun, x_vector, y_vector):  # step only
    for index, y in enumerate(y_vector, 0):
        if y == len(redshift_vector):
            upperlimit_index = index
            break
    x_step_vector = x_vector[:upperlimit_index+1]
    y_step_vector = y_vector[:upperlimit_index+1]
    N_th_vector = []
    non_zero_comp_num = 0
    for index, x in enumerate(x_step_vector, 0):
        if fun(x, N_o, Rate_function) > 0.0:
            non_zero_comp_num += 1
            N_th_vector.append(fun(x, N_o, Rate_function))
    # the degree of freedom is (the # of user preferd step) + (# 2 parts of flat) - 1 - (# of parameters), which is equivalent to (# of non zero data points) - 1 - (# of constraint)
    print (x_step_vector)
    print (y_step_vector)
    from scipy.stats import chisquare
    chisq, pvalue = chisquare(y_step_vector, np.array(N_th_vector), ddof=1)
    reduced_chisq = chisq/(len(y_step_vector)-2)
    return chisq, reduced_chisq, pvalue


#++++++++++++++++++++++++++++++++++++++++++++++
# write out the results to a file             |
# this is used for putting constraint on      |
# a set of the beam paramters                 |
#++++++++++++++++++++++++++++++++++++++++++++++



with open(result_file_name, mode = 'a') as fh:
    label_str = str(label)
    N_o_str = str('%03.6f' % N_o)
    N80Mpc_str = str('%03.6E' %N80Mpc)
    N120Mpc_str = str('%03.6E' %N120Mpc)
    if Fchosen_str == 'F_pow':
        theta_c_str = str('%.f' % (180./np.pi*theta_c))  # degrees
        s_str = str('%f' %s )
        theta_obsGW_str = str('%.f' % (180./np.pi*theta_obsGW))  # degrees
        fh.write(label_str + " " + N_o_str + " " + N80Mpc_str + " " + N120Mpc_str + " " + theta_c_str + " " + s_str + " " + theta_obsGW_str)
    elif Fchosen_str == 'F_gau':
        sigma_str = str('%.6E' %sigma)
        fh.write(label_str + " " + N_o_str + " " + N80Mpc_str + " " + N120Mpc_str + " " + sigma_str)
    fh.write('\n')
    fh.close()





