import numpy as np
from scipy import integrate
from math import exp, pi, log, log10, pow, gamma, sqrt, cos, pi, sin
import sys
from scipy.special import gammainc
from scipy.stats import norm
import matplotlib.pyplot as plt
import scipy.optimize
from pynverse import inversefunc
from matplotlib.ticker import MultipleLocator

file_id = sys.argv[1]                   # import data file
binwidth = round(float(sys.argv[2]), 1)  # user choosed parameter about the width of the bins

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# making data series for cumulative histogram about luminosity and redshift  |
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def make_luminohist_zhist():
    '''
    description:
        1. import a metadata comprising SGRB observational data
        2. split the variables into array
        3. calculate the cumulative number of the histogrm
        4. calculate the center values of the bins
    :return:
        center_of_bin_1og10lumino_vector: center values of the bins
        cum_hist_log10lumino_vector: cumulative number
        log10_luminosity_vector: log10 of luminosities of SGRBs
        center_of_bin_redshift_vector: center values of redshift histogram
        cum_hist_redshift_vector, redshift_vector: cumulative number of redshift histogram
    '''
    file_id_str = str(file_id)                                      # import data file which contains a sample
    raw_data_matrix = np.loadtxt(file_id_str, dtype = float)        # taking data matrix from the import data file
    duration_vector = raw_data_matrix[:,0]  #[s]                    # taking a vector about duration of the sample
    fluence_vector = raw_data_matrix[:,1] # [10^(-7)erg/cm^2]       # taking a vector about fluence of the sample
    redshift_vector = raw_data_matrix[:,2]                          # taking a vector about redshift af the sample
    log10_luminosity_vector = [0 for row in range(len(redshift_vector))]
    
    
    
    #+++++++++++++++++++++++++++++++++
    #                 I(z)           |
    #+++++++++++++++++++++++++++++++++
    c = 3.0*10**8   # [m/s]
    H_o = 67.89/(3.0857*10**19) # [km/s/Mpc]*[Mpc/km]= [1/s]
    omega_M = 0.308
    omega_A = 0.692
    
    def I(z):
        '''
        :param z: redshift
        :return: z dependent function for luminosity distance
        '''
        value, a = integrate.quad(lambda t: 1.0/pow(omega_M*pow(1.0+t, 3.0)+omega_A, 0.5), 0.0, z)
        return value
    
    
    #+++++++++++++++++++++++++++++++
    #              d_L             |
    #+++++++++++++++++++++++++++++++
    
    def d_L(z):
        '''
        luminosity distance
        :param z: redshift
        :return: luminosity distance in meter
        '''
        #c = 3.0*10**8   # [m/s]
        #H_o = 23*pow(10, -19)   # [1/s]
        d_L_value = c*(z+1)*I(z)/H_o    #[m]
        return d_L_value #[m]

    # calculate luminosities of SGRBs and fill them in to the list
    for index in range(len(redshift_vector)):
        log10_luminosity_vector[index] += log10(4*pi*pow(d_L(redshift_vector[index]), 2)*fluence_vector[index]*pow(10, -3)/duration_vector[index])
    

    


    #++++++++++++++++++++++++++++++++
    #   luminosity                  |
    #++++++++++++++++++++++++++++++++

    # find the histogram of log10(L) where L is luminosity
    # frequency and the list of edges of the bins
    (hist_log10lumino_vector, binedge_log10lumino_vector) = np.histogram(log10_luminosity_vector, bins=np.arange(-250., 250., binwidth), density=False)        # fill data of luminoisty of the sample into the bins of a histogram
    # round the list to reduce the memory
    binedge_log10lumino_vector = [round(i, 2) for i in binedge_log10lumino_vector]              # taking values of the edge of the bins
    # the empty list for the center value of the bins
    center_of_bin_1og10lumino_vector = [0 for row in range(len(binedge_log10lumino_vector)-1)]
    #  calculate the center values of the bins
    for index in range(len(binedge_log10lumino_vector)-1):
        center_of_bin_1og10lumino_vector[index] += (binedge_log10lumino_vector[index] + binedge_log10lumino_vector[index+1])/2      # taking value at the center of the bins

    center_of_bin_1og10lumino_vector = [round(i, 2) for i in center_of_bin_1og10lumino_vector]          # round them up for computational convenience
    # make the cumulative number
    cum_hist_log10lumino_vector = [0 for row in range(len(hist_log10lumino_vector))]
    for index in range(len(hist_log10lumino_vector)):
        for sum_num in range(0, index+1):
            cum_hist_log10lumino_vector[index] += hist_log10lumino_vector[sum_num]          # making cumulative number of luminosity



    #+++++++++++++++++++++++++++++++++++++
    #++++++++++++++++++++++++++++++++
    #   redshift                    |
    #++++++++++++++++++++++++++++++++
    bin_num_red = len(redshift_vector)
    (hist_redshift_vector, binedge_redshift_vector) = np.histogram(redshift_vector, bins=bin_num_red, density=False)
    
    
    center_of_bin_redshift_vector = [0 for row in range(len(binedge_redshift_vector)-1)]
    for index in range(len(binedge_redshift_vector)-1):
        center_of_bin_redshift_vector[index] += (binedge_redshift_vector[index] + binedge_redshift_vector[index+1])/2


    cum_hist_redshift_vector = [0 for row in range(len(hist_redshift_vector))]
    for index in range(len(hist_redshift_vector)):
        for sum_num in range(0, index+1):
            cum_hist_redshift_vector[index] += hist_redshift_vector[sum_num]

    return center_of_bin_1og10lumino_vector, cum_hist_log10lumino_vector, log10_luminosity_vector, center_of_bin_redshift_vector, cum_hist_redshift_vector, redshift_vector

#++++++++++++++++++++++++++++++++++++++++++
#  launching the function to get values   |
#++++++++++++++++++++++++++++++++++++++++++

# run the function defined above
(center_of_bin_1og10lumino_vector, cum_hist_log10lumino_vector, log10_luminosity_vector, center_of_bin_redshift_vector, cum_hist_redshift_vector, redshift_vector) = make_luminohist_zhist() # get the luminoisty vector and the cumulative number of luminoisty


#+++++++++++++++++++++++++++++++++++++++
# making cumulative hist for plot      |
#+++++++++++++++++++++++++++++++++++++++
def make_lumino_cumhist_value():
    '''
    this is used for making a plot
    :return: the list of cumulative number of histgram of SRGB luminosities
    '''
    loglumino = []
    for index, ceter_of_bin in enumerate(center_of_bin_1og10lumino_vector, 0):
        for number in range(0, cum_hist_log10lumino_vector[index]):
            loglumino.append(ceter_of_bin)
    return loglumino



def make_redshift_cumhist_value():
    '''
    this is used for making a plot
    :return: the list of cumulative number of histgram of SRGB redshifts
    '''
    redshift = []
    for index, ceter_of_bin in enumerate(center_of_bin_redshift_vector, 0):
        for number in range(0, cum_hist_redshift_vector[index]):
            redshift.append(ceter_of_bin)
    return redshift


#+++++++++++++++++++++++++++++++++
# change the name for fitting    |
#+++++++++++++++++++++++++++++++++
log10L_vector = center_of_bin_1og10lumino_vector        # rename it
N_ex_vector =  cum_hist_log10lumino_vector              # rename it



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#     shcechter function with python module of gammainc and gamma   |
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




def N_th_sch3logscale(log10L, C, log10L_o, a, log10delta): # convert L -> log10(L/L_o)
    '''
    calculate by means of log10L
    :param log10L: the list of log10 of L
    :param C: a LF parameter
    :param log10L_o: a LF parameter
    :param a: a LF parameter
    :return: cumulative number of the LF
    '''
    if log10L - log10L_o < -log10delta:
        output = 0.0
    else:
        value, s = integrate.quad(lambda x: log(10)*C*pow(10., (1-a)*x)*exp(-pow(10., x)), -log10delta, log10L - log10L_o)
        output = value
    return output



def N_th_sche3(log10L, C, log10L_o, a, log10delta):
    '''
    :param log10L: the list of log10 of L
    :param C: a LF parameter
    :param log10L_o: a LF parameter
    :param a: a LF parameter
    :return: cumulative number of the LF
    '''
    L = pow(10., log10L)
    L_o = pow(10., log10L_o)
    delta = pow(10., log10delta)
    if L < L_o/delta:
        output = 0.
    else:
        output = C*gamma(1-a)*(gammainc(1-a, L/L_o)-gammainc(1-a, 1.0/delta))
    return output




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#      Plot of  schechter  distribution                                     |
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def plot(C, L_o, a):
    '''
    :param C: a LF parameter
    :param L_o: a LF parameter
    :param a: a LF parameter
    :return: 0
    '''

    logL_vector = np.arange(43, 60.0, 0.1)
    L_vector = [0 for row in range(len(logL_vector))]
    y_vector = [0 for row in range(len(logL_vector))]
    for index, logL in enumerate(logL_vector, 0):
        L_vector[index] += pow(10, logL)
        #y_vector[index] += N_th_GP(L_vector[index], N_o, L_o, a, b)
        y_vector[index] += N_th_sche1(L_vector[index], C, L_o, a)
    plt.plot(logL_vector, y_vector)
    plt.hist(make_lumino_cumhist_value(), bins=bin_num_lumino, range=(startpoint, endpoint), color='red', alpha=0.5)
    plt.xlabel('logL')
    plt.ylabel('number')
    plt.title('schechter')
    plt.grid(True)
    plt.savefig("test.png")
    plt.show()
    return 0


def N_th_sche2(log10L, C, log10L_o, a, log10delta): # free delta
    '''
    the cumulative value of the LF
    :param log10L:
    :param C: a LF parameter
    :param log10L_o: a LF parameter
    :param a: a LF parameter
    :param log10delta: a LF parameter
    :return: cumulative number of the LF
    '''
    L = pow(10., log10L)
    L_o = pow(10., log10L_o)
    delta = pow(10., log10delta)
    if L < L_o/delta:
        output = 0.
    else:
        output = C*gamma(1-a)*(gammainc(1-a, L/L_o)-gammainc(1-a, 1.0/delta))
    return output




# def intN(log10L, C ,log10L_0, theta_0_deg, s):
#     L = 10**log10L
#     L_0 = 10**log10L_0
#     theta_0 = theta_0_deg/180.*pi
#     Lmin = L_0*(40./theta_0_deg)**(-s)
#     sigma = 0.5
#     t = L/L_0
#     if t < 1.:
#         #value1, a = integrate.quad(lambda k: C*1000.*sin(theta_0*(k)**(-1/s))*(k)**(0.5-1/s), Lmin/L_0, L/L_0)
#         value1, a = integrate.quad(lambda k: C*(k)**(0.5-1/s), Lmin/L_0, L/L_0)
#         output = value1
#     if t >= 1.:
#         #value2, a = integrate.quad(lambda k: C*1000.*sin(theta_0*(k)**(-1/s))*(k)**(0.5-1/s), Lmin/L_0, 1.)
#         #value3, a = integrate.quad(lambda x: C*1000.*sin(theta_0)*exp(1.5*x)*exp(-x**2/(2*sigma**2)), 0., log(10.)*(log10L - log10L_0))
#         value2, a = integrate.quad(lambda k: C*(k)**(0.5-1/s), Lmin/L_0, 1.)
#         value3, a = integrate.quad(lambda x: C*exp(1.5*x)*exp(-x**2/(2*sigma**2)), 0., log(10.)*(log10L - log10L_0))
#         output = value2 + value3
#     return output

def find_fitted_parameters_by_curve_fit():
    '''
    fitting the LF against the data
    :return: fitted parameters and the errors
    '''
    paramater_initial = np.array([5., 50., 0.1, 2.])
    param_bounds=([5.0 , 46., 0.01, 0.1],[100.0, 52., 0.99, 7.])
    paramater_optimal, covariance = scipy.optimize.curve_fit(np.vectorize(N_th_sche2), log10L_vector, N_ex_vector, p0=paramater_initial, bounds=param_bounds)

    std_dv = np.sqrt(np.diag(covariance))
    return paramater_optimal, std_dv

# print out the fitted parameters
(parameter_optimal, std_dv) = find_fitted_parameters_by_curve_fit()
print ("paramater (C, L_o, a, delta)=", parameter_optimal[0], parameter_optimal[1], parameter_optimal[2], parameter_optimal[3])
print ("standard deviation (C, L_o, a, delta)=", std_dv)




def plot2(C, log10L_o, a, log10delta):
    '''

    :param C: a LF parameter
    :param log10L_o: a LF parameter
    :param a: a LF parameter
    :param log10delta: a LF parameter
    :return: 0
    '''
    log10L_vector = np.arange(43, 60.0, 0.1)
    y_vector = [0 for row in range(len(log10L_vector))]
    for index, log10L in enumerate(log10L_vector, 0):
        y_vector[index] += N_th_sche2(log10L, C, log10L_o, a, log10delta)
    plt.plot(log10L_vector, y_vector)
    plt.hist(make_lumino_cumhist_value(), bins=np.arange(-250., 250., binwidth), color='red', alpha=0.5)
    plt.xlabel('logL')
    plt.ylabel('number')
    plt.xlim([46, 54])
    plt.title('schechter')
    plt.xticks([46, 48, 50, 52, 54], ['10$^{46}$', '10$^{48}$', '10$^{50}$', '10$^{52}$', '10$^{54}$'])
    #plt.yticks([len(log10_luminosity_vector)*i for i in [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]], [0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
    plt.grid(True)
    plt.savefig("test.png")
    plt.show()
    return 0



#+++++++++++++++++++++++++++++++++++
# find Luminosity of a typical GRB |
#+++++++++++++++++++++++++++++++++++

def FindTypicalGRB(C, log10L_o, a, log10delta):
    '''
        Find the typical on axis L_iso,
        that is defined as the one giving the 50% of the cumulative number,
        i.e., the middle of luminosity
        '''
    # total number of the cumulative number
    totalnum = len(redshift_vector)
    # cumulative number with optimal parameters
    def cumNwithOptV(log10L):
        '''
            This is used for finding the a typical L_iso
            '''
        return N_th_sche2(log10L, C, log10L_o, a, log10delta)
    #inverse function of cumNwithOptV, giving a typical L_iso
    typL_iso = inversefunc(cumNwithOptV, domain=[48.5, 52.])
    return typL_iso(totalnum/2.)




# plot
plot2(parameter_optimal[0], parameter_optimal[1], parameter_optimal[2], parameter_optimal[3])




#+++++++++++++++++++++++++++++++++++++++++++++
#  chi-square test                           |
#+++++++++++++++++++++++++++++++++++++++++++++


def chi_sqr_test(fun, x_vector, y_vector):
    '''
    :param fun: a given function
    :param x_vector:
    :param y_vector:
    :return:
        chisq: the vlaue of chi-square
        non_zero_comp_num: the number non-zero value used to compute chi-square
        pvalue: p-value
    '''
    for index, y in enumerate(y_vector, 0):
        if y > 0.:
            starting_non_zero_index = index
            break
    for index, y in enumerate(y_vector, 0):
        if y == len(log10_luminosity_vector):
            upperlimit_index = index
            break
    N_th_vector = []
    non_zero_comp_num = 0
    x_step_vector = x_vector[starting_non_zero_index-1:upperlimit_index+1]
    y_step_vector = y_vector[starting_non_zero_index-1:upperlimit_index+1]
    #print (x_step_vector)
    #print (y_step_vector)
    for index, x in enumerate(x_step_vector, 0):
        if fun(x, parameter_optimal[0], parameter_optimal[1], parameter_optimal[2], parameter_optimal[3]) > 0.0:
            non_zero_comp_num += 1
            N_th_vector.append(fun(x, parameter_optimal[0], parameter_optimal[1], parameter_optimal[2], parameter_optimal[3]))

    #print (N_th_vector)
    #print (y_step_vector[len(y_step_vector)-len(N_th_vector):])
    #print (4 -(len(y_step_vector) - non_zero_comp_num))
    # the degree of freedom is (the # of user preferd step) + (# 2 parts of flat) - 1 - (# of parameters), which is equivalent to (# of non zero data points) - 1 - (# of constraint)
    
    from scipy.stats import chisquare
    chisq, pvalue = chisquare(y_step_vector[len(y_step_vector)-len(N_th_vector):], np.array(N_th_vector), ddof=4 - (len(y_step_vector)- non_zero_comp_num))
    return chisq, non_zero_comp_num, pvalue


print("the value of the chi square, the number of data points with non-zero value (N_th), p-value")
chisquare, NonzeroNth, p = chi_sqr_test(N_th_sche2, log10L_vector, N_ex_vector)
print(chisquare, NonzeroNth, p)

# find the luminosity of typical GRB
print('Typical on-axis isotropic equivalent log10L is %.3f' %FindTypicalGRB(parameter_optimal[0], parameter_optimal[1], parameter_optimal[2], parameter_optimal[3]))


print('cumulataive number below 10**50')
print(N_th_sche2(50., parameter_optimal[0], parameter_optimal[1], parameter_optimal[2], parameter_optimal[3]))
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Write out the optimazied values                                      |
# if they pass a chi-square two-tailed test at 90% confidnece level    |
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#
# if parameter_optimal[0]>std_dv[0] and parameter_optimal[1]>std_dv[1] and parameter_optimal[2]>std_dv[2] and parameter_optimal[3]>std_dv[3]:
#     if p>0.05 and p<0.95 and p!='nan':
#         with open("optimized.txt", mode = 'a') as fh:
#             fh.write('%.3f %.3f %.2f' %(parameter_optimal[1], parameter_optimal[2], parameter_optimal[3]))
#             fh.write('\n')
#             fh.close()


