import numpy as np
from scipy import integrate, special
from math import exp, pi, log, log10, gamma, sqrt, cos, pi
import sys
import matplotlib.pyplot as plt
from scipy.stats import norm
import scipy.optimize
from pynverse import inversefunc
from matplotlib.ticker import MultipleLocator


# user defined input, that is supposed to be SGRB meta data
file_id = sys.argv[1]
# user defined input, that is the bin of the cumulative histogram
binwidth = round(float(sys.argv[2]), 1)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# making cumulative hist data about logarism of luminosity and redshift  |
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
    # the name of a meta data file
    file_id_str = str(file_id)
    # import the file
    raw_data_matrix = np.loadtxt(file_id_str, dtype = float)
    # the list of T90 of SGRBs
    duration_vector = raw_data_matrix[:,0]  #[s]
    # the list of fluence
    fluence_vector = raw_data_matrix[:,1] # [10^(-7)erg/cm^2]
    # the list of redshift
    redshift_vector = raw_data_matrix[:,2]
    # the empty list for log10 of luminosity
    log10_luminosity_vector = [0 for row in range(len(redshift_vector))]
    
    
    
    #+++++++++++++++++++++++++++++++++
    #                 I(z)           |
    #+++++++++++++++++++++++++++++++++
    # speed of light
    c = 3.0*10**8   # [m/s]
    # the Hubble constant
    H_o = 67.89/(3.0857*10**19) # [km/s/Mpc]*[Mpc/km]= [1/s]
    # the matter density
    omega_M = 0.308
    # the dark energy density
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
    (hist_log10lumino_vector, binedge_log10lumino_vector) = np.histogram(log10_luminosity_vector, bins=np.arange(-250., 250., binwidth), density=False)
    # round the list to reduce the memory
    binedge_log10lumino_vector = [round(i, 2) for i in binedge_log10lumino_vector]
    # the empty list for the center value of the bins
    center_of_bin_1og10lumino_vector = [0 for row in range(len(binedge_log10lumino_vector)-1)]
    # round the values of the list
    #center_of_bin_1og10lumino_vector = [round(i, 2) for i in center_of_bin_1og10lumino_vector]
    
    #  calculate the center values of the bins
    for index in range(len(binedge_log10lumino_vector)-1):
        center_of_bin_1og10lumino_vector[index] += (binedge_log10lumino_vector[index] + binedge_log10lumino_vector[index+1])/2
    
    # make the cumulative number
    cum_hist_log10lumino_vector = [0 for row in range(len(hist_log10lumino_vector))]
    for index in range(len(hist_log10lumino_vector)):
        for sum_num in range(0, index+1):
            cum_hist_log10lumino_vector[index] += hist_log10lumino_vector[sum_num]




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

(center_of_bin_1og10lumino_vector, cum_hist_log10lumino_vector, log10_luminosity_vector, center_of_bin_redshift_vector, cum_hist_redshift_vector, redshift_vector) = make_luminohist_zhist()


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
log10L_vector = center_of_bin_1og10lumino_vector
N_ex_vector =  cum_hist_log10lumino_vector # extract a set of cumulative number





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#      Theroritical nubemr with Guetta & Pirau  distribution                |
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# def N_th_GP(t, N_o, L_o, a, b):
#     if t < L_o:
#         output = N_o*pow((t/L_o), (-1.0)*a)
#     else:
#         output = N_o*pow((t/L_o), (-1.0)*b)
#     return output
#


#
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# #      Plot of  Guetta & Pirau  distribution                                |
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#
#
#
# def plot(C, L_o, a, b):
#     logL_vector = np.arange(43., 60.0, 0.1)
#     L_vector = [0 for row in range(len(logL_vector))]
#     y_vector = [0 for row in range(len(logL_vector))]
#     for index, logL in enumerate(logL_vector, 0):
#         L_vector[index] += pow(10, logL)
#         #y_vector[index] += N_th_GP(L_vector[index], N_o, L_o, a, b)
#         y_vector[index] += cum_N_th_GP_ana(L_vector[index], C, L_o, a, b)
#     plt.plot(logL_vector, y_vector)
#     plt.hist(make_lumino_cumhist_value(), bins=bin_num_lumino, range=(cutoff_minLumino, cutoff_maxLumino), color='red', alpha=0.5)
#     plt.xlabel('logL')
#     plt.ylabel('number')
#     plt.title('broken power')
#     plt.grid(True)
#     plt.savefig("test.png")
#     plt.show()
#     return 0
#
#


#++++++++++++++++++++++++++++++++++++++++++++++++
#       experimental cumulative number          |
#++++++++++++++++++++++++++++++++++++++++++++++++

def N_ex(i):
    return float(N_ex_vector[i])


def cum_N_th_GP_ana1(log10L, C, log10L_o ,a, b): # delta1 and delta2 are constant
    '''
    :param log10L: the list of log10 of L
    :param C: a LF parameter
    :param log10L_o: a LF parameter
    :param a: a LF parameter
    :param b: a LF parameter
    :return: cumulative number of the LF
    '''
    L = pow(10., log10L)
    L_o = pow(10., log10L_o)
    delta1 = pow(10., log10L_o - cutoff_minLumino)
    delta2 = pow(10., cutoff_maxLumino - log10L_o)

    if L < L_o/delta1:
        output = 0.
    elif L_o/delta1 <= L < L_o:
        output = C/(1.0-a)*(pow(L/L_o, 1.0-a)-pow(1./delta1, 1-a))
    elif L_o <= L < delta2*L_o:
        output = C/(1.0-a)*(1.-pow(1./delta1, 1.0-a))+C/(1.0-b)*(pow(L/L_o, 1.0-b) - 1.0)
    elif delta2*L_o <= L:
        output = C/(1.0-a)*(1.0-pow(1/delta1, 1.0-a))+C/(1.0-b)*(pow(delta2, 1.0-b) - 1.0)
    return output


def cum_N_th_GP_ana2(log10L, C, log10L_o ,a, b, log10delta1, log10delta2): # delta1 and delta2 are free
    '''
    :param log10L: the list of log10 of L
    :param C: a LF parameter
    :param log10L_o: a LF parameter
    :param a: a LF parameter
    :param b: a LF parameter
    :param log10delta1: a LF parameter
    :param log10delta2: a LF parameter
    :return: cumulative number of the LF
    '''
    L = pow(10., log10L)
    L_o = pow(10., log10L_o)
    delta1 = pow(10., log10delta1)
    delta2 = pow(10., log10delta2)


    if L < L_o/delta1:
        output = 0.
    elif L_o/delta1 <= L < L_o:
        output = C/(1.0-a)*(pow(L/L_o, 1.0-a)-pow(1./delta1, 1-a))
    elif L_o <= L < delta2*L_o:
        output = C/(1.0-a)*(1.-pow(1./delta1, 1.0-a))+C/(1.0-b)*(pow(L/L_o, 1.0-b) - 1.0)
    elif delta2*L_o <= L:
        output = C/(1.0-a)*(1.0-pow(1/delta1, 1.0-a))+C/(1.0-b)*(pow(delta2, 1.0-b) - 1.0)
    return output




def cum_N_th_GP_ana3(log10L, C, log10L_o ,a, b, log10delta1): # delta1 is free, delta2 is constant
    '''
    :param log10L: he list of log10 of L
    :param C: a LF parameter
    :param log10L_o: a LF parameter
    :param a: a LF parameter
    :param b: a LF parameter
    :param log10delta1: a LF parameter
    :return: cumulative number of the LF
    '''
    log10delta2 = 3.
    if log10L < log10L_o - log10delta1:
        output = 0.
    elif log10L_o - log10delta1 <= log10L < log10L_o:
        output = C/(1.0-a)*(10.**((log10L-log10L_o)*(1.0-a))-10.**((a-1.)*(log10delta1)))
    elif log10L_o <= log10L < log10delta2 + log10L_o:
        output = C/(1.0-a)*(1.-10.**((a-1.)*log10delta1))+C/(1.0-b)*(10.**((log10L-log10L_o)*(1.-b)) - 1.0)
    elif log10delta2+log10L_o <= log10L:
        output = C/(1.0-a)*(1.0-10.**((a-1.)*log10delta1))+C/(1.0-b)*(10.**((1-b)*log10delta2) - 1.0)
    return output


#==================================
def cum_N_th_GP_ana4(log10L, C, log10L_min, g, a, b):
    '''
    :param log10L: he list of log10 of L
    :param C: a LF parameter
    :param log10L_o: a LF parameter
    :param a: a LF parameter
    :param b: a LF parameter
    :param log10delta1: a LF parameter
    :return: cumulative number of the LF
    '''
    L = 10**log10L
    L_min = 10**log10L_min
    log10L_s = 52.30
    log10L_ss = 49.70
    log10L_max = 54.30
    L_s = 10**log10L_s
    L_ss = 10**log10L_ss
    L_max = 10**log10L_max
    if log10L < log10L_min:
        output = 0.
    elif log10L_min <= log10L < log10L_ss:
        output = C/(1-g)*(L**(1-g)-L_min**(1-g))
    elif log10L_ss <= log10L < log10L_s:
        value1 = C/(1-g)*(L_ss**(1-g)-L_min**(1-g))
        value2 = C*L_ss**a*L_s**(-g)/(1-a)*(L**(1-a)-L_ss**(1-a))
        output = value1 + value2
    elif log10L_s <= log10L < log10L_max:
        value1 = C/(1-g)*(L_ss**(1-g)-L_min**(1-g))
        value2 = C*L_ss**a*L_s**(-g)/(1-a)*(L_s**(1-a)-L_ss**(1-a))
        value3 = C*L_ss**a*L_s**(-a-g+b)/(1-b)*(L**(1-b)-L_s**(1-b))
        output = value1 + value2 + value3
    elif log10L_s+2 <= log10L:
        value1 = C/(1-g)*(L_ss**(1-g)-L_min**(1-g))
        value2 = C*L_ss**a*L_s**(-g)/(1-a)*(L_s**(1-a)-L_ss**(1-a))
        value3 = C*L_ss**a*L_s**(-a-g+b)/(1-b)*(L_max**(1-b)-L_s**(1-b))
        output = value1 + value2 + value3
    return output



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  delta1 is a free, delta2 is fixed parameter         |
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def find_fitted_parameters_by_curve_fit():
    '''
    fitting of the LF against data
    :return: fitted parameter and their errors
    '''
    paramater_initial = np.array([5., 51., 0.5, 1.5, 2.])
    param_bounds = ([1., 46., -0.5, 0.1, 0.1], [ 12.0 , 52., 3.0, 10.0, 20.])
    paramater_optimal, covariance = scipy.optimize.curve_fit(np.vectorize(cum_N_th_GP_ana3), log10L_vector, N_ex_vector, p0=paramater_initial, bounds=param_bounds)
    std_dv = np.sqrt(np.diag(covariance))
    return paramater_optimal, std_dv
   

# run the function defined above
parameter_optimal, std_dv = find_fitted_parameters_by_curve_fit()
#+++++++++++++++++++++++++++++++++++++++++++++++++
# print out the fitted parameters
print (parameter_optimal, 3.)
print (std_dv)




#++++++++++++++++++++++++++++++++++++
#               plot                |
#++++++++++++++++++++++++++++++++++++
def plot1(C, log10L_o, a, b): # delta1 and delta2 are constant
    '''
    :param C: a LF parameter
    :param log10L_o: a LF parameter
    :param a: a LF parameter
    :param b: a LF parameter
    :return: 0
    '''
    log10L_vector = np.arange(43., 60.0, 0.1)
    y_vector = [0 for row in range(len(log10L_vector))]
    for index, log10L in enumerate(log10L_vector, 0):
        y_vector[index] += cum_N_th_GP_ana1(log10L, C, log10L_o, a, b)
    plt.plot(log10L_vector, y_vector)
    plt.hist(make_lumino_cumhist_value(), bins=bin_num_lumino, range=(-250, 250), color='red', alpha=0.5)
    plt.xlabel('logL')
    plt.ylabel('number')
    plt.title('broken power')
    plt.xlim([45, 55])
    plt.grid(True)
    plt.savefig("test.png")
    plt.show()
    return 0

def plot2(C, log10L_o, a, b, log10delta1, log10delta2): # delta1 and delta2 are free parameters
    '''
    :param C: a LF parameter
    :param log10L_o: a LF parameter
    :param a: a LF parameter
    :param b: a LF parameter
    :param log10delta1: a LF parameter
    :param log10delta2: a LF parameter
    :return: 0
    '''
    log10L_vector = np.arange(43., 60.0, 0.1)
    y_vector = [0 for row in range(len(log10L_vector))]
    for index, log10L in enumerate(log10L_vector, 0):
        y_vector[index] += cum_N_th_GP_ana2(log10L, C, log10L_o, a, b, log10delta1, log10delta2)
    plt.plot(log10L_vector, y_vector)
    plt.hist(make_lumino_cumhist_value(), bins=bin_num_lumino, range=(-250, 250), color='red', alpha=0.5)
    plt.xlabel('logL')
    plt.ylabel('number')
    plt.title('broken power')
    plt.xlim([45, 55])
    plt.grid(True)
    plt.savefig("test.png")
    plt.show()
    return 0




def plot3(C, log10L_o, a, b, log10delta1): #  delta2 is constant
    '''

    :param C: a LF parameter
    :param log10L_o: a LF parameter
    :param a: a LF parameter
    :param b: a LF parameter
    :param log10delta1: a LF parameter
    :return: 0
    '''
    log10L_vector = np.arange(43., 60.0, 0.1)
    y_vector = [0 for row in range(len(log10L_vector))]
    for index, log10L in enumerate(log10L_vector, 0):
        y_vector[index] += cum_N_th_GP_ana3(log10L, C, log10L_o, a, b, log10delta1)
    plt.plot(log10L_vector, y_vector)
    plt.hist(make_lumino_cumhist_value(), bins=np.arange(-250., 250., binwidth), color='red', alpha=0.5)
    plt.xlabel('logL')
    plt.ylabel('number')
    plt.title('broken power')
    plt.xlim(46, 53)
    plt.grid(True)
    #plt.savefig("test.png")
    plt.show()
    return 0




#+++++++++++++++++++++++++++++++++++
# find Luminosity of a typical GRB |
#+++++++++++++++++++++++++++++++++++
def FindTypicalGRB(C, log10L_o, a, b, log10delta1):
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
        return cum_N_th_GP_ana3(log10L, C, log10L_o, a, b, log10delta1)
    #inverse function of cumNwithOptV, giving a typical L_iso
    typL_iso = inversefunc(cumNwithOptV, domain=[48.5, 52.])
    return typL_iso(totalnum/2.)



# plot it
plot3(parameter_optimal[0], parameter_optimal[1], parameter_optimal[2], parameter_optimal[3], parameter_optimal[4]) # delta2 is constant



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
        if fun(x, parameter_optimal[0], parameter_optimal[1], parameter_optimal[2], parameter_optimal[3], parameter_optimal[4]) > 0.0:
            non_zero_comp_num += 1
            N_th_vector.append(fun(x, parameter_optimal[0], parameter_optimal[1], parameter_optimal[2], parameter_optimal[3], parameter_optimal[4]))

    #print (N_th_vector)
    #print (y_step_vector[len(y_step_vector)-len(N_th_vector):])
    #print (5 -(len(y_step_vector) - non_zero_comp_num))
    # the degree of freedom is (the # of user preferd step) + (# 2 parts of flat) - 1 - (# of parameters), which is equivalent to (# of non zero data points) - 1 - (# of constraint)
    
    from scipy.stats import chisquare
    chisq, pvalue = chisquare(y_step_vector[len(y_step_vector)-len(N_th_vector):], np.array(N_th_vector), ddof=5 - (len(y_step_vector)- non_zero_comp_num))
    return chisq, non_zero_comp_num, pvalue


# print out the chi-square test
print ("the value of the chi square, the number of data points with non-zero value (N_th)")
chisquare, NonzeroNth, p = chi_sqr_test(cum_N_th_GP_ana3, log10L_vector, N_ex_vector)
print(chisquare, NonzeroNth, p)



def find_coverage_of_sample(log10Lmin, log10Lmax):
    '''
    find the how much LF covers the SGRB samples
    :param log10Lmin: the minimum value of log10 of L
    :param log10Lmax: the maximum value of log10 of L
    :return:
        value1: the minimum threshold of the luminosity with hist-fitted normal dist
        value2: the maximum threshold of the luminosity with hist-fitted normal dist
        value2 - value1: the value between
    '''
    (hist_log10lumino_vector1, binedge_log10lumino_vector1) = np.histogram(log10_luminosity_vector, bins=np.arange(-250., 250., binwidth) ,density=False) # extract the the # of the hist  and the edge of the bin
    binedge_log10lumino_vector1 = [round(i, 2) for i in binedge_log10lumino_vector1]
    center_of_bin_1og10lumino_vector1 = [0 for row in range(len(binedge_log10lumino_vector1)-1)]
    for index in range(len(binedge_log10lumino_vector1)-1):
        center_of_bin_1og10lumino_vector1[index] += (binedge_log10lumino_vector1[index] + binedge_log10lumino_vector1[index+1])/2 # compute the value of the center of the bin

    center_of_bin_1og10lumino_vector1 = [round(i, 2) for i in center_of_bin_1og10lumino_vector1]

    for index, y in enumerate(hist_log10lumino_vector1, 0):
        if y > 0.:
            starting_non_zero_index = index
            break
    for index, y in enumerate(hist_log10lumino_vector1[starting_non_zero_index:], 0):
        if y == 0.:
            ending_zero_index = index+starting_non_zero_index
            break
    hist_log10lumino_vector1_nonzero_density = [i/float(len(log10_luminosity_vector)) for i in hist_log10lumino_vector1[starting_non_zero_index:ending_zero_index]]

    # find the normal distribution by least square method applied the histogram
    def func(x, mu, std):
        return norm.pdf(x, mu, std)

    parameter_initial_hist = np.array([49.0, 0.5])
    parameter_optimal_hist, covariance = scipy.optimize.curve_fit(func, center_of_bin_1og10lumino_vector1[starting_non_zero_index:ending_zero_index], hist_log10lumino_vector1_nonzero_density, p0=parameter_initial_hist)


    # find the minimum and maximum threshold of the luminosity with hist-fitted normal dist
        
    value1, a1 = integrate.quad(lambda t: norm.pdf(t, parameter_optimal_hist[0], parameter_optimal_hist[1]), -np.inf, log10Lmin)
    value2, a2 = integrate.quad(lambda t: norm.pdf(t, parameter_optimal_hist[0], parameter_optimal_hist[1]), -np.inf, log10Lmax)
    
   
    plt.plot(center_of_bin_1og10lumino_vector1[starting_non_zero_index:ending_zero_index], hist_log10lumino_vector1_nonzero_density)
    plt.xlim([43., 58])
    a = np.linspace(45., 55, 100)
    plt.plot(a, norm.pdf(a, parameter_optimal_hist[0],parameter_optimal_hist[1]))
    plt.scatter([log10Lmin, log10Lmax], [norm.pdf(log10Lmin, parameter_optimal_hist[0], parameter_optimal_hist[1]), norm.pdf(log10Lmax, parameter_optimal_hist[0], parameter_optimal_hist[1])], color='r')
    plt.show()
    return value1, value2, value2 - value1

print ("coverage of the sample")
print (find_coverage_of_sample(parameter_optimal[1] - parameter_optimal[4], parameter_optimal[1] + 3.))
# find the luminosity of a typical GRB 
print('Typical on-axis isotropic equivalent logl10 is %.3f' % FindTypicalGRB(parameter_optimal[0], parameter_optimal[1], parameter_optimal[2], parameter_optimal[3], parameter_optimal[4]))

print('cumulataive number below 10**50')
print(cum_N_th_GP_ana3(50., parameter_optimal[0], parameter_optimal[1], parameter_optimal[2], parameter_optimal[3], parameter_optimal[4]))


