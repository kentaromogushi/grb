import numpy as np
#import os


#rate_function_list = ["delayed_hernquist", "delayed_porciani", "delayed_fardal", "delayed_cole", "delayed_hopkins", "delayed_wilkins"]



#++++++++++++++
# RF and tmin |
#++++++++++++++
rate_function_list = ["delayed_hernquist", "delayed_hopkins"]
t_min_list = ["1.0", "0.1", "0.02"] # 1Gyr, 100Myr, 20Myr


# #==========check for the theta_c distribution ========
# rate_function_list = ["delayed_hernquist"]
# t_min_list = ["1.0"] # 1Gyr
# #=====================================================


#power_time_pro_list = ["1.0", "1.2", "1.5"] # these are negative power
#matterEng_density_list = ["0.308" + " " + "0.692" , "0.32" + " " + "0.68", "0.296" + " " + "0.704"] # without error, matter up, matter down


#bin_num_list = ["5", "6", "7", "8", "9", "10"]
bin_num_list = binwidth = ["0.3"]

#++++++++++++++
# theta_obsGW |
#+++++++++++++
theta_obsGW_list = ['%.5f' %(np.pi/180*i) for i in np.arange(15., 37., 1.)] # radian



# #=========== check for the theta_c distribution ========
# theta_obsGW_list = ['%.5f' %(np.pi/180*i) for i in [5.]] # radian
# #=======================================================

# #========== add theta_obsGW = [2, 14] degs ==================
# theta_obsGW_list = ['%.5f' %(np.pi/180*i) for i in np.arange(2., 15., 1.)] # radian
# #======================================================

# #=========== add theta_obsGW = [12, 13, 14, 37, 38, 39, 40] ========
# theta_obsGW_list = ['%.5f' %(np.pi/180*i) for i in [12., 13., 14., 37., 38., 39., 40.]] # radian
# #====================================================================

# #=========== add theta_obsGW = [37, 38, 39, 40] ========
# theta_obsGW_list = ['%.5f' %(np.pi/180*i) for i in [37., 38., 39., 40.]] # radian
# #====================================================================

#++++++++++
# theta_c |
#++++++++++
def theta_c_list_make(theta_obsGW):
    theta_obsGW = round(180./np.pi*float(theta_obsGW)) # convert to radian to degree
    return ['%.5f' %(np.pi/180.*i) for i in np.arange(1., theta_obsGW, 1.)]

# # ==============check for the theta_c distribution ========
# def theta_c_list_make(theta_obsGW):
#     theta_obsGW = round(180./np.pi*float(theta_obsGW)) # convert to radian to degree
#     return ['%.5f' %(np.pi/180.*i) for i in np.arange(1., theta_obsGW, 1.)]
# #===========================================================

# #============= add theta_c = 1, 2 degs ===========
# theta_c_list = ['%.8f' %(np.pi/180.*i) for i in [0.1]]
# #================================

# #================= check theta_obsGW = 5 degs=================
# def theta_c_list_make(theta_obsGW):
#     theta_obsGW = round(180./np.pi*float(theta_obsGW)) # convert to radian to degree
#     return ['%.5f' %(np.pi/180.*i) for i in np.arange(1., theta_obsGW, 0.2)]
# #=============================================================

# #============= add theta_obsGW = [2, 14] degs ===========
# def theta_c_list_make(theta_obsGW):
#     theta_obsGW = round(180./np.pi*float(theta_obsGW)) # convert to radian to degree
#     return ['%.5f' %(np.pi/180.*i) for i in np.linspace(1., theta_obsGW, 18., endpoint=False)]
# #===================================================


#++++++++++++++++++++++++++++++++++++
# Karelle schechter with free delta |++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++


Karelle_revised_schechter_list = ["51.474" + " " + "0.724" + " " + "3.46",
                                  "51.194" + " " + "0.669" + " " + "3.61",
                                  "51.278" + " " + "0.689" + " " + "3.57",
                                  "51.301" + " " + "0.745" + " " + "3.54",
                                  "51.203" + " " + "0.682" + " " + "3.90",
                                  "51.349" + " " + "0.796" + " " + "3.44"] # witdth = 0.5, 0.6, 0.7, 0.8, 0.9, 1.0




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





#++++++++++++++++++++++++++++++++++++++++++++++++++++++
#   broken power with free delta1 and fixed delta2(1000)         |
#++++++++++++++++++++++++++++++++++++++++++++++++++++++

Karelle_revised_broken_list = ["51.236" + " " + "0.783" + " " + "2.50" + " " + "3.14" + " " + "3.",
                              "50.836" + " " + "0.702" + " " + "2.11" + " " + "3.24" + " " + "3.",
                              "51.112" + " " + "0.762" + " " + "2.71" + " " + "3.28" + " " + "3."] # witdth = 0.5, 0.6, 0.7


#March 16th 2018, I ran it again then, I got the different values below!!!!!!
# However I found this is wrong on April 8th 2018. So the above one is correct!!!!
# Karelle_revised_broken_list = ["49.991" + " " + "0.758" + " " + "2.12" + " " + "2.87" + " " + "3.",
#                                "49.944" + " " + "0.747" + " " + "2.19" + " " + "2.95" + " " + "3.",
#                                "49.884" + " " + "0.702" + " " + "2.40" + " " + "3.13" + " " + "3."] # witdth = 0.5, 0.6, 0.7







#+++++++++++++++++++++++++++++++++++++++++
#         schechter & Karelle revised       |
#+++++++++++++++++++++++++++++++++++++++++




def schechterKarellerevised():
    for j, t_min in enumerate(t_min_list, 1):
        for k, rate_function in enumerate(rate_function_list, 2):
            for l, lumino_para in enumerate(Karelle_revised_schechter_list, 1):
                for n, bin_num in enumerate(bin_num_list, 1):  # fix the number of bin for redshift in this case
                    for m, theta_obsGW in enumerate(theta_obsGW_list, 1):
                    #for m, theta_obsGW in enumerate(theta_obsGW_list, -30):
                        for p, theta_c in enumerate(theta_c_list_make(theta_obsGW), 1): 
                        #for p, theta_c in enumerate(theta_c_list_make(theta_obsGW), -60): 

                            if k == 3:
                                k = 6
                            if t_min == "1.0":
                                t_min_label = "1Gyr"
                            if t_min == "0.1":
                                t_min_label = "100Myr"
                            if t_min == "0.02":
                                t_min_label = "20Myr"
                            launching = "python allobtain_for_swift_schechter_jet.py Karelle_revised_schechterOUTPUT.txt " + str(j)+","+str(k)+","+str(l)+","+str(n)+","+str(m)+","+str(p) + " " + "Karelle_list_revised_updated.txt " + rate_function + " " + lumino_para + " " + "0.308" + " " + "0.692" + " " + t_min + " " + "1.0" + " " + bin_num + " " + theta_obsGW + " " + "F_pow" +  " " + theta_c
                            with open("Karelle_revised_schechter.sh", mode = 'a') as fh:
                                fh.write(launching)
                                fh.write('\n')
                            fh.close()












#+++++++++++++++++++++++++++++++++++++++++
#         broken power & Karelle revised |
#+++++++++++++++++++++++++++++++++++++++++

def brokenKarellerevised():
    for j, t_min in enumerate(t_min_list, 1):
        for k, rate_function in enumerate(rate_function_list, 2):
            for l, lumino_para in enumerate(Karelle_revised_broken_list, 1):
                for n, bin_num in enumerate(bin_num_list, 1):  # fix the number of bin for redshift in this case
                    for m, theta_obsGW in enumerate(theta_obsGW_list, 1):
                    #for m, theta_obsGW in enumerate(theta_obsGW_list, -30):
                        for p, theta_c in enumerate(theta_c_list_make(theta_obsGW), 1): 
                        #for p, theta_c in enumerate(theta_c_list_make(theta_obsGW), -60): 
                            if k == 3:
                                k = 6
                            if t_min == "1.0":
                                t_min_label = "1Gyr"
                            if t_min == "0.1":
                                t_min_label = "100Myr"
                            if t_min == "0.02":
                                t_min_label = "20Myr"
                            launching = "python allobtain_for_swift_GP_jet.py Karelle_revised_brokenOUTPUT.txt " + str(j)+","+str(k)+","+str(l)+","+str(n)+","+str(m)+","+str(p) + " " + "Karelle_list_revised_updated.txt " + rate_function + " " + lumino_para + " " + "0.308" + " " + "0.692" + " " + t_min + " " + "1.0" + " " + bin_num + " " + theta_obsGW + " " + "F_pow" + " " + theta_c
                            #with open("Karelle_revised_broken.sh", mode = 'a') as fh:
                            with open("Karelle_revised_broken.sh", mode = 'a') as fh:
                                fh.write(launching)
                                fh.write('\n')
                            fh.close()







if __name__=="__main__":

#    schechterKarellerevised()
#    brokenKarellerevised()

