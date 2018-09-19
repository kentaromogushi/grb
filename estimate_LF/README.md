# Estimate luminosity function 

## discription 
1. chi_sqare_schechter-v2.py is the code to estimate the Schechter LF 
2. chi_sqare_GueetaPiran-v2.py is the code to estimate the broken power LF 
3. Karelle_list_revised_updated.txt is the meta data about SGRBs



### the way to run a code 
1. for the Schechter LF
    1. python chi_sqare_schechter-v2.py Karelle_list_revised_updated.txt {BINNUM}
2. for the broken power LF 
    1. python chi_sqare_GueetaPiran-v2.py Karelle_list_revised_updated.txt {BINNUM}
    
2. {BINNUM} = (0.5, 0.6, 0.7, 0.8, 0.9, 1.)is the user defined the number of bin in logL space 


 ### ouput 
 1. fitted paramters
 2. uncertainty of each parameter
 3. chi-sqare, the number of data point excluding zero values, and p-value 
 4. median value of L 