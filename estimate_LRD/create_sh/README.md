#make .sh file to run the code allobtain_for_swift_schechter_jet.py or allobtain_for_swift_GP_jet.py with different sets of parameters

### description 
1. edit run_get_multi_fixedwidth.py 
    1. if one wants to study with a different list of \theta_c and \theta_obsGW where \theta_c is the half opening angle of SGRB jet and \theta_obsGW is GW inclination angel 
    2. uncomment schechterKarellerevised() for the Schecher LF, brokenKarellerevised() for the broken power LF 
2. run by typing run_get_multi_fixedwidth.py 
3. This creates a file called Karelle_revised_schechter.sh or Karelle_revised_broken.sh
4. The each line of the file is for running allobtain_for_swift_schechter_jet.py or allobtain_for_swift_GP_jet.py in the above directory 
5. chmod -x Karelle_revised_schechter.sh or chmod -x Karelle_revised_broken.sh 
6. type ./Karelle_revised_schechter.sh or ./Karelle_revised_broken.sh

### requiremets 
1. To run these bash scripts, one needs allobtain_for_swift_schechter_jet.py or  allobtain_for_swift_GP_jet.py along with Karelle_list_revised_updated.txt in a same directory
