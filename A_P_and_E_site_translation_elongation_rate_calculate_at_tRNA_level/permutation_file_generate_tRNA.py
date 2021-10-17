dataset = ['Jan','Williams','Young','Weinberg','Nissley1','Nissley2','Pool']

import os
#import cPickle as Pickle
import pickle
import pandas as pd
import sys
import statistics


tRNA_list = ['A-t1', 'A-t2', 'R-t1', 'R-t2', 'R-t3', 'R-t4', 'D-t1', 'N-t1', 'C-t1', 'E-t1', 'E-t2', 'Q-t1', 'Q-t2', 'G-t1', 'G-t2', 'G-t3', 'H-t1', 'I-t1', 'I-t2', 'L-t1', 'L-t2', 'L-t3', 'L-t4', 'K-t1', 'K-t2', 'M-t1', 'F-t1', 'P-t1', 'P-t2', 'S-t1', 'S-t2', 'S-t3', 'S-t4', 'T-t1', 'T-t2', 'T-t3', 'W-t1', 'Y-t1', 'V-t1', 'V-t2', 'V-t3', '*-t1']


def coopravity_and_no_cooprativity_effect(file,input,filename,outpath):
    time_file = pd.read_csv(input+file, sep='\t')
    pvalue_list = []    
    count = 0    
    LIST_FOR_PERMUTATION = {}
    COOPRATIVITY_LIST = {}
    for asite_aa in tRNA_list:
        if asite_aa == '*-t1':
            continue
        LIST_FOR_PERMUTATION[asite_aa] = {}
        COOPRATIVITY_LIST[asite_aa] = {}
        for psite_aa in tRNA_list:
            if psite_aa == '*-t1':
                continue
            LIST_FOR_PERMUTATION[asite_aa][psite_aa] = {}
            COOPRATIVITY_LIST[asite_aa][psite_aa] = {}
            for esite_aa in tRNA_list:
                if esite_aa == '*-t1':
                    continue
                LIST_FOR_PERMUTATION[asite_aa][psite_aa][esite_aa] = []
                COOPRATIVITY_LIST[asite_aa][psite_aa][esite_aa] = []
                a_site = time_file.Asite_tRNA == asite_aa
                p_site = time_file.psite_tRNA == psite_aa
                alt_p_site = time_file.psite_tRNA != psite_aa
                e_site = time_file.esite_tRNA == esite_aa
                alt_e_site = time_file.esite_tRNA != esite_aa
                a_and_p = a_site & p_site & alt_e_site
                a_and_e = a_site & e_site & alt_p_site
                a_p_and_e = a_site & p_site & e_site
                alt_a_site = a_site & alt_e_site & alt_p_site

               
                # (X,Y,Z) are in the E-, P- and A-sites.
                # (x,..,z) = a_e_pair
                # (..,y,z) = a_p_pair
                # (..,..,z) = a_pair
                # (x,y,z) = a_p_e_pair

                a_p_pair_dataframe = time_file[a_and_p].Translation_time_Asite_codon
                a_e_pair_dataframe = time_file[a_and_e].Translation_time_Asite_codon
                a_pair_dataframe = time_file[alt_a_site].Translation_time_Asite_codon
                a_p_e_pair_dataframe = time_file[a_p_and_e].Translation_time_Asite_codon

                try:
                    a_p_pair = statistics.median(time_file[a_and_p].Translation_time_Asite_codon)
                    a_e_pair = statistics.median(time_file[a_and_e].Translation_time_Asite_codon)
                    a_pair = statistics.median(time_file[a_site].Translation_time_Asite_codon)
                    a_p_e_pair = statistics.median(time_file[a_p_and_e].Translation_time_Asite_codon)
                except statistics.StatisticsError:
                    continue

                    
                # delta1 = median p(X,Y,Z) - median p(...,...,Z)
                # delta2 = median p(...,Y,Z) - median p(...,...,Z)
                # delta3 = median p(X,...,Z) - median p(...,...,Z)
 
                
                delta1 = a_p_e_pair- a_pair
                delta2 = a_p_pair- a_pair
                delta3 = a_e_pair- a_pair
                delta4 = delta2+delta3                
                delta = delta1-delta2-delta3
                
                                
                #cooperetivity_condition
                # cooperetivity_condition = 1 - if delta1 = delta2+delta3 (Non-cooperative)
                # cooperetivity_condition = 2 - if delta1 > delta2+delta3 (Positive cooperativity)
                # cooperetivity_condition = 3 - if delta1 < delta2+delta3 (Negative cooperativity)
                
                if delta1 == delta4:
                    cooperetivity_condition = 1
                elif delta1 > delta4:
                    cooperetivity_condition = 2
                elif delta1 < delta4:
                    cooperetivity_condition = 3
                
                a_p_pair_list = a_p_pair_dataframe.values.tolist()
                a_p_pair_list_for_permutation = a_p_pair_list

                a_e_pair_list = a_e_pair_dataframe.values.tolist()
                a_e_pair_list_for_permutation = a_e_pair_list

                a_pair_list = a_pair_dataframe.values.tolist()
                a_pair_list_for_permutation = a_pair_list

                a_p_e_pair_list = a_p_e_pair_dataframe.values.tolist()
                a_p_e_pair_list_for_permutation = a_p_e_pair_list
                
                LIST_FOR_PERMUTATION[asite_aa][psite_aa][esite_aa] = (a_p_pair_list_for_permutation,a_e_pair_list_for_permutation,a_pair_list_for_permutation,a_p_e_pair_list_for_permutation,delta)
            
                COOPRATIVITY_LIST[asite_aa][psite_aa][esite_aa] = (delta1,delta2,delta3,len(a_p_pair_dataframe),len(a_e_pair_dataframe),len(a_p_e_pair_dataframe),len(a_pair_dataframe))
                
    pickle.dump(LIST_FOR_PERMUTATION, open(outpath+filename+'_LIST_FOR_PERMUTATION_dic.p', 'wb'))
    pickle.dump(COOPRATIVITY_LIST, open(outpath+filename+'_COOPRATIVITY_LIST.p', 'wb'))
    
    return



data_name = sys.argv[1]
cwd = os.getcwd()
path = cwd+'/'

input1 = path+'Result/esite_asite_matrix_tRNA/'
print(input1)
files = os.listdir(input1)
print(files)
jobs = []

outpath = path+'Permutation_test/permutation_file_generate/'

try:
    os.makedirs(outpath)
except OSError:
    print("Creation of the directory %s failed")
else:
    print("Successfully created the directory %s ")

times_coordinates = []
for file in files:
    if 'Asite_psite_esite_times_coordinates' in file:
         if data_name in file:
            times_coordinates += [file]
    else:
        pass

print('times_coordinates',times_coordinates)
count = 0
file = times_coordinates[0]
filename = times_coordinates[0]
print(filename)
filename = filename.split('_')
print(filename)
filename = filename[:3]
filename = '_'.join(filename)

coopravity_and_no_cooprativity_effect(file,input1,filename,outpath)
    
print('Done')

