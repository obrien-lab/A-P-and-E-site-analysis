import os
import pickle
import pandas as pd
import sys
import statistics

CODON_TYPES = ['UUU', 'UUC', 'UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG', 'AUU', 'AUC', 'AUA', 'AUG', 'GUU', 'GUC', 'GUA',
              'GUG', 'UCU', 'UCC', 'UCA', 'UCG', 'CCU', 'CCC', 'CCA', 'CCG', 'ACU', 'ACC', 'ACA', 'ACG', 'GCU', 'GCC',
              'GCA', 'GCG', 'UAU', 'UAC', 'CAU', 'CAC', 'CAA', 'CAG', 'AAU', 'AAC', 'AAA', 'AAG', 'GAU', 'GAC', 'GAA',
              'GAG', 'UGU', 'UGC', 'UGG', 'CGU', 'CGC', 'CGA', 'CGG', 'AGU', 'AGC', 'AGA', 'AGG', 'GGU', 'GGC', 'GGA',
              'GGG', 'UAA', 'UAG', 'UGA']


def coopravity_and_no_cooprativity_effect(file,input,filename,outpath):

    time_file = pd.read_csv(input+file, sep='\t')
    pvalue_list = []    
    count = 0
    
    pkl_file = 'pkl_file/'
    LIST_FOR_PERMUTATION = {}
    COOPRATIVITY_LIST = {}
    for asite_aa in CODON_TYPES:
        if asite_aa in ['UAA', 'UAG', 'UGA']:
            continue
        COOPRATIVITY_LIST[asite_aa] = {}
        for psite_aa in CODON_TYPES:
            if psite_aa in ['UAA', 'UAG', 'UGA']:
                continue
            LIST_FOR_PERMUTATION[psite_aa] = {}
            COOPRATIVITY_LIST[asite_aa][psite_aa] = {}
            for esite_aa in CODON_TYPES:
                if esite_aa in ['UAA', 'UAG', 'UGA']:
                    continue
                LIST_FOR_PERMUTATION[psite_aa][esite_aa] = []
                COOPRATIVITY_LIST[asite_aa][psite_aa][esite_aa] = []
                a_site = time_file.Asite == asite_aa
                p_site = time_file.psite == psite_aa
                alt_p_site = time_file.psite != psite_aa
                e_site = time_file.esite == esite_aa
                alt_e_site = time_file.esite != esite_aa
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

                # (x,..,z) = a_e_pair
                a_e_pair_list = a_e_pair_dataframe.values.tolist()
                a_e_pair_list_for_permutation = a_e_pair_list

                # (..,..,z) = a_pair
                a_pair_list = a_pair_dataframe.values.tolist()
                a_pair_list_for_permutation = a_pair_list

                # (x,y,z) = a_p_e_pair
                a_p_e_pair_list = a_p_e_pair_dataframe.values.tolist()
                a_p_e_pair_list_for_permutation = a_p_e_pair_list
                
                
                LIST_FOR_PERMUTATION[psite_aa][esite_aa] = (a_p_pair_list_for_permutation,a_e_pair_list_for_permutation,a_pair_list_for_permutation,a_p_e_pair_list_for_permutation,delta)            
                COOPRATIVITY_LIST[asite_aa][psite_aa][esite_aa] = (delta1,delta2,delta3,len(a_p_pair_dataframe),len(a_e_pair_dataframe),len(a_p_e_pair_dataframe),len(a_pair_dataframe))
                count+=1
                print(count)
        pickle.dump(LIST_FOR_PERMUTATION, open(outpath+filename+'_LIST_FOR_PERMUTATION_dic_'+asite_aa+'.p', 'wb'))
    pickle.dump(COOPRATIVITY_LIST, open(path+'/Permutation_test/'+filename+'_COOPRATIVITY_LIST.p', 'wb'))
    return




data_name = sys.argv[1]
cwd = os.getcwd()
path = cwd

dataset = ['Jan','Williams','Young','Weinberg','Nissley1','Nissley2','Pool']

input1 = path+'/Result/esite_asite_matrix_codon/'
print(input1)
files = os.listdir(input1)
print(files)
jobs = []

outpath = path+'/Permutation_test/pkl_file/'+data_name+'/'

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
