from __future__ import division
import math
import os
import pickle
from optparse import OptionParser
import sys
import random
import statistics
from sklearn.utils import resample
import timeit
import numpy as np


def vector(a,b):
    a_sq = a**2
    b_sq = b**2
    ab_total = a_sq+b_sq
    vector_of_ab = math.sqrt(ab_total)

    return vector_of_ab


def permutation_test(a_p_pair_list_for_permutation,a_e_pair_list_for_permutation,a_pair_list_for_permutation,a_p_e_pair_list_for_permutation,delta,asite_aa, psite_aa, esite_aa):

    # (X,Y,Z) are in the E-, P- and A-sites.
    # (x,..,z) = a_e_pair
    # (..,y,z) = a_p_pair
    # (..,..,z) = a_pair
    # (x,y,z) = a_p_e_pair
   

    #delta1_permuation_list = []
    #delta2_permuation_list = []
    #delta3_norm_permuation_list = []
    #elta_norm_permuation_list = []
    
    cooperetivity_count = 0
    pos_cooperetivity_count = 0
    neg_cooperetivity_count = 0

    #length calculate
    # (x,..,z) = a_e_pair
    length_of_e_and_a_list = len(a_e_pair_list_for_permutation)
    #print(length_of_e_and_a_list)
    # (..,y,z) = a_p_pair
    length_of_a_p_list = len(a_p_pair_list_for_permutation)
    #print(length_of_a_p_list)
    # (..,..,z) = a_pair
    length_of_a_list = len(a_pair_list_for_permutation)
    #print(length_of_a_list)
    # (x,y,z) = a_p_e_pair
    length_of_a_p_e_list = len(a_p_e_pair_list_for_permutation)
    #print(length_of_a_p_e_list)


    # delta1 = median p(X,Y,Z) - median p(...,...,Z)
    # delta2 = median p(...,Y,Z) - median p(...,...,Z)
    # delta3 = median p(X,...,Z) - median p(...,...,Z)
        
    #all 4 list shuffle and create new list
    combined_list = a_p_e_pair_list_for_permutation + a_p_pair_list_for_permutation + a_e_pair_list_for_permutation + a_pair_list_for_permutation
    list_len = len(combined_list)


    for i in range(10000):
        
        combined_list_shuffled = combined_list
        np.random.shuffle(combined_list_shuffled)
                
        
        new_a_p_e_pair_shuffle_list = combined_list_shuffled[:length_of_a_p_e_list]
        #print('new_a_p_e_pair_shuffle_list',len(new_a_p_e_pair_shuffle_list))
        new_a_p_pair_shuffle_list = combined_list_shuffled[length_of_a_p_e_list:(length_of_a_p_list+length_of_a_p_e_list)]
        #print('new_a_p_pair_shuffle_list',len(new_a_p_pair_shuffle_list))
        new_a_e_pair_shuffle_list = combined_list_shuffled[(length_of_a_p_list+length_of_a_p_e_list):(length_of_a_p_list+length_of_a_p_e_list+length_of_e_and_a_list)]
        #print('new_a_e_pair_shuffle_list',len(new_a_e_pair_shuffle_list))
        new_a_pair_shuffle_list = combined_list_shuffled[(length_of_a_p_list+length_of_a_p_e_list+length_of_e_and_a_list):]
        #print('new_a_pair_shuffle_list',len(new_a_pair_shuffle_list))
        
        delta1_perm = statistics.median(new_a_p_e_pair_shuffle_list)-statistics.median(new_a_pair_shuffle_list)
        delta2_perm = statistics.median(new_a_p_pair_shuffle_list)-statistics.median(new_a_pair_shuffle_list)
        delta3_perm = statistics.median(new_a_e_pair_shuffle_list)-statistics.median(new_a_pair_shuffle_list)
        delta_perm = delta1_perm-delta2_perm-delta2_perm
        
        if delta_perm > delta:
            #elta_perm_val = 1
            pos_cooperetivity_count+=1
                
        if delta_perm < delta:
            #elta_perm_val = 1
            neg_cooperetivity_count+=1
        
        if abs(delta_perm) > abs(delta):
            #elta_perm_val = 1
            cooperetivity_count+=1
        
            
        
        


        
    #non_cooperetivity_count = 0
    #pos_cooperetivity_count = 0
    #neg_cooperetivity_count = 0
    cooperativity_pvalue = cooperetivity_count/10000
    pos_cooperativity_pvalue = pos_cooperetivity_count/10000
    neg_cooperativity_pvalue = neg_cooperetivity_count/10000
    
    
    
    #print('non_cooperativity_pvalue',cooperativity_pvalue)
    #print('pos_cooperativity_pvalue',pos_cooperativity_pvalue)
    #print('neg_cooperativity_pvalue',neg_cooperativity_pvalue)
    

    
    #g_v = len([key for key in delta_norm_permuation_list if key == 0])

    if cooperativity_pvalue < 0.05:
        if pos_cooperativity_pvalue < neg_cooperativity_pvalue:
            pvalue = pos_cooperativity_pvalue
            perm_cooperativity = 'pos_cooperativity'
        elif neg_cooperativity_pvalue < pos_cooperativity_pvalue:
            pvalue = neg_cooperativity_pvalue
            perm_cooperativity = 'neg_cooperativity'
    else:
        pvalue = cooperativity_pvalue
        perm_cooperativity = 'non_cooperativity'
        
    
    
    #print('pvalue',pvalue)
    #print('perm_cooperativity',perm_cooperativity)
    

    #print('Done2')
    #sys.exit()

    return asite_aa, psite_aa, esite_aa, pvalue, perm_cooperativity, cooperativity_pvalue, pos_cooperativity_pvalue, neg_cooperativity_pvalue





data_name = sys.argv[1]
#print(data_name)

cwd = os.getcwd()
asite = sys.argv[2]
#print(asite)
psite = sys.argv[3]
#print(psite)
esite = sys.argv[4]
#print(esite)

path1 = 'Permutation_test/'+data_name+"/"

try:
	os.makedirs(path1)
except OSError:
	pass
	#print("Creation of the directory %s failed")
else:
	pass
	#print("Successfully created the directory %s ")


path = cwd+'/'
pvalue_save = open(path+path1+data_name+'_'+asite+'_'+psite+'_'+esite+'.tab','w')
start = timeit.default_timer()


files = os.listdir(path+'Permutation_test/')
#print('files',files)

PERMUTATION_LIST_FILE = []
for file in files:
	if 'LIST_FOR_PERMUTATION_dic.p' in file:
		PERMUTATION_LIST_FILE.append(file)
	else:
		pass

#print(PERMUTATION_LIST_FILE)

file_name = ''
n = 0
for data_file in PERMUTATION_LIST_FILE:
	n+=1
	#print(data_file)
	if data_name in data_file:
		file_name = data_file

#print(file_name)
PERMUTATION_LIST_file = open(path+'Permutation_test/'+file_name,'rb')
PERMUTATION_LIST_load = pickle.load(PERMUTATION_LIST_file)

#print(PERMUTATION_LIST_load[asite][psite][esite][0])

a_p_pair_list_for_permutation = PERMUTATION_LIST_load[asite][psite][esite][0]
#print(a_p_pair_list_for_permutation)
#print(len(a_p_pair_list_for_permutation))
a_e_pair_list_for_permutation = PERMUTATION_LIST_load[asite][psite][esite][1]
#print(len(a_e_pair_list_for_permutation))
a_pair_list_for_permutation = PERMUTATION_LIST_load[asite][psite][esite][2]
#print(len(a_pair_list_for_permutation))
a_p_e_pair_list_for_permutation = PERMUTATION_LIST_load[asite][psite][esite][3]
#print(len(a_p_e_pair_list_for_permutation))
v = PERMUTATION_LIST_load[asite][psite][esite][4]
#print(v)
PERMUTATION_LIST_file.close()
asite_aa, psite_aa, esite_aa, pvalue, perm_cooperativity, cooperativity_pvalue, pos_cooperativity_pvalue, neg_cooperativity_pvalue = permutation_test(a_p_pair_list_for_permutation, a_e_pair_list_for_permutation, a_pair_list_for_permutation, a_p_e_pair_list_for_permutation, v, asite, psite, esite)
pvalue_save = open(path+path1+data_name+'_'+asite+'_'+psite+'_'+esite+'.tab','a')

pvalue_save.write(str(asite_aa)+'\t'+str(psite_aa)+'\t'+str(esite_aa)+'\t'+str(pvalue)+'\t'+str(perm_cooperativity)+'\t'+str(cooperativity_pvalue)+'\t'+str(pos_cooperativity_pvalue)+'\t'+str(neg_cooperativity_pvalue)+'\n')
#stop = timeit.default_timer()
#total_time = stop-start
#pvalue_save.write(str(total_time))
#print(files)
#print('done')