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
import time

def vector(a,b):
    a_sq = a**2
    b_sq = b**2
    ab_total = a_sq+b_sq
    vector_of_ab = math.sqrt(ab_total)

    return vector_of_ab

#this is a temporary work around
def returnshuf(x):
    new_x = list(x)  # make a copy
    np.random.shuffle(x)
    return(new_x)


def permutation_test(a_p_pair_list_for_permutation,a_e_pair_list_for_permutation,a_pair_list_for_permutation,a_p_e_pair_list_for_permutation,delta,asite_aa, psite_aa, esite_aa):

    # (X,Y,Z) are in the E-, P- and A-sites.
    # (x,..,z) = a_e_pair
    # (..,y,z) = a_p_pair
    # (..,..,z) = a_pair
    # (x,y,z) = a_p_e_pair
    cooperetivity_count = 0
    pos_cooperetivity_count = 0
    neg_cooperetivity_count = 0

    #length calculate
    # (x,..,z) = a_e_pair
    length_of_e_and_a_list = len(a_e_pair_list_for_permutation)
    # (..,y,z) = a_p_pair
    length_of_a_p_list = len(a_p_pair_list_for_permutation)
    # (..,..,z) = a_pair
    length_of_a_list = len(a_pair_list_for_permutation)
    # (x,y,z) = a_p_e_pair
    length_of_a_p_e_list = len(a_p_e_pair_list_for_permutation)


    # delta1 = median p(X,Y,Z) - median p(...,...,Z)
    # delta2 = median p(...,Y,Z) - median p(...,...,Z)
    # delta3 = median p(X,...,Z) - median p(...,...,Z)
        
    #all 4 list shuffle and create new list
    combined_list = a_p_e_pair_list_for_permutation + a_p_pair_list_for_permutation + a_e_pair_list_for_permutation + a_pair_list_for_permutation
    list_len = len(combined_list)


    #*#*#*#start for loop replacement using numpy 'vectorization'
    nloop=10000 #chunk if required by mem
    listofshufs=[returnshuf(combined_list) for x in range(nloop)]
    medianvec_new_a_pair_shuffle_list = np.median([x[(length_of_a_p_list+length_of_a_p_e_list+length_of_e_and_a_list):] for x in listofshufs],axis=1)
    delta1_perm_vec = np.median([x[:length_of_a_p_e_list] for x in listofshufs], axis=1) - medianvec_new_a_pair_shuffle_list
    delta2_perm_vec = np.median([x[length_of_a_p_e_list:(length_of_a_p_list+length_of_a_p_e_list)] for x in listofshufs], axis=1) - medianvec_new_a_pair_shuffle_list
    delta3_perm_vec = np.median([x[(length_of_a_p_list+length_of_a_p_e_list):(length_of_a_p_list+length_of_a_p_e_list+length_of_e_and_a_list)] for x in listofshufs], axis=1) - medianvec_new_a_pair_shuffle_list
    delta_perm = delta1_perm_vec-delta2_perm_vec-delta2_perm_vec

    for element in delta_perm:
        if element > delta:
            pos_cooperetivity_count+=1    
        if element < delta:
            neg_cooperetivity_count+=1
        if abs(element) > abs(delta):
            cooperetivity_count+=1
    #*#*#*end for loop replacement   


    cooperativity_pvalue = cooperetivity_count/10000
    pos_cooperativity_pvalue = pos_cooperetivity_count/10000
    neg_cooperativity_pvalue = neg_cooperetivity_count/10000
    
    

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
        
    
    return asite_aa, psite_aa, esite_aa, pvalue, perm_cooperativity, cooperativity_pvalue, pos_cooperativity_pvalue, neg_cooperativity_pvalue



data_name = sys.argv[1]

asite = sys.argv[2]
psite = sys.argv[3]
esite = sys.argv[4]

path1 = data_name+"/"

try:
	os.makedirs(path1)
except OSError:
	pass
	#print("Creation of the directory %s failed")
else:
	pass
	#print("Successfully created the directory %s ")

path='/gpfs/group/epo2/default/nks5512/find_cooperativity_and_no_cooperativity/A_P_Site_matrix_tRNA_level_find_cooperativity_and_no_cooperativity_with_75%_coverage_test/'
pvalue_save = open(path+path1+'/'+data_name+'_'+asite+'_'+psite+'_'+esite+'.tab','w')


picklepath=path+'pkl_file/'
tic=time.time()
PERMUTATION_LIST_file = open(picklepath+data_name+'/A-site_profiles_'+data_name+'_LIST_FOR_PERMUTATION_dic_'+asite+'.p','rb')
PERMUTATION_LIST_load = pickle.load(PERMUTATION_LIST_file)

a_p_pair_list_for_permutation = PERMUTATION_LIST_load[psite][esite][0]
a_e_pair_list_for_permutation = PERMUTATION_LIST_load[psite][esite][1]
a_pair_list_for_permutation = PERMUTATION_LIST_load[psite][esite][2]
a_p_e_pair_list_for_permutation = PERMUTATION_LIST_load[psite][esite][3]
v = PERMUTATION_LIST_load[psite][esite][4]
PERMUTATION_LIST_file.close()
del PERMUTATION_LIST_load
asite_aa, psite_aa, esite_aa, pvalue, perm_cooperativity, cooperativity_pvalue, pos_cooperativity_pvalue, neg_cooperativity_pvalue = permutation_test(a_p_pair_list_for_permutation, a_e_pair_list_for_permutation, a_pair_list_for_permutation, a_p_e_pair_list_for_permutation, v, asite, psite, esite)
print('For loop time:', time.time()-tic)

pvalue_save.write(str(asite_aa)+'\t'+str(psite_aa)+'\t'+str(esite_aa)+'\t'+str(pvalue)+'\t'+str(perm_cooperativity)+'\t'+str(cooperativity_pvalue)+'\t'+str(pos_cooperativity_pvalue)+'\t'+str(neg_cooperativity_pvalue)+'\n')

