path= '/gpfs/group/epo2/default/nks5512/find_cooperativity_and_no_cooperativity/A_P_Site_matrix_Amino_acid_level_find_cooperativity_and_no_cooperativity_with_75%_coverage/'


AMINO_ACIDS = ['A', 'R', 'D', 'N', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']

import numpy as np
import math
import os
#import cPickle as Pickle
import pickle
print('Imported numpy, math, os and cPickle')
from optparse import OptionParser
#import plotter
print('Imported optionparser, plotter and misc')
from scipy import stats

import pandas as pd
import sys




def odds_speed_change_aa(densities_list1, densities_list2):
    list_diff = []
    for density1 in densities_list1:
        for density2 in densities_list2:
            diff = density2-density1
            list_diff.append(diff)
    median_diff = np.median(densities_list2) - np.median(densities_list1)
    # If the median of second list is higher, then  Density_list1 is fast and Density_list2 is slow and hence we calculate the odds of slowdown upon mutation from 1-->2
    if median_diff > 0:
        try:
            odds = float(sum(i > 0 for i in list_diff))/float(sum(i < 0 for i in list_diff))
        except ZeroDivisionError:
            odds = -1
    # If the median of second list is lower, then  Density_list1 is slow and Density_list2 is fast and hence we calculate the odds of speedup upon mutation from 1-->2
    else:
        try:
            odds = float(sum(i < 0 for i in list_diff))/float(sum(i > 0 for i in list_diff))
        except ZeroDivisionError:
            odds = -1

    return odds


def Identify_significant_triplet(input_file):    
    file_read = open(input_file, 'r')
    next(file_read)
    significnat_pair_list = []
    for line in file_read:
        line1 = line.split('\t')
        if line1[4] == 'Non_Robust':
            continue
        significnat_pair_list.append(line1[0]+line1[1]+line1[2])

    return significnat_pair_list

def Total_triplets_in_robust_file(input_file):    
    file_read = open(input_file, 'r')
    next(file_read)
    significnat_pair_list = []
    for line in file_read:
        line1 = line.split('\t')
        #if line1[4] == 'Non_Robust':
            #continue
        significnat_pair_list.append(line1[0]+line1[1]+line1[2])

    return significnat_pair_list

def Total_triplets():
    pair_list = []
    for aa in AMINO_ACIDS:
        if aa == '*':
            continue
        for psite_aa in AMINO_ACIDS:            
            if psite_aa == '*':
                continue            
            for esite_aa in AMINO_ACIDS:                
                if esite_aa == '*':
                    continue
                pair_list.append(aa+psite_aa+esite_aa)

    return pair_list


def directory_create(out_put_path):
    
    outpath = out_put_path
    try:
        os.makedirs(outpath)
    except OSError:
        print("Creation of the directory %s failed")
    else:
        print("Successfully created the directory %s ")
        
    return



def combination_dic(out_path, compare_triplets_list):    
    triplet_compare_dic = {}
    triplet_compare_result_dic = {}
    compare_triplets_list1 = compare_triplets_list
    triplet_list = []
    
    for triplet in compare_triplets_list:
        triplet_list.append(triplet)
        triplet_compare_dic[triplet] = []
        triplet_compare_result_dic[triplet] = {}
        for triplet1 in compare_triplets_list1:
            if triplet == triplet1:
                continue
            triplet_compare_dic[triplet].append(triplet1)
            triplet_compare_result_dic[triplet][triplet1] = []
        compare_triplets_list1 = compare_triplets_list1[1:]
    print(triplet_compare_dic)
    print(triplet_compare_result_dic)
    pickle.dump(triplet_compare_dic, open(out_path + 'triplets_dic.p', 'wb'))
    pickle.dump(triplet_compare_result_dic, open(out_path + 'result_triplets_dic.p', 'wb'))
    
    return triplet_list



def list_sort_by_number_of_instances(triplet, out_path, Input_file_path_pickel, input_file_for_significant_triplet):
    robust_file = pd.read_csv(input_file_for_significant_triplet, sep = '\t')
    dict_triplet = pickle.load(open(out_path+'triplets_dic.p', 'rb'))
    ribsome_density_dictionary = pickle.load(open(Input_file_path_pickel, 'rb'))
    dict_triplet_list = dict_triplet[triplet]
    dict_triplet_list.append(triplet)
    
    # comepare triplet with all triplets and calculate odds, p-value and percentage change
    directory_create(out_path+'result/')
    out_path_result = out_path+'result/'
    result_output_file = open(out_path_result+'each_significant_triplet_ribosome_denity_compare_with_'+str(triplet)+'_odds, pvalue_and_perc_change.tab', 'w')
    result_output_file.write('Asite\tPsite\tEsite\talt_Asite\talt_Psite\talt_Esite\tpercentage_change_Speed\tdelta1_Speed\talt_delta1_Speed\tsize\talt_size\tp_value\tperc_change\todds\n')
    A_site = triplet[0]
    P_site = triplet[1]
    E_site = triplet[2]
    asite_list = robust_file.Asite == A_site
    psite_list = robust_file.Psite == P_site
    esite_list = robust_file.Esite == E_site
    speed_of_triplet = asite_list & psite_list & esite_list    
    try:
        speed_status1 = robust_file[speed_of_triplet].Speed.values.tolist()[0]
    except IndexError:
        speed_status1 = 'None'
    
    ribosome_density_list = ribsome_density_dictionary[triplet][0]

    size = len(ribosome_density_list)

    for key in dict_triplet_list:
        if key == triplet:
            continue
        alt_ribosome_density_list = ribsome_density_dictionary[key][0]
        alt_size = len(alt_ribosome_density_list)
        A_site = key[0]
        P_site = key[1]
        E_site = key[2]
        asite_list = robust_file.Asite == A_site
        psite_list = robust_file.Psite == P_site
        esite_list = robust_file.Esite == E_site
        speed_of_triplet = asite_list & psite_list & esite_list    
        try:
            speed_status2 = robust_file[speed_of_triplet].Speed.values.tolist()[0]
        except IndexError:
            speed_status2 = 'None'
        
        u, p = stats.mannwhitneyu(ribosome_density_list, alt_ribosome_density_list)
        #perc_change = ((np.median(ribosome_density_list) - np.median(alt_ribosome_density_list)) / np.median(alt_ribosome_density_list)) * 100
        perc_change = (np.median(ribosome_density_list) - np.median(alt_ribosome_density_list)) / ((np.median(ribosome_density_list)+np.median(alt_ribosome_density_list))/2)* 100
        
        if perc_change > 0:
            speed_status = 'Slow_Down'
        if perc_change < 0:
            speed_status = 'Speed_up'
        
        odds = odds_speed_change_aa(ribosome_density_list, alt_ribosome_density_list)
        result_output_file.write(triplet[0]+'\t'+triplet[1]+'\t'+triplet[2]+'\t'+key[0]+'\t'+key[1]+'\t'+key[2]+'\t'+speed_status+'\t'+str(speed_status1)+'\t'+str(speed_status2)+'\t'+str(size)+'\t'+str(alt_size)+'\t'+str(p)+'\t'+str(perc_change)+'\t'+str(odds)+'\n')
   
    return






triplet = sys.argv[1]
print(triplet)
cwd = os.getcwd()

out_path = cwd+'/Result/Analysis/Mutational_prediction_test/'
#out_path = path+'A_site_profile_file/matrix_file/Test4_compare_triplets_with_other_triplets/'
Input_file_path_pickel = cwd+'/Result/esite_asite_matrix_for_amino_acid/A-site_profiles_Pool_Translation_timesssAsite_psite_esite_times_coordinates.p'
input_file_for_significant_triplet = cwd+'/Result/robust_amino_acid_matrix/Robust_aminoacid_pairs_Not_controlled_for_any_factors_amino_acid_matrix_pvalue_1_delta1_direction_7_cooperativity_effect_direction_threshold_4_datasets_100%.tab'

directory_create(out_path)

list_sort_by_number_of_instances(triplet, out_path, Input_file_path_pickel, input_file_for_significant_triplet)
print('Done')