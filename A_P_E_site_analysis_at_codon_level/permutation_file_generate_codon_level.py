path= '/gpfs/group/epo2/default/nks5512/find_cooperativity_and_no_cooperativity/A_P_Site_matrix_codon_level_find_cooperativity_and_no_cooperativity_with_75%_coverage_test/'

dataset = ['Jan','Williams','Young','Weinberg','Nissley1','Nissley2','Pool']
#dataset = ['Jan']

#Intronic gene list
Intronic_gene_list = 'data/intronic_genes.p'

#Overlap gene list
Overlap_gene_list = 'data/overlap_genes.p'

# convert gene nucleotide to codon for every gene
codon_type_dic = 'data/codon_type_dict.p'


# 358 gene list to see robustness in data, This 358 gene data will be use to generate final amino acid metrix and tRNA matrix
Genelist_Williams_data = 'data/Genelist_Williams_data.tab'

# Threshold value for data (eg. Threshold value 4 for 6 data)
threshold_value = 4

#gene coverage thrushold
coverage_value = 75

print('\n Done')


import numpy as np
import math
import os
import pickle
print('Imported numpy, math, os and cPickle')
from optparse import OptionParser
print('Imported optionparser, plotter and misc')
import matplotlib.backends.backend_pdf as pdf
from scipy import stats
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.table import Table
import sys
import statsmodels.sandbox.stats.multicomp as mc
import operator as op
from matplotlib.ticker import FormatStrFormatter
from time import localtime, strftime
import itertools
from matplotlib import rcParams
import shutil
import matplotlib as mpl
mpl.rc('figure', max_open_warning = 0)
import statistics
import math
import random
import timeit
import multiprocessing
from multiprocessing import Pool
import subprocess


CODON_TYPES = ['UUU', 'UUC', 'UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG', 'AUU', 'AUC', 'AUA', 'AUG', 'GUU', 'GUC', 'GUA',
               'GUG', 'UCU', 'UCC', 'UCA', 'UCG', 'CCU', 'CCC', 'CCA', 'CCG', 'ACU', 'ACC', 'ACA', 'ACG', 'GCU', 'GCC',
               'GCA', 'GCG', 'UAU', 'UAC', 'CAU', 'CAC', 'CAA', 'CAG', 'AAU', 'AAC', 'AAA', 'AAG', 'GAU', 'GAC', 'GAA',
               'GAG', 'UGU', 'UGC', 'UGG', 'CGU', 'CGC', 'CGA', 'CGG', 'AGU', 'AGC', 'AGA', 'AGG', 'GGU', 'GGC', 'GGA',
               'GGG', 'UAA', 'UAG', 'UGA']

genetic_code = {'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C', 'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',
                'UUA': 'L', 'UCA': 'S', 'UAA': '*', 'UGA': '*', 'UUG': 'L', 'UCG': 'S', 'UAG': '*', 'UGG': 'W',
                'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R', 'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', 'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S', 'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
                'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G', 'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}

# In the following dict, synonymous codons for each amino acid are grouped in list such that they are decoded by similar tRNA.
# For example, in amino acid 'A', GCU and GCC codons are decoded by one type of tRNA while GCA and GCG are decoded by another kind of tRNA
synonymous = {'A': [['GCU', 'GCC'], ['GCA', 'GCG']],
              'C': [['UGU', 'UGC']],
              'D': [['GAU', 'GAC']],
              'E': [['GAA'], ['GAG']],
              'F': [['UUU', 'UUC']],
              'G': [['GGU', 'GGC'], ['GGA'], ['GGG']],
              'H': [['CAU', 'CAC']],
              'I': [['AUU', 'AUC'], ['AUA']],
              'K': [['AAG'], ['AAA']],
              'L': [['UUG'], ['UUA'], ['CUC', 'CUU'], ['CUA', 'CUG']],
              'M': [['AUG']],
              'N': [['AAU', 'AAC']],
              'P': [['CCA', 'CCG'], ['CCU', 'CCC']],
              'Q': [['CAA'], ['CAG']],
              'R': [['AGA'], ['CGU', 'CGC'], ['CGG', 'CGA'], ['AGG']],
              'S': [['UCU', 'UCC'], ['AGU', 'AGC'], ['UCA'], ['UCG']],
              'T': [['ACU', 'ACC'], ['ACA'], ['ACG']],
              'V': [['GUU', 'GUC'], ['GUG'], ['GUA']],
              'W': [['UGG']],
              'Y': [['UAU', 'UAC']],
              '*': [['UAA', 'UAG', 'UGA']]
              }

properties = {'A': 'NON-POLAR',
              'C': 'NON-POLAR',
              'D': 'Negative',
              'E': 'Negative',
              'F': 'NON-POLAR',
              'G': 'NON-POLAR',
              'H': 'Positive',
              'I': 'NON-POLAR',
              'K': 'Positive',
              'L': 'NON-POLAR',
              'M': 'NON-POLAR',
              'N': 'POLAR',
              'P': 'NON-POLAR',
              'Q': 'POLAR',
              'R': 'Positive',
              'S': 'POLAR',
              'T': 'POLAR',
              'V': 'NON-POLAR',
              'W': 'NON-POLAR',
              'Y': 'POLAR',
              }

hydropathy = {'A': '1.8',
              'C': '2.5',
              'D': '-3.5',
              'E': '-3.5',
              'F': '2.8',
              'G': '-0.4',
              'H': '-3.2',
              'I': '4.5',
              'K': '-3.9',
              'L': '3.8',
              'M': '1.9',
              'N': '-3.5',
              'P': '-1.6',
              'Q': '-3.5',
              'R': '-4.5',
              'S': '-0.8',
              'T': '-0.7',
              'V': '4.2',
              'W': '-0.9',
              'Y': '-1.3',
              }

AMINO_ACIDS = ['A', 'R', 'D', 'N', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']

tRNA_list = ['A-t1', 'A-t2', 'R-t1', 'R-t2', 'R-t3', 'R-t4', 'D-t1', 'N-t1', 'C-t1', 'E-t1', 'E-t2', 'Q-t1', 'Q-t2', 'G-t1', 'G-t2', 'G-t3', 'H-t1', 'I-t1', 'I-t2', 'L-t1', 'L-t2', 'L-t3', 'L-t4', 'K-t1', 'K-t2', 'M-t1', 'F-t1', 'P-t1', 'P-t2', 'S-t1', 'S-t2', 'S-t3', 'S-t4', 'T-t1', 'T-t2', 'T-t3', 'W-t1', 'Y-t1', 'V-t1', 'V-t2', 'V-t3', '*-t1']

base_pairing = {'A': {'Wobble': ['GCC', 'GCG'], 'Watson-Crick': ['GCA', 'GCU']}, 'R': {'Wobble': ['CGA', 'CGC'], 'Watson-Crick': ['AGA', 'AGG', 'CGG', 'CGU']},
                'D': {'Wobble': ['GAU'], 'Watson-Crick': ['GAC']},
                'N': {'Wobble': ['AAU'], 'Watson-Crick': ['AAC']}, 'C': {'Wobble': ['UGU'], 'Watson-Crick': ['UGC']}, 'E': {'Wobble': [], 'Watson-Crick': ['GAA', 'GAG']},
                'Q': {'Wobble': [], 'Watson-Crick': ['CAA', 'CAG']}, 'G': {'Wobble': ['GGU'], 'Watson-Crick': ['GGA', 'GGC', 'GGG']},
                'H': {'Wobble': ['CAU'], 'Watson-Crick': ['CAC']},
                'I': {'Wobble': ['AUC'], 'Watson-Crick': ['AUA', 'AUU']}, 'L': {'Wobble': ['CUG', 'CUU'], 'Watson-Crick': ['CUA', 'CUC', 'UUA', 'UUG']},
                'K': {'Wobble': [], 'Watson-Crick': ['AAA', 'AAG']},
                'M': {'Wobble': [], 'Watson-Crick': ['AUG']}, 'F': {'Wobble': ['UUU'], 'Watson-Crick': ['UUC']}, 'P': {'Wobble': ['CCC', 'CCG'], 'Watson-Crick': ['CCA', 'CCU']},
                'S': {'Wobble': ['UCC', 'AGU'], 'Watson-Crick': ['UCA', 'UCG', 'UCU', 'AGC']}, 'T': {'Wobble': ['ACC'], 'Watson-Crick': ['ACA', 'ACU', 'ACG']},
                'W': {'Wobble': [], 'Watson-Crick': ['UGG']}, 'Y': {'Wobble': ['UAU'], 'Watson-Crick': ['UAC']}, 'V': {'Wobble': ['GUC'], 'Watson-Crick': ['GUA', 'GUG', 'GUU']},
                '*': {'Wobble': [], 'Watson-Crick': ['UAA', 'UAG', 'UGA']}}

# Optimal codons selected based on their corresponding tRNA abundance (measured by RNA-Seq in Weinberg et al). Wobble only pairs are measured by 0.64*cognate tRNA concentration.
# Corrected mistake for G. Earlier it was 'G': {'Non-optimal': ['GGC', 'GGG'], 'Optimal': ['GGA', 'GGU']},
optimal_codon_usage = {'A': {'Non-optimal': ['GCC', 'GCG'], 'Optimal': ['GCA', 'GCU']}, 'C': {'Non-optimal': ['UGU'], 'Optimal': ['UGC']}, 'D': {'Non-optimal': ['GAU'], 'Optimal': ['GAC']},
                       'E': {'Non-optimal': ['GAG'], 'Optimal': ['GAA']}, 'F': {'Non-optimal': ['UUU'], 'Optimal': ['UUC']}, 'G': {'Non-optimal': ['GGA', 'GGG'], 'Optimal': ['GGC', 'GGU']},
                       'H': {'Non-optimal': ['CAU'], 'Optimal': ['CAC']}, 'I': {'Non-optimal': ['AUA'], 'Optimal': ['AUC', 'AUU']}, 'K': {'Non-optimal': ['AAA'], 'Optimal': ['AAG']},
                       'L': {'Non-optimal': ['CUA', 'CUC', 'CUG', 'CUU'], 'Optimal': ['UUA', 'UUG']}, 'M': {'Non-optimal': [], 'Optimal': ['AUG']}, 'N': {'Non-optimal': ['AAU'], 'Optimal': ['AAC']},
                       'P': {'Non-optimal': ['CCC', 'CCU'], 'Optimal': ['CCA', 'CCG']}, 'Q': {'Non-optimal': ['CAG'], 'Optimal': ['CAA']},
                       'R': {'Non-optimal': ['AGG', 'CGG', 'CGA', 'CGC'], 'Optimal': ['AGA', 'CGU']}, 'S': {'Non-optimal': ['UCA', 'UCG', 'AGU', 'AGC'], 'Optimal': ['UCC', 'UCU']},
                       'T': {'Non-optimal': ['ACA', 'ACG'], 'Optimal': ['ACC', 'ACU']}, 'V': {'Non-optimal': ['GUA', 'GUG'], 'Optimal': ['GUC', 'GUU']}, 'W': {'Non-optimal': [], 'Optimal': ['UGG']},
                       'Y': {'Non-optimal': ['UAU'], 'Optimal': ['UAC']}}

# Most optimal codon for every amino acid
most_optimal_codon = {'A': 'GCU', 'C': 'UGC', 'D': 'GAC', 'E': 'GAA', 'F': 'UUC', 'G': 'GGC', 'H': 'CAC', 'I': 'AUU', 'K': 'AAG', 'L': 'UUG', 'M': 'AUG', 'N': 'AAC', 'P': 'CCA', 'Q': 'CAA', 'R': 'AGA',
                      'S': 'UCU', 'T': 'ACU', 'V': 'GUU', 'W': 'UGG', 'Y': 'UAC', '*': 'UAA'}

# Optimal and non-optimal codons based on Penchman, Frydman, tAI cutoff of 0.47 as well as used for codon optimality in Jeff Coller's paper.
optimal_dict = {'Optimal': ['GCU', 'GCC', 'GAC', 'GAA', 'UUC', 'GGC', 'AUU', 'AUC', 'AAG', 'UUG', 'AUG', 'AAC', 'CCA', 'CAA', 'AGA', 'UCU', 'UCC', 'ACU', 'ACC', 'GUU', 'GUC', 'UAC'],
                'Non-optimal': ['GCA', 'GCG', 'UGC', 'UGU', 'GAU', 'GAG', 'UUU', 'GGU', 'GGA', 'GGG', 'CAC', 'CAU', 'AUA', 'AAA', 'UUA', 'CUA', 'CUC', 'CUG', 'CUU', 'AAU', 'CCG', 'CCU', 'CCC',
                                'CAG', 'CGU', 'AGG', 'CGC', 'CGG', 'CGA', 'UCA', 'AGC', 'UCG', 'AGU', 'ACA', 'ACG', 'GUG', 'GUA', 'UGG', 'UAU']}

CHROMOSOMES = ['chrI', 'chrII', 'chrIII', 'chrIV', 'chrV', 'chrVI', 'chrVII', 'chrVIII', 'chrIX', 'chrX', 'chrXI', 'chrXII', 'chrXIII', 'chrXIV', 'chrXV', 'chrXVI', 'chrM']

print ('\nDone')













def vector(a,b):
    a_sq = a**2
    b_sq = b**2
    ab_total = a_sq+b_sq
    vector_of_ab = math.sqrt(ab_total)

    return vector_of_ab


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
        pickle.dump(LIST_FOR_PERMUTATION, open(path+pkl_file+filename+'_LIST_FOR_PERMUTATION_dic_'+asite_aa+'.p', 'wb'))
    pickle.dump(COOPRATIVITY_LIST, open(path+filename+'_COOPRATIVITY_LIST.p', 'wb'))
    return




data_name = sys.argv[1]


from pexecute.process import ProcessLoom
loom = ProcessLoom(max_runner_cap=8000)
import multiprocessing 
import time
from sklearn.utils import resample
#Clculate cooprativity and no cooprativity effect
input1 = path+'A_site_profile_file/esite_asite_matrix_for_codon_level/'
print(input1)
files = os.listdir(input1)
print(files)
jobs = []

outpath = path+'A_site_profile_file/cooperativity_and_cooperativity_effect/'

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
