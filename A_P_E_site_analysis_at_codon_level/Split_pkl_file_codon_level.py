import os
import pickle
dataset = ['Williams','Young','Weinberg','Nissley1','Nissley2']
path= '/gpfs/group/epo2/default/nks5512/find_cooperativity_and_no_cooperativity/A_P_Site_matrix_codon_level_find_cooperativity_and_no_cooperativity_with_75%_coverage_test/'

for data in dataset:
    outpath = path+'pkl_file/'+data
    count = 0

    try:
        os.makedirs(outpath)
    except OSError:
        print("Creation of the directory %s failed")
    else:
        print("Successfully created the directory %s ")
    
    PERMUTATION_LIST_file = open(path+'dic/A-site_profiles_'+data+'_LIST_FOR_PERMUTATION_dic.p','rb')
    PERMUTATION_LIST_load = pickle.load(PERMUTATION_LIST_file)
    
    for key in PERMUTATION_LIST_load.keys():
        outputfile=open(path+'pkl_file/'+data+'/A-site_profiles_'+data+'_LIST_FOR_PERMUTATION_dic_'+key+'.p','wb')
        'A-site_profiles_Jan_LIST_FOR_PERMUTATION_dic_'+key+'.p'
        pickle.dump(PERMUTATION_LIST_load[key],outputfile,protocol=pickle.HIGHEST_PROTOCOL)
        outputfile.close()
        print(count)
        count+=1
    PERMUTATION_LIST_file.close()
    del PERMUTATION_LIST_load

    
print('Done')