
import os
import pickle

path= '/gpfs/group/epo2/default/nks5512/find_cooperativity_and_no_cooperativity/A_P_Site_matrix_tRNA_level_find_cooperativity_and_no_cooperativity_with_75%_coverage_test/'





dataset = ['Jan','Williams','Young','Weinberg','Nissley1','Nissley2','Pool']
for data in dataset:
    outpath = path+'pkl_file/'+data
    count = 0

    try:
        os.makedirs(outpath)
    except OSError:
        print("Creation of the directory %s failed")
    else:
        print("Successfully created the directory %s ")
    PERMUTATION_LIST_load = pickle.load(open(path+'A-site_profiles_'+data+'_LIST_FOR_PERMUTATION_dic.p','rb'))
    for key in PERMUTATION_LIST_load.keys():
        outputfile=open(path+'pkl_file/'+data+'/A-site_profiles_'+data+'_LIST_FOR_PERMUTATION_dic_'+key+'.p','wb')
        'A-site_profiles_Jan_LIST_FOR_PERMUTATION_dic_'+key+'.p'
        pickle.dump(PERMUTATION_LIST_load[key],outputfile,protocol=pickle.HIGHEST_PROTOCOL)
        outputfile.close()
    
print('Done')