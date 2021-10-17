
import os
import pickle

cwd = os.getcwd()
path = cwd+'/'

dataset = ['Jan','Williams','Young','Weinberg','Nissley1','Nissley2','Pool']
input_path= path+'Permutation_test/permutation_file_generate/'
for data in dataset:
    outpath = path+'Permutation_test/seperate_cooperativity_and_cooperativity_effect_pkl_file/'+data
    count = 0
    try:
        os.makedirs(outpath)
    except OSError:
        print("Creation of the directory %s failed")
    else:
        print("Successfully created the directory %s ")
    PERMUTATION_LIST_file = open(input_path+'A-site_profiles_'+data+'_LIST_FOR_PERMUTATION_dic.p','rb')
    PERMUTATION_LIST_load = pickle.load(PERMUTATION_LIST_file)
    for key in PERMUTATION_LIST_load.keys():
        outputfile=open(outpath+'/A-site_profiles_'+data+'_LIST_FOR_PERMUTATION_dic_'+key+'.p','wb')
        'A-site_profiles_Jan_LIST_FOR_PERMUTATION_dic_'+key+'.p'
        pickle.dump(PERMUTATION_LIST_load[key],outputfile,protocol=pickle.HIGHEST_PROTOCOL)
        outputfile.close()
        print(count)
        count+=1
    PERMUTATION_LIST_file.close()
    del PERMUTATION_LIST_load
    
print('Done')
