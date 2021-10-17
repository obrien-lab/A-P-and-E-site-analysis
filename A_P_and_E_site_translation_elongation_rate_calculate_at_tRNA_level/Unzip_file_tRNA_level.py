# Unzip 7z file
# install pachage  - python -m pip install py7zr
import py7zr
import os

File1 = 'Permutation_test/cooprative_effect_result_dic_abs.7z'
with py7zr.SevenZipFile(File1, 'r') as archive:
    archive.extractall(path="Permutation_test/")
os.remove(File1)

File2 = 'Permutation_test/Permutation_test_result_dic.7z'
with py7zr.SevenZipFile(File2, 'r') as archive:
    archive.extractall(path="Permutation_test/")
os.remove(File2)

File3 = 'Result/robust_amino_acid_matrix/Robust_tRNA_pairs_Not_controlled_for_any_factors_amino_acid_matrix_0_datasets.7z'
with py7zr.SevenZipFile(File3, 'r') as archive:
    archive.extractall(path="Result/robust_amino_acid_matrix/")
os.remove(File3)

print('Unizip completed')