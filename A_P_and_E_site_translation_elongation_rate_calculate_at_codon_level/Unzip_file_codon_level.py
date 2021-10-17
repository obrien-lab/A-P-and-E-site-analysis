# Unzip 7z file
# install pachage  - python -m pip install py7zr
import py7zr

File1 = 'Permutation_test/cooprative_effect_result_dic_abs.7z'
with py7zr.SevenZipFile(File1, 'r') as archive:
    archive.extractall(path="Permutation_test/")

File2 = 'Permutation_test/Permutation_test_result_dic.7z'
with py7zr.SevenZipFile(File2, 'r') as archive:
    archive.extractall(path="Permutation_test/")
File3 = 'Result/robust_amino_acid_matrix/Merged_pair_statistics_Not_controlled_for_any_factors_amino_acid_matrix.7z'
with py7zr.SevenZipFile(File3, 'r') as archive:
    archive.extractall(path="Result/robust_amino_acid_matrix/")
File4 = 'Result/robust_amino_acid_matrix/Robust_aminoacid_pairs_Not_controlled_for_any_factors_amino_acid_matrix_0_datasets.7z'
with py7zr.SevenZipFile(File4, 'r') as archive:
    archive.extractall(path="Result/robust_amino_acid_matrix/")

print('Unizip completed')