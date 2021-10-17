# A-P-and-E-site-analysis
Translation elongation rate analysis at amino acid, tRNA, and Codon level by ribosome profiling dataset

Here 3 diffrent analysis script generate to analysis of amino acid, tRNA and codon level.
Every Jupyter notebook split in 3 steps. 
1. Read raw reads density file than calculate cooperativity effect and p-value through permutation effects all 7 diffrent dataset separately.
2. Use cooperativity effect and p-value Calculated file to calculate robustness of triplets.
3. Robust file use for further visualtion of result. 

## Important Note
Before using tRNA and Codon level Jupyter notebook. Run the Unzip_file_tRNA_level.py script for tRNA and Unzip_file_codon_level.py for codon level to unzip some files. 

# Point of action cover in amino acid level analysis
1. Cooperativity effect calculated for E and P-site effect on A-site
2. 7 ribo-seq data used for this analysis
3. Normalised ribsome density calculated for each gene
4. Only 75% coverage gene selected for analysis
5. Statistically significance calculate for each triplet at E-, P- and A-site by permutation test
6. Speed define by delta1
7. Cooperativity effect calculated through delta1, delta2, and delta3
8. Corrected p-value calculated for each triplet by Benjamini Hochberg correction method to control false discovery
9. Only top 70% triplets selected based on instances of each triplet.
10. The robustness of each triplet defines by the three filter criteria. For each triplet, a minimum of one dataset should be statistically significant (check corrected p-value). All seven dataset delta1 directions should be in the same direction. Any four datasets should show the same cooperativity effect.
11. analysis of calculated data
12. RiBO-seq data correlation at codon label
13. Create bar plot to show E-site effect (Average of delta1 medium)
14. A, P and E site matrix (Any single site fix)
15. Compare A, P & A, P and E site, and present how E site play import role in regulating the speed of translation
16. Speed change (A and P site fix and E site change), triplet selected only from Robust file (Total 588 robust triplet)
17. Plotting the normalized ribosome density distributions between two triplets (triplet pair example show how mutated at E-site change the translation elongation rate)
18. Identify E-site and P-site effect triplets from robust triplets
19. Plotting the normalized ribosome density distributions between E-site and P-site selected triplets
20. Switching the P and E site amino acid position resulting get opposite speed changes. We have selected the triplet pairs from robust triplet for this analysis. Plotting the normalized ribosome density distributions between selected triplets.
21. Mutating at the E- and P-site amino acid is predicted to alter the translation elongation rate at the A-site. Compare all 8000 amino acid triplet with each other in complete proteome and calculate odds, p-value, and percentage difference
22. Plotting the normalized ribosome density distributions between E-site and P-site selected triplets. Best triplet pair selected based on best odds value of each amino acid at A-site. In total, between 20 triplet pairs, normalized ribosome density distributions were plotted. 
23. Signatures of evolutionary selection for fast and slow translating amino acid triplet
24. Function definitions for the Linker vs Domain analysis. Demonstrates that slow pairs are enriched in linker regions with respect to domain regions
25. Compare AP and APE site, and present how E site play import role in regulating the speed of translation
26. The combined average effect of a molecular factor on translation elongation rate at the A, P and E-site
27. Find out predicted mutant triplet pair significance. We show that mutational experiments are consistent with amino acid identified influencing translation elongation rates



# Point of action cover at tRNA level analysis
1. Cooperativity effect calculated for E and P-site effect on A-site
2. 7 ribo-seq data used for thsis analysis
3. Normalised ribsome density calculated for each gene
4. Only 75% coverage gene selected for analysis
5. Statistically significance calculate for each triplet at E-, P- and A-site by permutation test
6. Speed define by delta1
7. Cooperativity effect calculated through delta1, delta2, and delta3
8. Corrected p-value calculated for each triplet by Benjamini Hochberg correction method to control false discovery
10. The robustness of each triplet defines by the two filter criteria. For each triplet, All seven dataset delta1 directions should be in the same direction. Any four datasets should show the same cooperativity effect.
11. analysis of calculated data
12. The combined average effect of a molecular factor on translation elongation rate at the A, P and E-site


# Point of action cover at codon level analysis
1. Cooperativity effect calculated for E and P-site effect on A-site
2. 7 ribo-seq data used for thsis analysis
3. Normalised ribsome density calculated for each gene
4. Only 75% coverage gene selected for analysis
5. Statistically significance calculate for each triplet at E-, P- and A-site by permutation test
6. Speed define by delta1
7. Cooperativity effect calculated through delta1, delta2, and delta3
8. Corrected p-value calculated for each triplet by Benjamini Hochberg correction method to control false discovery
10. The robustness of each triplet defines by the two filter criteria. For each triplet, All seven dataset delta1 directions should be in the same direction. Any four datasets should show the same cooperativity effect.
11. analysis of calculated data
12. The combined average effect of a molecular factor on translation elongation rate at the A, P and E-site


## Installation
Start by cloning this repository using git:
```
git clone https://github.com/obrien-lab/A-P-and-E-site-analysis
cd A-P-and-E-site-analysis
```

Create a Python 3.8 conda environment to manage the required software packages:
```
conda create --name py38 python=3.8 -y requirement.txt
```

Activate the environment:
```
source activate py38
```

Begin installing packages with conda and pip:
```
conda install -n <env_name> requirements.txt
or install manualy
conda install -y -q 'Package name'
```

## Run the Jupyter notebooks
Activate the py38 environment and run Jupyter Notebook.  Within Jupyter Notebook, navigate to one of the dataset folders (A_P_E_site_analysis_at_amino_acid_level, A_P_E_site_analysis_at_tRNA_level, A_P_E_site_analysis_at_codon_level), and open the .ipynb file.
```
source activate py38
jupyter notebook
```
