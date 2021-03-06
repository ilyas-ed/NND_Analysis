# NND_Analysis

Code by Ilyas Ed-daoui

Data by Álvaro Castells-García

Python version 3.7


These Python scripts can be used to perform the NND analysis. To generate the NNDs for different molecules, the scripts use the centroids obtained from the cluster analysis that has been published in this paper:


‘ Ricci MA, Manzo C, García-Parajo MF, Lakadamyali M, Cosma MP. Chromatin fibers are formed by heterogeneous groups of nucleosomes in vivo. Cell. 2015 Mar 12;160(6):1145-58. doi: 10.1016/j.cell.2015.01.054. PMID: 25768910. ‘

---------------------------------------------------------------------------------------------

NND-clusterArea.py : This scripts finds cluster pairs between clusters of two different types. The cluster pairs are then binned by cluster size and its frequency is represented in a heat map. Values are normalized to sum to unit.


Input: two .mat file with the characteristics of each cluster type, including its number of localizations per cluster and its area. It is the output of the cluster analysis.


Output: heat map figure of the cluster pair distribution


Data: in the folder “RNA-Protein NND example Dataset” there are the following files BJFIB_RNA.mat: .mat file with the output of the Carlo Manzo cluster analysis file. This should be one of the inputs for the python script 
BJFIB_RNAPII.mat: .mat file with the output of the Carlo Manzo cluster analysis file. This should be the other input for the python script.






---------------------------------------------------------------------------------------------

NND-GAPDH: This scripts finds cluster pairs between clusters of two different types. It was used to calculate the size of the nucleoplasmic RNA clutches.

Input: two .mat file with the characteristics of each cluster type, including its number of localizations per cluster and its area. It is the output of the cluster analysis.

Output: a .csv file with the characteristics of the clusters comprising each cluster pair.

Data: in the folder “NND-GAPDH example dataset” there are the following files 

FISH_clusters_cellXXX: .mat file with the output of the Carlo Manzo cluster analysis file, used to analyze GAPDH smRNA-FISH signal on the cytoplasm. RNA_clusters_cellXXX: .mat file with the output of the Carlo Manzo cluster analysis file, used to analyze nascent RNA super resolved signal on the cytoplasm.



---------------------------------------------------------------------------------------------

MUC4-NND: This scripts uses as a basis for the NND the MUC4 locus, equivalent to a single, well defined cluster. It calculates the distances between the MUC4 locus and its 10 nearest neighbor RNA nanodomains. Once the NNDs are collected for all the cells, the mean and the standard deviation are calculated between the distances in the same closeness rank.


Input: 

.txt file: matrix in which the x, y position are in the columns 4 and 5, respectively. It is used to calculate the arithmetic center of the MUC4 locus 

.mat file: a .mat file with the characteristics of RNA nanodomains. It the output from the cluster analysis


Output: figure with the NND ranks distances and standard deviation


Data: in the folder “Locus-NND example dataset” there are the following files


GP220_+IL6_MUC4_20minEU.txt / GP220_NOIL6_MUC4_20minEU.txt: .txt file with the original localization list of the MUC4 locus in .txt format, for cells of two different conditions. This will be one of the input files for the script 

GP220_+IL6_RNA_20minEU.mat / GP220_NOIL6_RNA_20minEU.mat: .mat file with the output from the cluster analysis of the super resolved RNA in .mat format, for cells of two different conditions.
