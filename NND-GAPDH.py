

'''

Written by Ilyas Ed-daoui, ilyas@gibh.ac.cn

'''

'''

This script uses as input the localizations of the Muc4 and the centroids obtained 
from the cluster analysis that has been published in this paper:
‘
Ricci MA, Manzo C, García-Parajo MF, Lakadamyali M, Cosma MP. 
Chromatin fibers are formed by heterogeneous groups of nucleosomes in vivo. Cell. 
2015 Mar 12;160(6):1145-58. doi: 10.1016/j.cell.2015.01.054. PMID: 25768910.
‘

This scripts finds cluster pairs between clusters of two different types. It was used to calculate the 
size of the nucleoplasmic RNA clutches. 

'''



import pandas as pd
from tqdm import tqdm
import glob
import os
from os import path
import scipy.io as sio
import math
from sklearn.neighbors import BallTree


#! **********     Global Parameters      #! **********

folder_path = ''

#! **********     Global Parameters      #! **********


subfolders_paths = [f.path for f in os.scandir(folder_path) if f.is_dir()]

for sub_path in tqdm(subfolders_paths):
    subsubfolders = [f.path for f in os.scandir(sub_path) if f.is_dir()]

    for sub in subsubfolders:
        fish_path = glob.glob(sub + "/FISH.mat")
        if fish_path:
            mat = sio.loadmat(fish_path[0])
            fish_xydata1 = mat['ClusterResults']
            fish_df = pd.DataFrame(fish_xydata1).iloc[:, 0:3]


        rna_path = glob.glob(sub + "/RNA.mat")
        if rna_path:
            mat = sio.loadmat(rna_path[0])
            rna_xydata1 = mat['ClusterResults']   
            rna_df = pd.DataFrame(rna_xydata1).iloc[:, 0:3] 


        fish_ = pd.DataFrame({
            'id_fish': k,
            'x_fish': e,
            'y_fish': f,
            'numberLocs_fish': g
        } for e, f, g, k in zip(fish_df.iloc[:, 0].values, fish_df.iloc[:, 1].values, fish_df.iloc[:, 2].values, range(len(fish_df)))
        )


        rna_ = pd.DataFrame({
            'id_rna': k,
            'x_rna': e,
            'y_rna': f,
            'numberLocs_rna': g
        } for e, f, g, k in zip(rna_df.iloc[:, 0].values, rna_df.iloc[:, 1].values, rna_df.iloc[:, 2].values, range(len(rna_df)))
        )


        tree = BallTree(rna_[['x_rna', 'y_rna']].values, leaf_size=2)
        fish_['distance'], fish_['id_rna'] = tree.query(fish_[['x_fish', 'y_fish']].values, k=1)
        merged_rna_fish = rna_.merge(fish_, on='id_rna')
        # print(merged_rna_fish)

        merged_rna_fish.to_csv(sub+"/rna_fish_nnd.csv")

   

