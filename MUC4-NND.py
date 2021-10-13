

'''

Written by Ilyas Ed-daoui, ilyas@gibh.ac.cn


This script uses as input the localizations of the Muc4 and the centroids obtained 
from the cluster analysis that has been published in this paper:
‘
Ricci MA, Manzo C, García-Parajo MF, Lakadamyali M, Cosma MP. 
Chromatin fibers are formed by heterogeneous groups of nucleosomes in vivo. Cell. 
2015 Mar 12;160(6):1145-58. doi: 10.1016/j.cell.2015.01.054. PMID: 25768910.
‘

It calculates the distances between the MUC4 locus and its 10 nearest neighbor RNA nanodomains. 
Once the NNDs are collected for all the cells, the mean and the standard deviation are calculated 
between the distances in the same closeness rank.

'''


from tqdm import tqdm

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm

import scipy.io as sio
import scipy
import math 

import seaborn as sns

#! **********     Global Parameters      #! **********

BINS=1
NUMBER_OF_CLUSTERS = 10 # This is the number of clusters we want to show their centers distances to the MUC4 centroids

#! **********     Global Parameters      #! **********



#*                                                                        
#*                Files reading                       
#*   

#!!!!   -IL6
#?                         
Muc4_localizations_1 = pd.read_csv() # Muc4 file path
mat1 = sio.loadmat() # mat file path
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_1 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_1.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_1.sort()
list_distances_1 = [i * 160 for i in list_distances_1]
#?

#?  
Muc4_localizations_1 = pd.read_csv()
mat1 = sio.loadmat()
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_2 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_2.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_2.sort()
list_distances_2 = [i * 160 for i in list_distances_2]
#?  

#?  
Muc4_localizations_1 = pd.read_csv()
mat1 = sio.loadmat()
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_3 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_3.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_3.sort()
list_distances_3 = [i * 160 for i in list_distances_3]
#?  

#?  
Muc4_localizations_1 = pd.read_csv()
mat1 = sio.loadmat()
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_4 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_4.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_4.sort()
list_distances_4 = [i * 160 for i in list_distances_4]
#?  

#?  
Muc4_localizations_1= pd.read_csv()
mat1 = sio.loadmat()
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_5 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_5.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_5.sort()
list_distances_5 = [i * 160 for i in list_distances_5]
#?  

#?  
Muc4_localizations_1 = pd.read_csv()
mat1 = sio.loadmat()
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_6 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_6.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_6.sort()
list_distances_6 = [i * 160 for i in list_distances_6]
#?  

#?  
Muc4_localizations_1 = pd.read_csv()
mat1 = sio.loadmat()
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_7 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_7.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_7.sort()
list_distances_7 = [i * 160 for i in list_distances_7]
#?  

#?  
Muc4_localizations_1 = pd.read_csv()
mat1 = sio.loadmat()
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_8 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_8.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_8.sort()
list_distances_8 = [i * 160 for i in list_distances_8]
#? 

#?  
Muc4_localizations_1 = pd.read_csv()
mat1 = sio.loadmat()
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_9 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_9.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_9.sort()
list_distances_9 = [i * 160 for i in list_distances_9]
#? 



lists_distances = [np.array(list_distances_1), np.array(list_distances_2), np.array(list_distances_3), np.array(list_distances_4), np.array(list_distances_5), np.array(list_distances_6), np.array(list_distances_7), np.array(list_distances_8), np.array(list_distances_9)]
ll = pd.DataFrame(lists_distances).mean()
ss = pd.DataFrame(lists_distances).sem()
ll = ll.groupby(ll.index // BINS).mean()
ss = ss.groupby(ss.index // BINS).mean()

ll = ll.loc[:NUMBER_OF_CLUSTERS]
ss = ss.loc[:NUMBER_OF_CLUSTERS]

fig, ax1 = plt.subplots(figsize=(7,5))
r1 = np.arange(len(ll))
r2 = [x + 0.3 for x in r1]
plt.bar(r1, ll, width = 0.3, color = 'red', edgecolor = 'black', yerr=ss, capsize=3, label='-IL6', alpha=0.8)



#!!!!  +IL6
#?                         
Muc4_localizations_1 = pd.read_csv()
mat1 = sio.loadmat()
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_1 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_1.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_1.sort()
list_distances_1 = [i * 160 for i in list_distances_1]
#?

#?
Muc4_localizations_1 = pd.read_csv()
mat1 = sio.loadmat()
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_2 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_2.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_2.sort()
list_distances_2 = [i * 160 for i in list_distances_2]
#?

#?
Muc4_localizations_1 = pd.read_csv()
mat1 = sio.loadmat()
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_3 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_3.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_3.sort()
list_distances_3 = [i * 160 for i in list_distances_3]
#?

#?
Muc4_localizations_1 = pd.read_csv()
mat1 = sio.loadmat()
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_4 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_4.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_4.sort()
list_distances_4 = [i * 160 for i in list_distances_4]
#?

#?
Muc4_localizations_1 = pd.read_csv()
mat1 = sio.loadmat()
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_5 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_5.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_5.sort()
list_distances_5 = [i * 160 for i in list_distances_5]
#?

#?
Muc4_localizations_1 = pd.read_csv()
mat1 = sio.loadmat()
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_6 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_6.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_6.sort()
list_distances_6 = [i * 160 for i in list_distances_6]
#?

#?
Muc4_localizations_1 = pd.read_csv()
mat1 = sio.loadmat()
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_7 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_7.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_7.sort()
list_distances_7 = [i * 160 for i in list_distances_7]
#?

#?
Muc4_localizations_1 = pd.read_csv()
mat1 = sio.loadmat()
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_8 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_8.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_8.sort()
list_distances_8 = [i * 160 for i in list_distances_8]
#?

#?
Muc4_localizations_1 = pd.read_csv()
mat1 = sio.loadmat()
centroid_1 = (Muc4_localizations_1.iloc[:, 0].sum()/len(Muc4_localizations_1.iloc[:, 0]) , Muc4_localizations_1.iloc[:, 1].sum()/len(Muc4_localizations_1.iloc[:, 1]))
mat_xydata1 = mat1['ClusterResults']
list_distances_9 = list()
df_1 = pd.DataFrame(mat_xydata1).iloc[:, 0:2]
for i in range(len(df_1)):
    xy_1 = tuple(df_1.iloc[i, 0:2])
    list_distances_9.append(math.sqrt((centroid_1[0] - xy_1[0])**2 + (centroid_1[1] - xy_1[1])**2))
list_distances_9.sort()
list_distances_9 = [i * 160 for i in list_distances_9]
#?
   

#*                                                                        
#*                Files reading                       
#*   


# Binning the results and plotting

lists_distances_2 = [np.array(list_distances_1), np.array(list_distances_2), np.array(list_distances_3), np.array(list_distances_4), np.array(list_distances_5), np.array(list_distances_6), np.array(list_distances_7), np.array(list_distances_8), np.array(list_distances_9)]
ll_2 = pd.DataFrame(lists_distances_2).mean()
ss_2 = pd.DataFrame(lists_distances_2).sem()
ll_2 = ll_2.groupby(ll_2.index // BINS).mean()
ss_2 = ss_2.groupby(ss_2.index // BINS).mean()

ll_2 = ll_2.loc[:NUMBER_OF_CLUSTERS]
ss_2 = ss_2.loc[:NUMBER_OF_CLUSTERS]

plt.bar(r2, ll_2, width = 0.3, color = 'blue', edgecolor = 'black', yerr=ss_2, capsize=3, label='+IL6', alpha=0.8)

plt.legend()
plt.title('Neighbors distances. Closest clusters')
plt.xticks(rotation=0)

plt.show()

