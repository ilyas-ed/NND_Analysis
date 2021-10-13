

'''

Written by Ilyas Ed-daoui, ilyas@gibh.ac.cn



The script loads the centroids obtained from the cluster analysis that has been published in this paper:
‘
Ricci MA, Manzo C, García-Parajo MF, Lakadamyali M, Cosma MP. 
Chromatin fibers are formed by heterogeneous groups of nucleosomes in vivo. Cell. 
2015 Mar 12;160(6):1145-58. doi: 10.1016/j.cell.2015.01.054. PMID: 25768910.
‘

Then for both molecules, the nanodomains centers are paired and the distance is calculated. The clusters pairs 
are binned by the cluster size and the frequency of each combination of cluster pairs was represented in a heat map.
The values are normalized to sum to unity.

'''


from tqdm import tqdm

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm

import scipy.io as sio

import seaborn as sns

from sklearn.neighbors import BallTree

from statistics import median



#! **********     Global Parameters      #! **********

THRESHOLD = 200 # The threshold on the NND.  
VALUE = np.nan # Value replacing all the other values above the threshold 
NBR_LOCS_BINS = list(range(5, 1500, 50))  # [0, 25, 50, 75, 100] # The bins for the first_molecule / second_molecule 's number of localization's heatmap 

#! **********     Global Parameters      #! **********


#*                                                                        
#*                first_molecule Files reading                       
#*   


paths_first_molecule_ = [] # the paths list to the mat files for the first molecule
paths_second_molecule_ = [] # the paths list to the mat files for the second molecule
   
merged_first_molecule_second_molecule_AllCells = pd.DataFrame()

for path_first_molecule_, path_second_molecule_ in tqdm(zip(paths_first_molecule_, paths_second_molecule_)):

#! ************************** first_molecule **************************

    matfirst_molecule = sio.loadmat(path_first_molecule_)
    matfirst_molecule_ClusterResults = matfirst_molecule['ClusterResults']
    matfirst_molecule_xynData = matfirst_molecule['xynData']
    clusterAreafirst_molecule_ = matfirst_molecule_xynData['clusterArea']
    clusterAreafirst_molecule = [float(elem) for elem in clusterAreafirst_molecule_[0][0]]

    X_first_molecule_Clust_centr=list()
    Y_first_molecule_Clust_centr=list()
    clusterAreafirst_molecule_list=list()

    for elem, nl, ca in zip(matfirst_molecule_ClusterResults, clusterAreafirst_molecule):
        X_first_molecule_Clust_centr.append(elem[0])
        Y_first_molecule_Clust_centr.append(elem[1])
        clusterAreafirst_molecule_list.append(ca)


    first_molecule_cluster_ = pd.DataFrame({
            'id_first_molecule': i,
            'x_first_molecule': a,
            'y_first_molecule': b,
            'clusterArea_first_molecule': d
        } for a, b, c, d, i in zip(X_first_molecule_Clust_centr, Y_first_molecule_Clust_centr, clusterAreafirst_molecule_list, range(len(clusterAreafirst_molecule_list)))
    )

#! ************************** first_molecule **************************


#? ************************** second_molecule **************************

    matsecond_molecule = sio.loadmat(path_second_molecule_)
    matsecond_molecule_ClusterResults = matsecond_molecule['ClusterResults']
    matsecond_molecule_xynData = matsecond_molecule['xynData']
    clusterAreasecond_molecule_ = matsecond_molecule_xynData['clusterArea']
    clusterAreasecond_molecule = [float(elem) for elem in clusterAreasecond_molecule_[0][0]]


    X_second_molecule_Clust_centr=list()
    Y_second_molecule_Clust_centr=list()
    clusterAreasecond_molecule_list=list()


    for elem, nl, ca in zip(matsecond_molecule_ClusterResults, clusterAreasecond_molecule):
        X_second_molecule_Clust_centr.append(elem[0])
        Y_second_molecule_Clust_centr.append(elem[1])
        clusterAreasecond_molecule_list.append(ca)


    second_molecule_cluster_ = pd.DataFrame({
            'id_second_molecule': k,
            'x_second_molecule': e,
            'y_second_molecule': f,
            'clusterArea_second_molecule': h
        } for e, f, g, h, k in zip(X_second_molecule_Clust_centr, Y_second_molecule_Clust_centr, clusterAreasecond_molecule_list, range(len(clusterAreasecond_molecule_list)))
    )

#? ************************** second_molecule **************************


#* ************************** first_molecule and second_moleculemerase NND calc **************************

    tree = BallTree(second_molecule_cluster_[['x_second_molecule', 'y_second_molecule']].values, leaf_size=2)
    first_molecule_cluster_['distance'], first_molecule_cluster_['id_second_molecule'] = tree.query(first_molecule_cluster_[['x_first_molecule', 'y_first_molecule']].values, k=1)
    merged_first_molecule_second_molecule = second_molecule_cluster_.merge(first_molecule_cluster_, on='id_second_molecule')

    merged_first_molecule_second_molecule_AllCells = merged_first_molecule_second_molecule_AllCells.append(merged_first_molecule_second_molecule)

#* ************************** first_molecule and second_moleculemerase NND calc **************************


print('Finished reading first_molecule and second_moleculemerase files ...')

merged_first_molecule_second_molecule_AllCells['distance'] = merged_first_molecule_second_molecule_AllCells['distance']*160
print(merged_first_molecule_second_molecule_AllCells['distance'])


#*                                                                        
#*                first_molecule Files reading                       
#*   



#!                                                                        
#!                first_molecule and second_moleculemerase cross analysis                       
#! 


clusterArea_map = merged_first_molecule_second_molecule_AllCells.groupby(['clusterArea_second_molecule', 'clusterArea_first_molecule']).size().reset_index(name="counts")
# clusterArea_map_map = clusterArea_map.pivot('clusterArea_second_molecule', 'clusterArea_first_molecule', 'counts') #! uncomment this line if results are not binned
# numberLocs_map = numberLocs_map.clip(0, 200)


#?******************************************************************************************
#*            Binning the results of the number of localizations

clusterArea_map_binned = clusterArea_map
clusterArea_map_binned['clusterArea_second_molecule_binned'] = pd.cut(clusterArea_map['clusterArea_second_molecule'], NBR_LOCS_BINS)
clusterArea_map_binned['clusterArea_first_molecule_binned'] = pd.cut(clusterArea_map['clusterArea_first_molecule'], NBR_LOCS_BINS)
clusterArea_map_binned_toplot = clusterArea_map_binned.groupby(['clusterArea_second_molecule_binned', 'clusterArea_first_molecule_binned'])['counts'].agg('sum').reset_index(name="counts_bin_sum")
clusterArea_map_binned_toplot = clusterArea_map_binned_toplot.pivot('clusterArea_second_molecule_binned', 'clusterArea_first_molecule_binned', 'counts_bin_sum')

#*            Binning the results of the number of localizations
#?******************************************************************************************


merged_first_molecule_second_molecule_AllCells.distance[merged_first_molecule_second_molecule_AllCells['distance'] > THRESHOLD] = VALUE
distances_map = merged_first_molecule_second_molecule_AllCells.groupby(['clusterArea_second_molecule', 'clusterArea_first_molecule'])['distance'].agg('mean').reset_index(name="distance_mean")
# distances_map['distance_mean'] = distances_map['distance_mean'] * 160
# distances_map_map = distances_map.pivot('clusterArea_second_molecule', 'clusterArea_first_molecule', 'distance_mean') #! uncomment this line if results are not binned


#?******************************************************************************************
#*            Binning the results of distances

distances_map_binned = distances_map
distances_map_binned['clusterArea_second_molecule_binned'] = pd.cut(distances_map['clusterArea_second_molecule'], NBR_LOCS_BINS)
distances_map_binned['clusterArea_first_molecule_binned'] = pd.cut(distances_map['clusterArea_first_molecule'], NBR_LOCS_BINS)
distances_map_binned_toplot = distances_map_binned.groupby(['clusterArea_second_molecule_binned', 'clusterArea_first_molecule_binned'])['distance_mean'].agg('mean').reset_index(name="distance_bin_mean")
distances_map_binned_toplot = distances_map_binned_toplot.pivot('clusterArea_second_molecule_binned', 'clusterArea_first_molecule_binned', 'distance_bin_mean')

#*            Binning the results of distances
#?******************************************************************************************


#?               Max normalization                
# print('Maximum number of localizations in a cluster is : ', numberLocs_map.max().max())
# numberLocs_map = numberLocs_map/numberLocs_map.max().max()
#?               Max normalization    


#?                normalization                
print('Sum of localizations in all clusters is : ', clusterArea_map.sum().sum())
clusterArea_map_binned_toplot = clusterArea_map_binned_toplot/clusterArea_map_binned_toplot.sum().sum()
#?                normalization    


#!                                                                        
#!                first_molecule and second_moleculemerase cross analysis                        
#! 


                                                       
#*                Plotting of the heatmap                       
#* 

ax1 = plt.subplot(121)
sns.heatmap(clusterArea_map_binned_toplot,cmap = 'jet',
            square=True, linewidths=.5)#, vmin=0, vmax=20

plt.title('Frequency of counts')
ax1.set_xlim(0, 51.5)
ax1.set_ylim(0, 105)

#*                                                                        
#*                Plotting of the heatmap                       


#?               KDE of distance    

plt.subplots()
sns.distplot(a=merged_first_molecule_second_molecule_AllCells['distance'], hist=True, kde=True, rug=False, 
                kde_kws={"color": "b", "alpha":0.25, "linewidth": 3, "shade":True }, bins=200, color='blue' 
            )
plt.title('NND distribution between second_molecule and first_molecule')
plt.xlim(xmin=0)

#?               KDE of distance   


plt.show()

