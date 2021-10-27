"""
save the intra-hemispheric (LL-RR) and inter-hemispheric (LR-RL) gradient
score differences for each subject in a solar-like *csv format
"""
import pandas as pd
import numpy as np
import os

codedir = os.path.dirname(os.path.abspath(__file__))

path = os.path.join(os.path.dirname(codedir), 
                       'data/data_results/gradient/')

path_list = os.listdir(os.path.join(path, 'LL-RR'))
path_list.sort()

solout = os.path.join(os.path.dirname(codedir), 
                       'data/data_results/solar')

# read-out intra-hemispheric (LL-RR) gradient score diff for each subject
# and save as *csv along primary, secondary and third gradients
G1 = np.array(pd.read_csv(os.path.join(path, 'LL-RR/100206.csv'), header=None))[:,0]
G2 = np.array(pd.read_csv(os.path.join(path, 'LL-RR/100206.csv'), header=None))[:,1]
G3 = np.array(pd.read_csv(os.path.join(path, 'LL-RR/100206.csv'), header=None))[:,2]

for i in path_list:
    G = np.array(pd.read_csv(os.path.join(path, 'LL-RR', str(i)), header=None))
    G_1 = G[:,0]
    G_2 = G[:,1]
    G_3 = G[:,2]
    G1 = np.vstack((G1, G_1))
    G2 = np.vstack((G2, G_2))
    G3 = np.vstack((G3, G_3))

# COMMENT
G_1 = G1[1:]
G_2 = G2[1:]
G_3 = G3[1:]

header = [None]*180
for i in range (180):
    header[i] = 'node_'+str(i+1)

final_1 = np.vstack((np.array(header), G_1)).T
final_2 = np.vstack((np.array(header), G_2)).T
final_3 = np.vstack((np.array(header), G_3)).T

np.savetxt(os.path.join(solout, 'LL-RR/G1/G1.csv'), final_1.T, 
           delimiter = ',',fmt = '%s')
np.savetxt(os.path.join(solout, 'LL-RR/G2/G2.csv'), final_2.T, 
           delimiter = ',',fmt = '%s')
np.savetxt(os.path.join(solout, 'LL-RR/G3/G3.csv'), final_3.T, 
           delimiter = ',',fmt = '%s')

# read-out inter-hemispheric (LR-RL) gradient score diff for each subject
# and save as *csv along primary, secondary and third gradients
G1 = np.array(pd.read_csv(os.path.join(path, 'LR-RL/100206.csv'), header=None))[:,0]
G2 = np.array(pd.read_csv(os.path.join(path, 'LR-RL/100206.csv'), header=None))[:,1]
G3 = np.array(pd.read_csv(os.path.join(path, 'LR-RL/100206.csv'), header=None))[:,2]

for i in path_list:
    G = np.array(pd.read_csv(os.path.join(path, 'LR-RL', str(i)), header=None))
    G_1 = G[:,0]
    G_2 = G[:,1]
    G_3 = G[:,2]
    G1 = np.vstack((G1, G_1))
    G2 = np.vstack((G2, G_2))
    G3 = np.vstack((G3, G_3))

G_1 = G1[1:]
G_2 = G2[1:]
G_3 = G3[1:]

final_1 = np.vstack((np.array(header), G_1)).T
final_2 = np.vstack((np.array(header), G_2)).T
final_3 = np.vstack((np.array(header), G_3)).T
np.savetxt(os.path.join(solout, 'LR-RL/G1/G1.csv'), final_1.T, 
           delimiter = ',',fmt = '%s')
np.savetxt(os.path.join(solout, 'LR-RL/G2/G2.csv'), final_2.T, 
           delimiter = ',',fmt = '%s')
np.savetxt(os.path.join(solout, 'LR-RL/G3/G3.csv'), final_3.T, 
           delimiter = ',',fmt = '%s')
