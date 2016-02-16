# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 10:46:42 2016

@author: jonyoung
"""

# import the stuff we need
import pandas as pd
from scipy.io import loadmat
import numpy as np
import os


# set directory to read and write data
dir = os.path.dirname(__file__)
metadata_dir = os.path.join(dir, 'Data/')


# map to convert from DXCHANGE to DXCURREN
DXCHANGE_to_DXCURREN = {1:1, 2:2, 3:3, 4:2, 5:3, 6:3, 7:1, 8:2, 9:1}

# read diagnoses into a pandas dataframe
diagnoses = pd.read_csv(metadata_dir + 'DXSUM_PDXCONV_ADNIALL_new.csv')

# read scan field strengths into pandas dataframe
scanty = pd.read_csv(metadata_dir + 'scanty.csv', header=None)
scanty.columns = ['RID','T']

# MCI conversion cutoff point
cutoff = 'm36'


# read in matlab .mat file kernel and list of RIDs
K = loadmat(metadata_dir + 'K_FDG.mat')
MRI_diags_goodreg_all = loadmat(metadata_dir + 'MRI_diags_goodreg_all.mat')
valid_registered_RID = MRI_diags_goodreg_all['diags'][:,0]

print K.keys()

kernel_RID = K['K_FDG'][0][0][0]
kernel_vals = K['K_FDG'][0][0][1]
print np.shape(kernel_RID)
print np.shape(kernel_vals)

# retain only the columns we are interested in
diagnoses = diagnoses[['Phase', 'ID', 'RID', 'SITEID', 'VISCODE', 'VISCODE2', 'EXAMDATE', 'DXCHANGE', 'DXCURREN']]

# get rid of sc & uns1 entries
diagnoses = diagnoses[diagnoses['VISCODE2'] != 'sc']
diagnoses = diagnoses[diagnoses['VISCODE2'] != 'uns1']

# group rows by RID to gather together all entries for a particular individual
diagnoses_grouped = pd.groupby(diagnoses, 'RID')

rows_list = []

# sort all entries within each group by date
for RID, group in diagnoses_grouped :
    
    row = {}
    diagnoses_dict = {'bl':0, 'm06':0, 'm12':0, 'm18':0, 'm24':0, 'm36':0, 'm48':0, 'm60':0, 'm72':0, 'm84':0, 'm96':0, 'm108':0, 'bl date':0, 'm06 date':0, 'm12 date':0, 'm18 date':0, 'm24 date':0, 'm36 date':0, 'm48 date':0, 'm60 date':0, 'm72 date':0, 'm84 date':0, 'm96 date':0, 'm108 date':0, 'last':0}
    
    group.sort_values(['EXAMDATE'])
    
    # Phase, ID, RID and SITEID all come from the first diagnosis
    first_diagnosis = group.iloc[0]
    row.update({'phase':first_diagnosis['Phase'], 'RID':first_diagnosis['RID'], 'ID':first_diagnosis['ID'], 'SITEID':first_diagnosis['SITEID']})
    rows_list.append(row)
            
    # now iterate through diagnoses
    for i in range(len(group)) :
        
        group_row = group.iloc[i]
        
        # get the VISCODE2
        VISCODE2 = group_row['VISCODE2']
        
        # get the date
        EXAMDATE = group_row['EXAMDATE']
        
        # get the diagnoses
        DXCURREN = group_row['DXCURREN']
        
        if pd.isnull(DXCURREN) :
            
            DXCURREN = DXCHANGE_to_DXCURREN[group_row['DXCHANGE']]
            
        diagnoses_dict[VISCODE2] = DXCURREN
        diagnoses_dict[VISCODE2 + ' date'] = EXAMDATE
        diagnoses_dict['last'] = VISCODE2
        
    row.update(diagnoses_dict)
scanty.ix[scanty['T'] == 2.9, 'T'] = 3        
subject_diagnoses = pd.DataFrame(rows_list)
subject_diagnoses = pd.merge(subject_diagnoses, scanty, on='RID')

#print(len(subject_diagnoses))
# filter out RIDs not in list of properly registered images
subject_diagnoses = subject_diagnoses[subject_diagnoses['RID'].isin(valid_registered_RID)]
#print(len(subject_diagnoses))

# filter out RIDs not in list of images in the kernel
subject_diagnoses = subject_diagnoses[subject_diagnoses['RID'].isin(np.squeeze(kernel_RID))]
#print(len(subject_diagnoses))

# reorder columns
subject_diagnoses = subject_diagnoses[['ID', 'RID', 'SITEID', 
'phase', 'T', 'bl', 'm06', 'm12', 'm18', 'm24', 'm36', 'm48', 'm60', 'm72', 'm84', 'm96', 'm108', 'bl date', 'm06 date', 'm12 date', 'm18 date', 'm24 date', 'm36 date', 'm48 date', 'm60 date', 'm72 date', 'm84 date', 'm96 date', 'm108 date', 'last']]

# label each subject as healthy, AD or various flavours of MCI
subject_diagnoses['diagnosis'] = 0
indices = subject_diagnoses.index

for i in range(len(subject_diagnoses)) :
    
    diagnosis_row =  subject_diagnoses.iloc[i]
    ind = indices[i]
    #print i
    #print ind
    #print diagnosis_row
    
    # AD and normal are easy
    # 5 different types of MCI so call normal & AD 1 and 6
    # dropout = 0
    # reverter = 2
    # stable = 3
    # converter = 4
    # late converter = 5
    if diagnosis_row['bl'] == 1 :
        
        subject_diagnoses.ix[ind, 'diagnosis'] = 1
        
    elif diagnosis_row['bl'] == 3 :
        
        subject_diagnoses.ix[ind, 'diagnosis'] = 6
        
    # MCI  is a bit more complicted...
    else :
        
        # dropouts have a final timepoint less than the cutoff, which is MCI

        if diagnosis_row['last'] == 'bl' or( int(diagnosis_row['last'][1:]) < int(cutoff[1:]) and diagnosis_row[diagnosis_row['last']] == 2) :
            
            subject_diagnoses.ix[ind, 'diagnosis'] = 0
            
        # if any timepoint later than bl is normal, it is a reverter
        elif any(timepoint == 1 for timepoint in diagnosis_row[['m06', 'm12', 'm18', 'm24', 'm36', 'm48', 'm60', 'm72', 'm84', 'm96', 'm108']].tolist())  :
            #print diagnosis_row['RID']
            #print diagnosis_row[['m06', 'm12', 'm18', 'm24', 'm36', 'm48', 'm60', 'm72', 'm84', 'm96', 'm108']].tolist()
            subject_diagnoses.ix[ind, 'diagnosis'] = 2
            
        # if not either, loop through all the timepoints
        else :
            
            timepoints = ['m06', 'm12', 'm18', 'm24', 'm36', 'm48', 'm60', 'm72', 'm84', 'm96', 'm108']            
            timepoint_diagnosis_list = diagnosis_row[timepoints].tolist()
            
            # assume stable by default
            current_MCI_group = 3 
            
            for timepoint in timepoints :
                
                # current group is stable MCI and next timepoint is MCI.
                # stay stable MCI and end at cutoff
                if current_MCI_group == 3 and diagnosis_row[timepoint] == 2 :
                    if diagnosis_row['RID'] == 54 :
                        print timepoint
                        print current_MCI_group
                    current_MCI_group = 3
                    if timepoint == cutoff :
                        break
                    
                # current group is stable MCI and next timepoint is AD and we are before cutoff
                # become converter but carry on to look for reversion
                elif current_MCI_group == 3 and diagnosis_row[timepoint] == 3 and timepoint <= cutoff :
                    if diagnosis_row['RID'] == 54 :
                        print timepoint
                        print current_MCI_group
                    current_MCI_group = 4
                    
                # current group is stable MCI and next timepoint is AD and we are after cutoff
                # become late converter but carry on to look for reversion
                elif current_MCI_group == 3 and diagnosis_row[timepoint] == 3 and timepoint > cutoff :
                    if diagnosis_row['RID'] == 54 :
                        print timepoint
                        print current_MCI_group
                    current_MCI_group = 5
                    
                # current group us is stable MCI and next timepoint is 0 (no diagnosis)
                # stay stable MCI but continue to look for late conversion
                elif current_MCI_group == 3 and diagnosis_row[timepoint] == 0 :
                    if diagnosis_row['RID'] == 54 :
                        print timepoint
                        print current_MCI_group
                    current_MCI_group = 3
                    
                # current group is converter and next timepoint is MCI
                # become reverter and break
                elif current_MCI_group == 4 and diagnosis_row[timepoint] == 2 :
                    if diagnosis_row['RID'] == 54 :
                        print timepoint
                        print current_MCI_group
                    current_MCI_group = 2
                    break
                
                # current group is converter and next timepoint is AD
                # stay converter but continue to look for reversion
                elif current_MCI_group == 4 and diagnosis_row[timepoint] == 3 :
                    if diagnosis_row['RID'] == 54 :
                        print timepoint
                        print current_MCI_group
                    current_MCI_group = 4
                    
                # current group is converter and next timepoint is 0 (no diagnosis)
                # stay converter but continue to look for reversion
                elif current_MCI_group == 4 and diagnosis_row[timepoint] == 0 :
                    if diagnosis_row['RID'] == 54 :
                        print timepoint
                        print current_MCI_group
                    current_MCI_group = 4
                    
                # current group is late converter and and timepoint is 0 (no diagnosis)
                # stay late converter but continue to look for reversion
                elif current_MCI_group == 5 and diagnosis_row[timepoint] == 0 :
                    if diagnosis_row['RID'] == 54 :
                        print timepoint
                        print current_MCI_group
                    current_MCI_group = 5
                    
                # current group is late converter and next timepoint is MCI
                # become reverter and break
                elif current_MCI_group == 5 and diagnosis_row[timepoint] == 2 :
                    if diagnosis_row['RID'] == 54 :
                        print timepoint
                        print current_MCI_group
                    current_MCI_group = 2
                    break
                
                # current group is late converter and next timepoint is AD
                # continue to look for reversion
                elif current_MCI_group == 5 and diagnosis_row[timepoint] == 3 : 
                    if diagnosis_row['RID'] == 54 :
                        print timepoint
                        print current_MCI_group
                    current_MCI_group = 5
                    
            subject_diagnoses.ix[ind, 'diagnosis'] = current_MCI_group

subject_diagnoses.to_csv(metadata_dir + 'subject_diagnoses_PET.csv')
