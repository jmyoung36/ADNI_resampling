# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 16:04:04 2016

@author: jonyoung
"""

# import the stuff we need
import pandas as pd
import numpy as np
import os

# set directory to read and write data
dir = os.path.dirname(__file__)
metadata_dir = os.path.join(dir, 'Data/')

# looking at which set of subjects? AD/normal or MCI-s/MCI-c
subjects = 'AD/Normal'

# import ADNI timecourses
subject_diagnoses = pd.read_csv(metadata_dir + 'subject_diagnoses_PET.csv')

# look at only the chosen subset of subjects
if subjects == 'AD/normal' :
    
    subject_diagnoses = subject_diagnoses[np.logical_or(subject_diagnoses['diagnosis'] == 1, subject_diagnoses['diagnosis'] == 6)]
    
elif subjects == 'MCI' :
    
    subject_diagnoses = subject_diagnoses[np.logical_or(subject_diagnoses['diagnosis'] == 1, subject_diagnoses['diagnosis'] == 6)]
    
# split by scan field strength
subject_diagnoses_3T = subject_diagnoses[subject_diagnoses['T'] == 3]
subject_diagnoses_1pt5T = subject_diagnoses[subject_diagnoses['T'] == 1.5]
    
# how many scans are at 3T?
print 'Number of scans at 3T = ' + str(len(subject_diagnoses_3T))

# across how many sites?
print 'Number of sites with scans at 3T = ' + str(len(subject_diagnoses_3T['SITEID'].unique()))

grouped_3T = subject_diagnoses_3T.groupby('SITEID')


# how many scans are at 1.5T?
print 'Number of scans at 1.5T = ' + str(len(subject_diagnoses_1pt5T))

# across how many sites?
print 'Number of sites with scans at 1.5T = ' + str(len(subject_diagnoses_1pt5T['SITEID'].unique()))

grouped_1pt5T = subject_diagnoses_1pt5T.groupby('SITEID')
print grouped_1pt5T.size()

print 'Number of sites with 1.5T scans only = ' + str(sum(subject_diagnoses.groupby('SITEID').apply(lambda x: all(x['T'] == 1.5)) == True))
print 'Number of sites with 3T scans only = ' + str(sum(subject_diagnoses.groupby('SITEID').apply(lambda x: all(x['T'] == 3)) == True))
print 'Number of sites with both 1.5T and 3T scans = ' + str(sum(subject_diagnoses.groupby('SITEID').apply(lambda x: len(x['T'].unique()) == 2)))

