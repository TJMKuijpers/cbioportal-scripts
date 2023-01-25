#########################################################################
# Python file to add a stable entry id to the mutational matrix count   #
# as well as the metadata column that specifies the mutation class      #
# T.J.M. Kuijpers   2023-01-25                                           #
#########################################################################

import os
import pandas as pd
from re import search,split


# Read the mutational matrix count file
data_file = pd.read_csv('data_mutational_signature_matrix_DBS.txt',sep='\t')
dbs_file = True
sbs_file = False
id_file = False


# Map for the ID legend
insertions_deletions_map={
    '1:Del:C':'1bp deletion (C)',
    '1:Del:T':'1bp deletion (T)',
    '1:Ins:C':'1bp insertion (C)',
    '1:Ins:T':'1bp insertion (T)',
    '2:Del:R':'2bp deletion at repeats',
    '3:Del:R':'3bp deletion at repeats',
    '4:Del:R':'4bp deletion at repeats',
    '5:Del:R':'5bp deletion at repeats ',
    '2:Ins:R':'2bp insertion at repeats',
    '3:Ins:R':'3bp insertion at repeats',
    '4:Ins:R':'4bp insertion at repeats',
    '5:Ins:R':'5bp insertion at repeats',
    '2:Del:M':'Microhomology (Del. length 2)',
    '3:Del:M':'Microhomology (Del. length 3)',
    '4:Del:M':'Microhomology (Del. length 4)',
    '5:Del:M':'Microhomology (Del. length 5)'
}

# For every file check iif the stable id column exists
if 'ENTITY_STABLE_ID' not in data_file.columns:
    stable_id = ['mutational_signature_matrix_'+data_file['MUTATION_TYPE'][x] for x in range(data_file['MUTATION_TYPE'].__len__())]
    data_file.insert(0,'ENTITY_STABLE_ID',stable_id)


if dbs_file:
    patterns_dbs = ['AC>', 'AT>', 'CC>', 'CG>', 'CT>', 'GC>', 'TA>', 'TC>', 'TG>', 'TT>']
    mutation_class = []
    # add the mutation class properties
    for x in patterns_dbs:
        sub_class = [str(x+'NN') for name in data_file.MUTATION_TYPE if search(x, name)]
        mutation_class.extend(sub_class)
    data_file.insert(1, 'MUTATION_CLASS', mutation_class)
    data_file.to_csv('data_mutational_signature_matrix_count_DBS.txt', sep='\t', index=False)

if sbs_file:
    # in the SBS file we can get the mutation class by performing a string split on the MUTATION_TYPE column
    mutation_class = []
    for x in data_file.MUTATION_TYPE:
        mutation_class.append(split('\[(.*?)\]',x)[1])
    data_file.insert(1, 'MUTATION_CLASS', mutation_class)
    data_file.to_csv('data_mutational_signature_matrix_count_SBS.txt', sep='\t', index=False)

if id_file:
    patterns_id = insertions_deletions_map.keys()
    mutation_class=[]
    for x in patterns_id :
        subclass= [insertions_deletions_map[x] for name in data_file.MUTATION_TYPE if search(x, name)]
        mutation_class.extend(subclass)
    data_file.insert(1, 'MUTATION_CLASS',  mutation_class)
    data_file.to_csv('data_mutational_signature_matrix_count_ID.txt', sep='\t', index=False)


