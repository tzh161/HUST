#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process data io
External APIs: 
    r_atom_loc_struct()
    r_feat_csv()
    gen_feats_rndImg()

Created on Wed Aug 22 10:58:18 PDT 2018
@author: lingmiao@lbl.gov at Prof. LinWang Wang's group
"""

import numpy as np
import pandas as pd

import os
import time
import sys

import parameters as pm
from func_io import r_atom_loc_struct, r_feat_csv
from func_struct import get_neighbor_struct_cp

if pm.ntypes == 1:
    from func_feat import get_feat_nbxyz, get_feat_dfeat_nbxyz
    print('pm.ntypes =1, from func_feat import get_feat_nbxyz, get_feat_dfeat_nbxyz')
else:
    from func_feat_ntype import get_feat_nbxyz, get_feat_dfeat_nbxyz
    print('pm.ntypes >1, from func_feat_ntype import get_feat_nbxyz, get_feat_dfeat_nbxyz')
#================================================================================

def get_neighbor_df_PFVEi(df_PFVEi, natoms, cell):
    itype = df_PFVEi[0:natoms,0]
    pos   = np.matmul(df_PFVEi[0:natoms, 1:4], cell.T)
    fors  = df_PFVEi[natoms*1:natoms*2, 1:4]
    engy  = df_PFVEi[natoms*3:natoms*4, 1]

    nNb, idxNb, nbxyz = get_neighbor_struct_cp(cell, pos)
    '''
    from func_struct_ase import get_neighbor_struct
    from ase.io import read
    atoms_str = read(pm.f_atoms)
    at_Pos   = df_PFVEi[0:natoms, 1:4]
    atoms_str.set_scaled_positions(at_Pos)
    itype, nNb, idxNb, nbxyz = get_neighbor_struct(atoms_str)
    '''
    nbxyz = nbxyz.reshape([-1, pm.maxNb*3])
    df_NL = np.concatenate([np.expand_dims(itype,1), \
            np.expand_dims(nNb,1), \
            np.expand_dims(engy,1), \
            fors, idxNb, nbxyz], axis=1)
    return df_NL

#================================================================================
def w_feats_file_PFVEi(f_atoms, f_PFVEi, idxPretrain, idxTrainVald, idxTest, f_feat):
    """read MOVEMENT_PFVEi, generate atomic local structure, and feature file
    """
    from ase.io import read
    atoms_str = read(f_atoms)
    natoms = len(atoms_str)
    cell = atoms_str.get_cell(complete=True)
    print('cell', cell)

    df_PFVEi = pd.read_csv(f_PFVEi, header=None,index_col=False,dtype="float", \
                    delim_whitespace=True)
    df_all = df_PFVEi.values.reshape([-1, natoms*4, 4])
    print("generate features R0, Rc, maxNb", pm.R0, pm.Rc, pm.maxNb, " at", time.ctime())
    print("checking, df_all ", df_all.shape, '\n', df_all[0:2,0:5,:], time.ctime())
    sys.stdout.flush()

    for i in idxTest:
        df_NL = get_neighbor_df_PFVEi(df_all[i], natoms, cell)
        df = pd.DataFrame(df_NL)
        w_feat_df2csv(df, f_feat+'_test')

    for i in idxTrainVald:
        df_NL = get_neighbor_df_PFVEi(df_all[i], natoms, cell)
        df = pd.DataFrame(df_NL)
        w_feat_df2csv(df, f_feat+'_train')
    
    nI_pre = idxPretrain.shape[0]
    os.system('head -' +str(nI_pre*natoms)+ ' ' +f_feat+'_train > ' +f_feat+'_pretrain')
    os.system('head -' +str(nI_pre*natoms)+ ' ' +f_feat+'_train.csv > ' +f_feat+'_pretrain.csv')
    
    print("Job finised at ", time.ctime())
    return 

#================================================================================
def MOVEMENT2f_PFVEi(f_MOVEMENT, f_PFVEi):
    """read MOVEMENT file, write to MOVEMENT_PFVEi
    NOTE: should be very caryfully used and double check f_PFVEi file!!!
       linux cmd like as:
       grep "^ *[1-9]" MOVEMENT | sed "/Iteration/d" | sed "s/    1  1  1//" | sed "s/    0  0  0//" | less
    """
    #cmd = 'grep "^ *[1-9]" ' +f_MOVEMENT+ ' | sed "/Iteration/d" | sed "s/    1  1  1//" | sed "s/    0  0  0//" > ' +f_PFVEi
    cmd = 'grep "^  14" ' +f_MOVEMENT+ ' | sed "/Iteration/d" | sed "s/     1  1  1//" | sed "s/     0  0  0//" > ' +f_PFVEi
    os.system(cmd)
    cmd = 'grep Iter ' +f_MOVEMENT+ ' > potentialEnergy'
    os.system(cmd)
    return
#================================================================================
# test
