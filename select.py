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


#================================================================================
def w_feats_file_PFVEi(f_PFVEi,f_feat):
    """read MOVEMENT_PFVEi, generate atomic local structure, and feature file
    """

    natoms = 18


    df_PFVEi = pd.read_csv(f_PFVEi, header=None,index_col=False,dtype="float", \
                    delim_whitespace=True)
    df_all = df_PFVEi.values.reshape([-1, natoms, 7])
    
    
    
    
    #ixyz_xyz = np.zeros((natoms,6))
    #ixyz_type= [['C'],['C'],['C'],['C'],['C'],['C'],['C'],['C'],['C'],['C'],['H'],['H'],['H'],['H'],['H'],['H'],['H'],['H']]
   # ixyz_type=np.array(ixyz_type)
    #f_test=open(f_feat+'_test.xyz','wb')
   # for i in range(500):
        #ixyz_xyz[:,0:6]=df_all[i,:,1:7]
        #ixyz = np.hstack((ixyz_type,ixyz_xyz))
        #np.savetxt(f_test,ixyz,fmt='%s %s %s %s %s %s %s',delimiter=' ',\
                       #header=str(natoms)+'\n',comments='')
   # f_test.close()
    #print("checking, df_all ", df_all.shape, '\n', df_all[0:2,0:18,:], time.ctime())   
    #sys.stdout.flush()



    n_atoms = df_all.shape[0]
    ixyz_type= [['C'],['C'],['C'],['C'],['C'],['C'],['H'],['H'],['H'],['H'],['H'],['H']]
    ixyz_type=np.array(ixyz_type)
    ixyz_xyz = np.zeros((12,6))
    
    
    
    
    f_1=open(f_feat+'_1.xyz','wb')
    for i in range(n_atoms):
        ixyz_xyz[0,:]=df_all[i,0,1:7]
        ixyz_xyz[1,:]=df_all[i,1,1:7]
        ixyz_xyz[2,:]=df_all[i,2,1:7]
        ixyz_xyz[3,:]=df_all[i,3,1:7]
        ixyz_xyz[4,:]=df_all[i,4,1:7]
        ixyz_xyz[5,:]=df_all[i,9,1:7]
        ixyz_xyz[6,:]=df_all[i,10,1:7]
        ixyz_xyz[7,:]=df_all[i,11,1:7]
        ixyz_xyz[8,:]=df_all[i,12,1:7]
        ixyz_xyz[9,:]=df_all[i,13,1:7]
        ixyz_xyz[10,:]=df_all[i,5,1:7]
        for j in range(3):
            ixyz_xyz[10,j] = ixyz_xyz[10,j]*0.77-ixyz_xyz[4,j]*0.77+ixyz_xyz[4,j]
        ixyz_xyz[11,:]=df_all[i,8,1:7]
        for j in range(3):
            ixyz_xyz[11,j] = ixyz_xyz[11,j]*0.77-ixyz_xyz[5,j]*0.77+ixyz_xyz[5,j]
        ixyz = np.hstack((ixyz_type,ixyz_xyz))
        np.savetxt(f_1,ixyz,fmt='%s %s %s %s %s %s %s',delimiter=' ',\
                       header=str(12)+'\n'+str(111111.111111),comments='')
    f_1.close()
    
    
    f_2=open(f_feat+'_2.xyz','wb')
    for i in range(n_atoms):
        ixyz_xyz[0,:]=df_all[i,8,1:7]
        ixyz_xyz[1,:]=df_all[i,9,1:7]
        ixyz_xyz[2,:]=df_all[i,4,1:7]
        ixyz_xyz[3,:]=df_all[i,5,1:7]
        ixyz_xyz[4,:]=df_all[i,6,1:7]
        ixyz_xyz[5,:]=df_all[i,7,1:7]
        ixyz_xyz[6,:]=df_all[i,17,1:7]
        ixyz_xyz[7,:]=df_all[i,0,1:7]
        for j in range(3):
            ixyz_xyz[7,j] = ixyz_xyz[7,j]*0.77-ixyz_xyz[1,j]*0.77+ixyz_xyz[1,j]
        ixyz_xyz[8,:]=df_all[i,3,1:7]
        for j in range(3):
            ixyz_xyz[8,j] = ixyz_xyz[8,j]*0.77-ixyz_xyz[2,j]*0.77+ixyz_xyz[2,j]
        ixyz_xyz[9,:]=df_all[i,14,1:7]
        ixyz_xyz[10,:]=df_all[i,15,1:7]
        ixyz_xyz[11,:]=df_all[i,16,1:7]
        ixyz = np.hstack((ixyz_type,ixyz_xyz))
        np.savetxt(f_2,ixyz,fmt='%s %s %s %s %s %s %s',delimiter=' ',\
                       header=str(12)+'\n'+str(111111.111111),comments='')
    f_2.close()
    
    
        
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
    cmd = 'grep "^[CHO]" ' +f_MOVEMENT+ ' | sed "/Iteration/d" | sed "s/C/6/" | sed "s/O/8/" | sed "s/H/1/" > ' +f_PFVEi
    os.system(cmd)
    #cmd = 'grep Iter ' +f_MOVEMENT+ ' > potentialEnergy'
    #os.system(cmd)
    return
#================================================================================
    
