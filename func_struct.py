#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 16:00:24 PDT 2018
@author: lingmiao@lbl.gov at Prof. LinWang Wang's group
"""


import numpy as cp
import numpy as np

#================================================================================
def get_neighbor_struct_cp(cell, pos, R0=0.2, Rc=8.0, maxNb=60):
    '''Analysis atomic structure to get local neighbor structure for each atom.
         If num_nb_atoms > maxNb in Rc, neighbors_list are sorted with dij, 
         and will remain nearest maxNb atoms, and may generate un-exact feature! 
       parameters:
         cell: 3x3 cell lattice, Ang
         pos: nx3 atomic postion, Ang
       return:
         nNb:  nx1, number of neighbor atoms in Rc
         idxNb: nxm, index (in cell) of neighbor atoms in Rc
         nbxyz: nxmx3, relative postion of neighbor atoms i
    '''
    natoms = pos.shape[0]

    # like in ase.atoms class
    # NOTE!!! need debug for cupy version
    #Rcry = cp.linalg.solve(cell.T, pos.T).T
    Rcry = cp.zeros_like(pos)
    Rcry[:,0] = pos[:,0]/cell[0,0]
    Rcry[:,1] = pos[:,1]/cell[1,1]
    Rcry[:,2] = pos[:,2]/cell[2,2]

    Rcry %= 1.0
    Rcry %= 1.0
    #_Rcry = Rcry.copy()
    #_Rcry[Rcry >1] = Rcry[Rcry >1] - cp.floor(Rcry[Rcry >1])
    #_Rcry[Rcry <0] = Rcry[Rcry <0] - cp.floor(Rcry[Rcry <0])
    #print('d_Rcry', cp.sum(cp.abs(Rcry -_Rcry)))

    Rd = Rcry[None,:,:] - Rcry[:,None,:]
    Rd[Rd> 0.5] = Rd[Rd> 0.5] - 1
    Rd[Rd<-0.5] = Rd[Rd<-0.5] + 1
    # to wraped in cell with Angstrom
    Rd = cp.matmul(Rd, cell.T)   
    
    Rd2 = cp.sum(Rd**2,axis=2)
    Rd[Rd2>Rc**2] = 0
    # idxMask also exclude i-th atom-self, comparing with RcMask
    idxMask = cp.sum(Rd**2,axis=2)>0
    nNb = idxMask.sum(axis=1)

    if cp.min(Rd2[Rd2>0]) < R0**2:
        iR0 = cp.where((0<Rd2) * (Rd2<R0**2))  # True*True = Ture
        print('===Warning in get_neighbor_struct_cp():', \
                iR0, 'th atoms have ', Rd2[iR0], 'distanse < R0', R0, \
                " Will generate wrong feature!" )
    if cp.max(nNb) > maxNb:
        # to sort as d and make a cutoff of maxNb atoms
        isort = cp.where(nNb>maxNb)
        print('===Warning in get_neighbor_struct_cp():', \
                cp.where(nNb>maxNb)[0], 'th atoms have ', nNb[nNb>maxNb], 'neighbors > maxNb', maxNb, \
                " Will remain nearest maxNb atoms, and generate un-exact feature!" )
        for i in cp.where(nNb>maxNb)[0]:
            di = cp.sum(Rd[i]**2,axis=1)
            arg_id    = cp.argsort(di)
            Rd[i, arg_id[natoms-nNb[i]+maxNb:]] = 0
            #print("Rd[i] !=0 ", cp.sum(Rd[i, arg_id[natoms-nNb[i]-1:]]**2,axis=1))
        idxMask = cp.sum(Rd**2,axis=2)>0
        nNb = idxMask.sum(axis=1)

    iidx, jidx = cp.where(idxMask)
    jidx2 = cp.concatenate([cp.arange(0, int(nNb[i])) for i in range(natoms)]).astype(int)

    idxNb = cp.zeros((natoms,maxNb)).astype(int)
    idxNb[(iidx,jidx2)] = jidx+1
    
    nbxyz = cp.zeros((natoms,maxNb,3))
    nbxyz[(iidx,jidx2)] = Rd[idxMask]

    return nNb, idxNb, nbxyz


