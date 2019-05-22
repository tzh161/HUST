import numpy as np
import pandas as pd

import os
import time
import sys

#import parameters as pm
from func_w_feat import  MOVEMENT2f_PFVEi,w_feats_file_PFVEi
#from func_struct import get_neighbor_struct_cp



MOVEMENT2f_PFVEi("naphthalene_test.xyz","naphthalene_test.PFV")
#将movement数据集转化为PFV文件#



#res_check = check_PFVEi_Eisum("MOVEMENT_H2O.PFV","potentialEnergy",288,0.1)




#res_rand = rand_Img1(res_check[1],96,96)




w_feats_file_PFVEi("naphthalene_test.PFV","benzene")



