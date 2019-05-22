import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

np.set_printoptions(threshold=np.inf)

ETHANOL= np.load('/Users/hejianping/sGDML/data/ethanol/ethanol_test.npz')
e0 = ETHANOL['E']
f0 = ETHANOL['F']
prin