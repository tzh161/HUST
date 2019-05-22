import numpy as np
import matplotlib.pyplot as plt
from sgdml.predict import GDMLPredict
from sgdml.utils import io
# from mpl_toolkits.mplot3d import Axes3D

,_ = io.read_xyz('examples/geometries/aspirin.xyz') # 9 atoms
print r.shape # (1,27)

model = np.load('models/aspirin.npz')
gdml = GDMLPredict(model)
e,f = gdml.predict(r)
print e.shape # (1,)
print f.shape # (1,27)

np.set_printoptions(threshold=np.inf)

ASPRN = np.load('/Users/sGDML/data/aspirin/aspirin_test.npz')
print(ASPRN['E'])
e0 = ASPRN['E']
e0 = np.reshape(e0,(1,-1))
aspirin_delta_e=e-e0
NX = aspirin_delta_e.shape+1 #500+1
x = np.arange(1,NX)
plt.scatter(x,aspirin_delta_e,c = 'r',marker = 'o')
plt.xlabel("Number of Aspirin Molecule Geometries")
plt.ylabel("Error of Energy Prediction")
plt.title("Error of Energy Prediction with Aspirin in sGDML")
plt.show()


ASPRN = np.load('/Users/sGDML/data/aspirin/aspirin_test.npz')
f0 = ASPRN['F']
f0.shape=(500,21,3)
f_0=np.zeros((500,6,3))
for i in range(500):
      for j in range(6):
         f_0[i,j,:]=f0[i,j,:]
f_0.resize(500,63)
f0=f_0
delta_aspirin_f=f-f0
delta_aspirin_f.resize(1,31500)
m=np.arange(1,31501)
plt.scatter(m,aspirin_delta_f,c = 'g',marker = 'o')
plt.xlabel("Number of Aspirin Molecule Geometries in xyz Directions")
plt.ylabel("Error of Force Prediction")
plt.title("Error of Force Prediction with Aspirin in sGDML")
plt.show()