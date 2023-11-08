
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import numpy as np


NX = 70                #number of radial grid point
NY = 70               #number of theta grid point
NZ = 70                #number of phi grid point

LX = 15.                 #maximun radius in Bohr unit
LY = 15.                 #maximun radius in Bohr unit
LZ = 15.                 #maximun radius in Bohr unit

X = np.logspace(start = np.log(0.2), stop = np.log(LX), num = int(NX/2), base = np.e, endpoint = True, dtype=float)
X = np.concatenate((np.flip(-X), X))
Y = np.logspace(np.log(0.2), np.log(LY), int(NY/2), base = np.e, endpoint = True, dtype=float)
Y = np.concatenate((np.flip(-Y), Y))
Z = np.logspace(np.log(0.2), np.log(LZ), int(NZ/2), base = np.e, endpoint = True, dtype=float)
Z = np.concatenate((np.flip(-Z), Z))

DX = np.zeros(NX-1)
DY = np.zeros(NY-1)
DZ = np.zeros(NZ-1)

for i in range(0, NX-1):
    DX[i] = X[i+1] - X[i]
for i in range(0, NY-1):
    DY[i] = Y[i+1] - Y[i]
for i in range(0, NZ-1):
    DZ[i] = Z[i+1] - Z[i]

ground_file = open("last_fail.npy", "rb")

psi = np.load(ground_file)

density = np.power(abs(psi), 2)
density/=np.max(density)
density = density[15:55,15:55, 15:55]
#print(density)
X, Y, Z = np.mgrid[-LX:LX:complex(0,NX),-LY:LY:complex(0,NY),-LZ:LZ:complex(0,NZ)]

X = X[15:55,15:55, 15:55]
Y = Y[15:55,15:55, 15:55]
Z = Z[15:55,15:55, 15:55]

fig = go.Figure(data=go.Volume(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=density.flatten(),
    isomin=0.2,
    isomax= 0.9,
    opacity=0.1, # needs to be small to see through all surfaces
    surface_count=20# needs to be a large number for good volume rendering
    ))
fig.show()

