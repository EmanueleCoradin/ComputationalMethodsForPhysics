
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import numpy as np
import os, sys
import stat
import plotly.express as px
'''

import plotly.io as pio

pio.renderers.default ="vscode"
'''
def radius(i, j, k):
    global DX, DY, DZ, NX, NY, NZ
    return np.sqrt(((i-NX/2.)*DX)**2 + ((j-NY/2.)*DY)**2 + ((k-NZ/2.)*DZ)**2)


NX = 151                #number of radial grid point
NY = 151               #number of theta grid point
NZ = 151                #number of phi grid point

LX = 27.                 #maximun radius in Bohr unit
LY = 27.                 #maximun radius in Bohr unit
LZ = 27.                 #maximun radius in Bohr unit

DX = 2.*LX/(NX-1.)        #X-grid [x1=-Lx, xN = Lx]
DY = 2.*LY/(NY-1.)        #Y-grid [y1=-Ly, yN = Ly]
DZ = 2.*LZ/(NZ-1.)        #Z-grid [z1=-Lz, zN = Lz]

ground_file = open("state_30.npy", "rb")

psi = np.load(ground_file)

density = np.power(abs(psi), 2)

#density = density[::2,::2,::2]
density = density[0:150:4,0:150:4,0:150:4]
density/=np.max(density)
#print(density)
X, Y, Z = np.mgrid[-LX:LX:complex(0,NX),-LY:LY:complex(0,NY),-LZ:LZ:complex(0,NZ)]

X = X[0:150:4,0:150:4,0:150:4]
Y = Y[0:150:4,0:150:4,0:150:4]
Z = Z[0:150:4,0:150:4,0:150:4]

fig = go.Figure(data=go.Volume(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=density.flatten(),
    isomin=0.1,
    isomax= 0.9,
    opacity=0.1, # needs to be small to see through all surfaces
    surface_count=30# needs to be a large number for good volume rendering
    ))
#fig.show()
#os.chmod("html\state_8.html", stat.S_IWGRP | stat.S_IRGRP | stat.S_IWUSR | stat.S_IRUSR)
fig.write_html("state_proof.html")
