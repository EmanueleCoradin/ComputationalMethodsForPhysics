
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import numpy as np

N = 81                #number of grid point for each direction
L = 10                 #grid length in Bohr radius [x1=-L, xN = L]
H = 2*L/(N-1)

ground_file = open("second_state1.txt", "rb")

psi = np.fromfile(ground_file, dtype=np.complex)

#visualization

density = np.power(abs(psi), 2)*1000
X, Y, Z = np.mgrid[-L:L:complex(0,N),-L:L:complex(0,N),-L:L:complex(0,N)]

fig = go.Figure(data=go.Volume(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=density.flatten(),
    isomin=0.1,
    isomax=0.8,
    opacity=0.1, # needs to be small to see through all surfaces
    surface_count=3 # needs to be a large number for good volume rendering
    ))
fig.show()

