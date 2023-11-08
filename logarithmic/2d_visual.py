import matplotlib.pyplot as plt
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

print(X)

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

density = np.real(np.power(abs(psi), 2))
print(density.shape)
density/=np.max(density)

X, Y, Z = np.mgrid[-LX:LX:complex(0,NX),-LY:LY:complex(0,NY),-LZ:LZ:complex(0,NZ)]

X = X[30:40,30:40,30:40]
Y = Y[30:40,30:40,30:40]
Z = Z[30:40,30:40,30:40]
density = density[30:40,30:40,30:40]

# generate 2 2d grids for the X & Y bounds

# X and Y are bounds, so density should be the value *inside* those bounds.
# Therefore, remove the last value from the density array.
density = density[5, :, :]
Z = Z[5, :, :]
Y = Y[5, :, :]

z_min, z_max = np.abs(density).min(), np.abs(density).max()

fig, ax = plt.subplots()

print(Z.shape)
print(Y.shape)
print(density.shape)

c = ax.pcolormesh(Z, Y, density, cmap='RdBu', vmin=z_min, vmax=z_max)

ax.set_title('pcolormesh')
# set the limits of the plot to the limits of the data
ax.axis([Z.min(), Z.max(), Y.min(), Y.max()])
fig.colorbar(c, ax=ax)

plt.show()