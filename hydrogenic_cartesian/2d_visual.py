import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
from numba import njit
from numba import jit
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings
@jit
def hamiltonian(psi):
    nx, ny, nz = psi.shape
    phi = np.zeros((nx,ny,nz), dtype=complex)
    for i in range(1,nx-1):
        for j in range(1, ny-1):
            for k in range(1, nz-1):
                der = (psi[i+1][j][k] - 2* psi[i][j][k] + psi[i-1][j][k])/DX/DX + (psi[i][j+1][k] - 2* psi[i][j][k] + psi[i][j-1][k])/DY/DY + (psi[i][j][k+1] - 2* psi[i][j][k] + psi[i][j][k-1])/DZ/DZ
                pho = radius(i,j,k)
                phi[i][j][k] = - der - 2. * psi[i][j][k]/pho
    return phi

@njit
def radius(i, j, k):
    global DX, DY, DZ, NX, NY, NZ
    return np.sqrt(((i-NX/2.)*DX)**2 + ((j-NY/2.)*DY)**2 + ((k-NZ/2.)*DZ)**2)

@jit
def vlm_integral(psi):
    global DX, DY, DZ
    nx, ny, nz = psi.shape
    intX = []
    for i in range(0,nx):
        intY = []
        for j in range(0, ny):
            #intgral in dz
            intZ = integrate.simpson(psi[i,j,:], dx = DZ)      #the number of points of the grid must be odd
            intY.append(intZ)
        #computing integral in dy
        intY = np.array(intY)
        integral = integrate.simpson(intY, dx =  DY)
        intX.append(integral)
    #computing integral in dx and return the result
    return integrate.simpson(intX, dx = DX)
   
NX = 151                #number of radial grid point
NY = 151               #number of theta grid point
NZ = 151                #number of phi grid point

LX = 27.                 #maximun radius in Bohr unit
LY = 27.                 #maximun radius in Bohr unit
LZ = 27.                 #maximun radius in Bohr unit

DX = 2.*LX/(NX-1.)        #X-grid [x1=-Lx, xN = Lx]
DY = 2.*LY/(NY-1.)        #Y-grid [y1=-Ly, yN = Ly]
DZ = 2.*LZ/(NZ-1.)        #Z-grid [z1=-Lz, zN = Lz]

ground_file = open("state_25.npy", "rb")

psi = np.load(ground_file)
phi = hamiltonian(psi)


psi_energy = vlm_integral(np.real(np.multiply(np.conjugate(psi) , phi)))


density = np.real(np.power(abs(psi), 2))
print(density.shape)

R2 = np.zeros((NX, NY, NZ))
for i in range(1, NX-1):
    for j in range(1, NY-1):
        for k in range(1, NZ-1):
            R2[i, j, k] = radius(i, j, k)

#extrapolation of the radial component
mr = vlm_integral(np.multiply(density, R2))

print("Psi_energy: ", psi_energy)
print("Mean radius of state: ", mr)

R2 = np.power(R2, 2)
smr = vlm_integral(np.multiply(density, R2))
print("Squared mean radius of state: ", smr)

density = np.multiply(density, R2)
X, Y, Z = np.mgrid[-LX:LX:complex(0,NX),-LY:LY:complex(0,NY),-LZ:LZ:complex(0,NZ)]

X = X[0:151,0:151, 0:151]
Y = Y[0:151,0:151, 0:151]
Z = Z[0:151,0:151, 0:151]
density = density[0:151,0:151, 0:151]

# generate 2 2d grids for the X & Y bounds

# X and Y are bounds, so density should be the value *inside* those bounds.
# Therefore, remove the last value from the density array.
density = density[75, :, :]
Z = Z[ 75, :, :]
Y = Y[ 75, :, :]

z_min, z_max = np.abs(density).min(), np.abs(density).max()

fig, ax = plt.subplots()

print(X.shape)
print(Y.shape)
print(density.shape)

c = ax.pcolormesh(Z, Y, density, cmap='viridis', vmin=z_min, vmax=z_max, norm='linear')

ax.set_title('pcolormesh')
# set the limits of the plot to the limits of the data
#ax.axis([Z.min(), Z.max(), Y.min(), Y.max()])
fig.colorbar(c, ax=ax)

plt.show()