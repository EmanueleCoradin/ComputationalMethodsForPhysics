import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
import scipy.special as spe

f = open("debug.txt", "w")
fh = open("debug_h.txt", "w")
fg = open("debug_g.txt", "w")

def hamiltonian(psi):
    global f
    nx, ny, nz = psi.shape
    phi = np.zeros((nx,ny,nz), dtype=complex)
    for i in range(1,nx-1):
        for j in range(1, ny-1):
            for k in range(1, nz-1):
                pho = radius(i,j,k)
                der = (psi[i+1][j][k] - 2. * psi[i][j][k] + psi[i-1][j][k])/DX/DX + (psi[i][j+1][k] - 2.* psi[i][j][k] + psi[i][j-1][k])/DY/DY + (psi[i][j][k+1] - 2.* psi[i][j][k] + psi[i][j][k-1])/DZ/DZ
                pot = - 2. * psi[i][j][k]/pho
                f.write(str("Cella: " + str(i) + ", " + str(j) + ", " + str(k) + "\n"))
                f.write(str("Kin: " + str(der)+ "\n"))
                f.write(str("Pot: " + str(pot)+ "\n"))
                
                phi[i][j][k] = - der + pot
    return phi

def radius(i, j, k):
    global DX, DY, DZ, NX, NY, NZ
    return np.sqrt(((i-NX/2.)*DX)**2 + ((j-NY/2.)*DY)**2 + ((k-NZ/2.)*DZ)**2)

def cartesianToSpherical(i, j, k):
    global DX, DY, DZ
    x = (i-NX/2.)*DX
    y = (j-NY/2.)*DY
    z = (k-NZ/2.)*DZ

    r = radius(i, j, k)
    theta = np.arccos(z/r)
    phi = np.arctan(y/x) 
    if(x < 0):
        phi += np.pi
    elif(x>0 and y < 0):
        phi += 2.*np.pi
    
    return r, theta, phi

def initialization(i, j, k, alfa):
    global DX, DY, DZ, NX, NY, NZ
    #if (i==0 or j==0 or k == 0 or i == NX-1 or j == NY-1 or k == NZ-1):
    #    return 0.
    #magnitude of the position vector
    rho, theta, phi = cartesianToSpherical(i, j, k)
    return np.exp(-rho/2) *rho* psi_ang(phi, theta, 1,1)

def psi_ang(phi,theta,l=0,m=0):
    return spe.sph_harm(m,l,phi,theta)
    
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
def trpz_integral(psi):
    global DX, DY, DZ
    nx, ny, nz = psi.shape
    integral = 0

    for i in range(1, nx-1):
        for j in range(1, ny-1):
            for k in range(1, nz-1):
                #internal points wheight is 1
                integral+=psi[i, j, k]
    integral*= (DX*DY*DZ)
    return integral
def L2_norm(psi):
    nx, ny, nz = psi.shape
    norm = trpz_integral(np.power(np.abs(psi), 2))
    return np.sqrt(norm)   
NX = 51                #number of radial grid point
NY = 51              #number of theta grid point
NZ = 51               #number of phi grid point

LX = 50.                 #maximun radius in Bohr unit
LY = 50.                 #maximun radius in Bohr unit
LZ = 50.                 #maximun radius in Bohr unit

DX = 2.*LX/(NX-1.)        #X-grid [x1=-Lx, xN = Lx]
DY = 2.*LY/(NY-1.)        #Y-grid [y1=-Ly, yN = Ly]
DZ = 2.*LZ/(NZ-1.)        #Z-grid [z1=-Lz, zN = Lz]

psi = np.zeros((NX, NY, NZ), dtype=complex)
for i in range(0, NX):
    for j in range(0, NY):
        for k in range(0, NZ):
            psi[i, j, k] = initialization(i, j, k, alfa= 1.)
            
psi = psi/L2_norm(psi) 
phi = hamiltonian(psi)

for i in range(0, NX):
    for j in range(0, NY):
        for k in range(0, NZ):
            fh.write(str("Cella: " + str(i) + ", " + str(j) + ", " + str(k) + "\n" + str(phi[i,j,k]) + "\n"))
            fg.write(str("Cella: " + str(i) + ", " + str(j) + ", " + str(k) + "\n" + str(psi[i,j,k]) + "\n"))

np.save(f'phi1.npy', phi)

psi_energy = vlm_integral(np.real(np.multiply(np.conjugate(psi) , phi)))

print("Psi_energy: ", psi_energy)
density = np.real(np.power(abs(phi), 2))
print(density.shape)

R2 = np.zeros((NX, NY, NZ))
for i in range(1, NX-1):
    for j in range(1, NY-1):
        for k in range(1, NZ-1):
            R2[i, j, k] = radius(i, j, k)

#extrapolation of the radial component
mr = vlm_integral(np.multiply(density, R2))
print("Mean radius of state: ", mr)

R2 = np.power(R2, 2)
smr = vlm_integral(np.multiply(density, R2))
print("Squared mean radius of state: ", smr)

#density = np.multiply(density, R2)
X, Y, Z = np.mgrid[-LX:LX:complex(0,NX),-LY:LY:complex(0,NY),-LZ:LZ:complex(0,NZ)]

X = X[:,:, :]
Y = Y[:,:, :]
Z = Z[:,:, :]
density = density[:,:, :]

# generate 2 2d grids for the X & Y bounds

# X and Y are bounds, so density should be the value *inside* those bounds.
# Therefore, remove the last value from the density array.
density = density[:, :, 5]
Z = Z[ :, :, 5]
Y = Y[ :, :, 5]
X = X[ :, :, 5]

z_min, z_max = np.abs(density).min(), np.abs(density).max()

fig, ax = plt.subplots()

print(X.shape)
print(Y.shape)
print(density.shape)

c = ax.pcolormesh(X, Y, density, cmap='viridis')

ax.set_title('pcolormesh')
# set the limits of the plot to the limits of the data
ax.axis([X.min(), X.max(), Y.min(), Y.max()])
fig.colorbar(c, ax=ax)

plt.show()