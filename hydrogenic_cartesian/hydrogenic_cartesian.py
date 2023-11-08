
import numpy as np
from scipy import integrate
import scipy.special as spe
import sys
from numba import njit
from numba import jit
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)


#------------------ CONSTANTS DECLARATION  -------------------
#Achtung: they must be odd
NX = 151
NY = 151
NZ = 151

LX = 27.                  #x limit in Bohr unit
LY = 27.                  #y limit in Bohr unit
LZ = 27.                  #z limit in Bohr unit

DX = 2.*LX/(NX-1.)        #X-grid [x1=-Lx, xN = Lx]
DY = 2.*LY/(NY-1.)        #Y-grid [y1=-Ly, yN = Ly]
DZ = 2.*LZ/(NZ-1.)        #Z-grid [z1=-Lz, zN = Lz]

LEARNING_RATE = 0.1       #set the learning_rate of the steepest descend method
EPSILON = 0.000001        #precision required for the eigenenergy
LOAD = True               #if True the program will read previous results

#--------------------- VARIABLES DECLARATION ---------------------

guess = np.zeros((NX,NY,NZ), dtype=complex)     #array containing the guess function for the steepest descend
psi = np.zeros((NX,NY,NZ), dtype=complex)
psi0 = np.zeros((NX,NY,NZ), dtype=complex)

eigenvalues = np.array([], dtype = float)       #array containing the eigenvalues
eigenstates = []                                #list of all the known eigeinstates 
norm = 0
guess_energy = 0
psi_energy = 0
values_file = open("eigenvalues.npy", "rb")     #file where to read/save the eigenvalues
stop = False
fail_count = 0
min_energy = 0
count = 0

#--------------------- FUNCTIONS ---------------------

def psi_ang(phi,theta,l=0,m=0):
    return spe.sph_harm(m,l,phi,theta)

@njit
def radius(i, j, k):
    global DX, DY, DZ, NX, NY, NZ
    return np.sqrt(((i-NX/2.)*DX)**2 + ((j-NY/2.)*DY)**2 + ((k-NZ/2.)*DZ)**2)

@njit
def cartesianToSpherical(i, j, k):
    global DX, DY, DZ
    x = (i-NX/2.)*DX
    y = (j-NY/2.)*DY
    z = (k-NZ/2.)*DZ

    r = radius(i, j, k)
    theta = np.arccos(z/r)
    phi = np.arctan(y/x)
    #in this frame of reference phi € [0,2pi]
    if(x < 0):
        phi += np.pi
    elif(x > 0 and y < 0):
        phi += 2.*np.pi
    return r, theta, phi

#change this function in order to change the guess function initialization
@jit
def initialization(i, j, k):
    global DX, DY, DZ, NX, NY, NZ
    if (i==0 or j==0 or k == 0 or i == NX-1 or j == NY-1 or k == NZ-1):
        return 0
    #magnitude of the position vector
    rho, theta, phi = cartesianToSpherical(i, j, k)

    return 1/rho * psi_ang(phi, theta, 3, 3)

#hamiltonian functional application on psi: phi = H psi
@jit(cache=True)
def hamiltonian(psi):
    nx, ny, nz = psi.shape
    phi = np.zeros((nx,ny,nz), dtype=complex)
    for i in range(1,nx-1):
        for j in range(1, ny-1):
            for k in range(1, nz-1):
                der = (psi[i+1][j][k] - 2* psi[i][j][k] + psi[i-1][j][k])/(DX**2.) + (psi[i][j+1][k] - 2* psi[i][j][k] + psi[i][j-1][k])(DY**2.) + (psi[i][j][k+1] - 2* psi[i][j][k] + psi[i][j][k-1])/(DZ**2.)
                rho = radius(i,j,k)
                phi[i][j][k] = - der - 2. * psi[i][j][k]/rho
    return phi

#3d integration with trapezoids
@jit(cache=True)
def trpz_integral(psi):
    global DX, DY, DZ
    nx, ny, nz = psi.shape
    integral = 0
    #surface and corner points are exluded because of the border condition
    for i in range(1, nx-1):
        for j in range(1, ny-1):
            for k in range(1, nz-1):
                #internal points wheight is 1
                integral+=psi[i, j, k]
    integral*= (DX*DY*DZ)
    return integral

#3d integration with simpson rule
@jit(cache=True)
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

@jit(cache=True)
def L2_norm(psi):
    norm = trpz_integral(np.power(np.abs(psi), 2))
    return np.sqrt(norm)

#apply the projector of the precomputed eigenstates
@jit(cache=True)
def projection(phi, list):
    for state in list:
        phi = phi - np.multiply(state, trpz_integral(np.multiply(np.conjugate(state), phi)))
    return phi

#--------------------- MAIN ---------------------

#load previous results if desired
print("------------ PROGRAM START -------------")
if LOAD:
    print("Loading calculated eigenstates")
    eigenvalues = np.load(values_file)
    for i in range(0,eigenvalues.size):
        print("E" + str(i) + " = ", eigenvalues[i])
        state_file = open("state_" + str(i+1) + ".npy", "rb")
        state = np.load(state_file)
        state = np.array(state)
        eigenstates.append(state)
        state_file.close()
    print("Terminated succesfully")
    print()

print("Starting algorithm")
print("First guess initialization")
for i in range(0, NX):
    for j in range(0, NY):
        for k in range(0, NZ):
            psi0[i, j, k] = initialization(i, j, k)
psi0 = psi0/L2_norm(psi0)
#np.save(f'guess.npy', psi0)

while not stop:
    print("Finding new Eigenstate")
    
    np.copyto(guess, psi0)
    guess = projection(guess, eigenstates)
    guess = guess/L2_norm(guess)
    phi = hamiltonian(guess)
    
    guess_energy = trpz_integral(np.real(np.multiply(np.conjugate(guess) , phi)))
    phi = projection(phi, eigenstates)
    print()
    print("Guess energy: ",  guess_energy)

    while True:    
        #steepest descend
        psi = guess + LEARNING_RATE * (guess_energy * guess - phi)
        
        psi = psi/L2_norm(psi)
        phi = hamiltonian(psi)
        psi_energy = trpz_integral(np.real(np.multiply(np.conjugate(psi) , phi)))
        print("Psi energy: ", count, psi_energy)
        count+=1

        if psi_energy < guess_energy:
            if guess_energy - psi_energy < EPSILON:
                print("The algorithm converged")
                break
            phi = projection(phi, eigenstates)
            guess_energy = psi_energy
            np.copyto(guess, psi)
        else:
            print("The algorithm diverged")
            if fail_count > 8:
                np.save(f"last_fail.npy", psi)
                sys.exit(1)
            LEARNING_RATE/=1.3
            fail_count+=1
            phi = hamiltonian(guess)
            
    print()
    
    print("State found: Energy ", psi_energy)    
    eigenvalues = np.append(eigenvalues, psi_energy)
    print("Saving state file: ")
    state_file = open("state_"+str(eigenvalues.size)+".npy", "wb")
    
    np.save(state_file, psi)
    
    eigenstates.append(psi)
    print("File saved")
    stop = input("Type q to stop ") == "q"
    count = 0

values_file.close()
values_file = open("eigenvalues.npy", "wb")
np.save(values_file, eigenvalues)

