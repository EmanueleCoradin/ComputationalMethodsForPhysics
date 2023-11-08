
import numpy as np
from scipy import integrate
import scipy.special as spe
import math
import sys
#--------- CONSTANTS DECLARATION ----------

ħc  = 1973269804E-15    #eV m
ME = 0.51099895000E06   #eV
MP = 938.27208816E06    #eV

NX = 70                #number of radial grid point
NY = 70               #number of theta grid point
NZ = 70                #number of phi grid point

LX = 15.                 #maximun radius in Bohr unit
LY = 15.                 #maximun radius in Bohr unit
LZ = 15.                 #maximun radius in Bohr unit

LEARNING_RATE = 0.02    #set the learning_rate of the steepest descend method
EPSILON = 0.0001         #precision required for the eigenenergy
LOAD = True             #if True the program will read previous results

#------------ VARIABLES DECLARATION ------------
#grid initialization
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

guess = np.zeros((NX,NY,NZ), dtype=complex)     #array containing the guess function for the steepest descend
psi = np.zeros((NX,NY,NZ), dtype=complex)
psi0 = np.zeros((NX,NY,NZ), dtype=complex)
backup = np.zeros((NX,NY,NZ), dtype=complex)

polynomial = spe.hermite(0) + spe.hermite(1) + spe.hermite(2)

eigenvalues = np.array([], dtype = float)       #array containing the eigenvalues
eigenstates = []                                #list of all the known eigeinstates 
norm = 0
guess_energy = 0
psi_energy = 0
values_file = open("eigenvalues.npy", "rb")     #file where to read/save the eigenvalues
stop = False
fail_count = 0
min_energy = 0

#---------- FUNCTIONS --------------

def psi_ang(phi,theta,l=0,m=0):
    return spe.sph_harm(m,l,phi,theta)
#change this function in order to change the guess function initialization

def radius(i, j, k):
    global X, Y, Z
    return np.sqrt(X[i]**2 + Y[j]**2 + Z[k]**2)

def cartesianToSpherical(i, j, k):
    global X, Y, Z
    r = radius(i, j, k)
    theta = np.arccos(Z[k]/r)
    phi = np.arctan(Y[j]/X[i])

    return r, theta, phi

def initialization(i, j, k, alfa):
    if (i==0 or j==0 or k == 0 or i == NX-1 or j == NY-1 or k == NZ-1):
        return 0. + 0.j
    #magnitude of the position vector
    rho, theta, phi = cartesianToSpherical(i, j, k)

    return (np.exp(-0.1*rho) - np.exp(-1.1*rho))/rho
    '''polynomial(rho) * (np.exp(-rho))*(
    2 + psi_ang(phi, theta, 1, -1) + psi_ang(phi, theta, 1, 0) + psi_ang(phi, theta, 1, 1) 
    + psi_ang(phi, theta, 2, -2) + psi_ang(phi, theta, 2, -1) + psi_ang(phi, theta, 2, 0) + psi_ang(phi, theta, 2, 1) + psi_ang(phi, theta, 2, 2) 
    + psi_ang(phi, theta, 3, 0))
    '''
#function calculating the coefficients in order to compute the laplacian in spherical coordinates

def hamiltonian(psi):
    global DX, DY, DZ
    nx, ny, nz = psi.shape
    phi = np.zeros((nx,ny,nz), dtype=complex)
    for i in range(1,nx-1):
        for j in range(1, ny-1):
            for k in range(1, nz-1):
                #TODO: change the derivative formula
                Xder = (DX[i-1] * psi[i+1][j][k] - (DX[i] + DX[i-1]) * psi[i][j][k] + DX[i] * psi[i-1][j][k]) * 2. / (DX[i] * DX[i-1] * (DX[i] + DX[i-1]))
                Yder = (DY[j-1] * psi[i][j+1][k] - (DY[j] + DY[j-1]) * psi[i][j][k] + DY[j] * psi[i][j-1][k]) * 2. / (DY[j] * DY[j-1] * (DY[j] + DY[j-1]))
                Zder = (DZ[k-1] * psi[i][j][k+1] - (DZ[k] + DZ[k-1]) * psi[i][j][k] + DZ[k] * psi[i][j][k-1]) * 2. / (DZ[k] * DZ[k-1] * (DZ[k] + DZ[k-1]))
                der = Xder + Yder + Zder
                pho = radius(i,j,k)
                phi[i][j][k] = - der - 2. * psi[i][j][k]/pho
    return phi

def vlm_integral(psi):
    global X, Y, Z
    nx, ny, nz = psi.shape
    intX = []
    for i in range(0,nx):
        intY = []
        for j in range(0, ny):
            #intgral in dz
            intZ = integrate.trapz(psi[i, j, :], x = Z)
            intY.append(intZ)
        #computing integral in dy
        intY = np.array(intY)
        integral = integrate.trapz(intY, x =  Y)
        intX.append(integral)
    #computing integral in dx and return the result
    return integrate.trapz(intX, x = X)

def L2_norm(psi):
    nx, ny, nz = psi.shape
    norm = vlm_integral(np.power(np.abs(psi), 2))
    return np.sqrt(norm)


def projection(phi, list):
    for state in list:
        phi = phi - np.multiply(state, vlm_integral(np.multiply(np.conjugate(state), phi)))
    return phi

#-------------------- MAIN -------------------
#load previous results if desired
print("------------ PROGRAM START -------------")
print("initializing the grid")

if LOAD:
    print("Loading calculated eigenstates")
    eigenvalues = np.load(values_file)
    for i in range(0,eigenvalues.size):
        print("E" + str(i) + " = ", eigenvalues[i])
        state_file = open("state_" + str(i) + ".npy", "rb")
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
            psi0[i, j, k] = initialization(i, j, k, alfa= 1.)
psi0 = psi0/L2_norm(psi0)

np.save(f'guess.npy', psi0)

while not stop:
    #Subtracting the ground state component from the vector
    print("Finding a new Eigenstate")
    
    np.copyto(guess, psi0)
    guess = projection(guess, eigenstates)
    guess = guess/L2_norm(guess)
    phi = hamiltonian(guess)
    
    guess_energy = vlm_integral(np.real(np.multiply(np.conjugate(guess) , phi)))
    phi = projection(phi, eigenstates)
    print()
    
    while True:    
        #steepest descend
        print("Guess energy: ",  guess_energy)
        psi = guess + LEARNING_RATE * (guess_energy * guess - phi)
        
        psi = psi/L2_norm(psi)
        phi = hamiltonian(psi)
        psi_energy = vlm_integral(np.real(np.multiply(np.conjugate(psi) , phi)))
        print("Psi energy: ", psi_energy)
        
        if psi_energy < guess_energy:
            if guess_energy - psi_energy < EPSILON:
                print("The algorithm converged")
                break
            phi = projection(phi, eigenstates)
            guess_energy = psi_energy
            np.copyto(guess, psi)
        else:
            print("The algorithm diverged")
            LEARNING_RATE/=2.
            fail_count+=1
            if fail_count > 5:
                np.save(f"last_fail.npy", psi)
                sys.exit(1)
    print()
    
    print("State found: Energy ", psi_energy)    
    eigenvalues = np.append(eigenvalues, psi_energy)
    print("Saving state file: ")
    state_file = open("state_" + str(eigenvalues.size - 1) + ".npy", "wb")
    np.save(state_file, psi)
    
    eigenstates.append(psi)
    print("File saved")
    stop = input("Type q to stop ") == "q"

values_file.close()
values_file = open("eigenvalues.npy", "wb")
np.save(values_file, eigenvalues)