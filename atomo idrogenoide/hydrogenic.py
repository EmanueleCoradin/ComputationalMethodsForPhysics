
import numpy as np
import pandas as pd
from itertools import product
import random

#--------- CONSTANTS DECLARATION ----------

ħc  = 1973269804E-15    #eV m
ME = 0.51099895000E06   #eV
MP = 938.27208816E06    #eV
N = 81                #number of grid point for each direction
L = 10                 #grid length in Bohr radius [x1=-L, xN = L]
H = 2*L/(N-1)
LEARNING_RATE = 0.015
EPSILON = 0.01
FAIL = 1

#------------ VARIABLES DECLARATION ------------

guess = np.zeros((N,N,N), dtype=complex)
psi = np.zeros((N,N,N), dtype=complex)
psi0 = np.zeros((N,N,N), dtype=complex)
eigenvalues = []
norm = 0
guess_energy = 0
psi_energy = 0
ground_file = open("ground_state1.txt", "bw")
second_file = open("second_state1.txt", "bw")

#---*------ FUNCTIONS --------------

def inizialization(i, j, k, N, H, alfa):
    if (i==0 or j==0 or k == 0 or i==N-1 or j == N-1 or k == N-1):
        return 0
    else:
        #magnitude of the position vector
        rho = np.sqrt((i-N/2)**2+(j-N/2)**2+(k-N/2)**2)*H
        theta = random.uniform(0, 2*np.pi)
        return rho**2*np.exp(-alfa*rho)*complex(1, 0)

def hamiltonian(psi, H):
    nx, ny, nz = psi.shape
    phi = np.zeros((nx,ny,nz), dtype=complex)
    for i in range(1,nx-1):
        for j in range(1, ny-1):
            for k in range(1, nz-1):
                der = (psi[i+1][j][k] - 2* psi[i][j][k] + psi[i-1][j][k]) + (psi[i][j+1][k] - 2* psi[i][j][k] + psi[i][j-1][k]) + (psi[i][j][k+1] - 2* psi[i][j][k] + psi[i][j][k-1])  
                der/=H**2
                pho = np.sqrt((i-N/2)**2+(j-N/2)**2+(k-N/2)**2)*H
                phi[i][j][k] = -der - 2 * psi[i][j][k]/pho
    return phi

def trapezoidal_integration(psi, H):
    nx, ny, nz = psi.shape
    integral = 0

    for i in range(1, nx-1):
        for j in range(1, ny-1):
            #external surfaces point wheight is 0.5
            integral += (psi[i][j][0] + psi [i][j][N-1] + psi[0][i][j] + psi[N-1][i][j] + psi[i][0][j] + psi[i][N-1][j])/2
            for k in range(1, nz-1):
                #internal points wheight is 1
                integral+=psi[i][j][k]
    integral+= (psi[0][0][0] + psi[0][0][N-1] + psi[0][N-1][0] + psi[N-1][0][0]+psi[N-1][N-1][0] + psi[N-1][0][N-1] + psi[0][N-1][N-1] + psi[N-1][N-1][N-1])/8
    integral*=H**3
    return integral

def trapezoidal_intregration_no_edges(psi, H):
    nx, ny, nz = psi.shape
    integral = 0

    for i in range(1, nx-1):
        for j in range(1, ny-1):
            for k in range(1, nz-1):
                #internal points wheight is 1
                integral+=psi[i][j][k]
    integral*=H**3
    return integral

def L2_norm(psi, H):
    nx, ny, nz = psi.shape
    norm = 0

    for i in range(1, nx-1):
        for j in range(1, ny-1):
            for k in range(1, nz-1):
                #internal points wheight is 1
                norm+= np.abs(psi[i][j][k])**2
    norm*=H**3
    return np.sqrt(norm)

#guess function initialization
for idx in product(*[range(N) for _ in range(3)]):
    i, j, k = idx
    guess[i, j, k] = inizialization(i, j, k, N, H, alfa= 2.5)
guess = guess/L2_norm(guess, H)

phi = hamiltonian(guess, H)
guess_energy = trapezoidal_intregration_no_edges(np.real(np.multiply(np.conjugate(guess) , phi)),H)
print("Initial energy: ", guess_energy)

while True:    
    #steepest descend
    psi = guess + LEARNING_RATE * (-phi + guess_energy * guess)
    psi = psi/L2_norm(psi, H)
    phi = hamiltonian(psi, H)
    psi_energy = trapezoidal_intregration_no_edges(np.real(np.multiply(np.conjugate(psi) , phi)),H)
    
    print("Psi energy: ", psi_energy)

    if np.abs(psi_energy-guess_energy) < EPSILON:
        print("The algorithm converged")
        break
    if np.abs(psi_energy-guess_energy) > FAIL or psi_energy > guess_energy:
        print("The algorithm diverged")
        np.copyto(psi, guess)
        psi_energy = guess_energy
        break
    guess_energy = psi_energy
    np.copyto(guess, psi)

#GROUND STATE FOUND
print("Saving ground state function")
psi.tofile(ground_file)
eigenvalues.append(psi_energy)
np.copyto(psi0, psi)

#-------- SECOND STATE ---------

for idx in product(*[range(N) for _ in range(3)]):
    i, j, k = idx
    guess[i, j, k] = inizialization(i, j, k, N, H, alfa= 2.5)
guess = guess/L2_norm(guess, H)
#Subtracting the ground state component from the vector
guess = guess - trapezoidal_intregration_no_edges(np.multiply(np.conjugate(psi0), guess), H) * psi0
guess = guess/L2_norm(guess, H)

phi = hamiltonian(guess, H)
guess_energy = trapezoidal_intregration_no_edges(np.real(np.multiply(np.conjugate(guess) , phi)),H)
print("Initial energy: ", guess_energy)
phi = phi - trapezoidal_intregration_no_edges(np.multiply(np.conjugate(psi0), phi), H) * psi0

while True:    
    #steepest descend
    psi = guess + LEARNING_RATE * (-phi + guess_energy * guess)
    psi = psi/L2_norm(psi, H)
    phi = hamiltonian(psi, H)
    psi_energy = trapezoidal_intregration_no_edges(np.real(np.multiply(np.conjugate(psi) , phi)),H)
    phi = phi - trapezoidal_intregration_no_edges(np.multiply(np.conjugate(psi0), phi), H) * psi0

    print("Psi energy: ", psi_energy)

    if np.abs(psi_energy-guess_energy) < EPSILON:
        print("The algorithm converged")
        break
    if np.abs(psi_energy-guess_energy) > FAIL or psi_energy > guess_energy:
        print("The algorithm diverged")
        np.copyto(psi, guess)
        psi_energy = guess_energy
        break
    guess_energy = psi_energy
    np.copyto(guess, psi)
print("Second state found")
psi.tofile(second_file)
eigenvalues.append(psi_energy)
np.copyto(psi0, psi)
print("File saved")
#df = pd.DataFrame(psi)
#df.to_csv("ground_state.csv")