import numpy as np
import matplotlib.pyplot as plt

# Define constants
hbar = 1.0   # Planck's constant / 2pi
m = 1.0      # particle mass
delta_x = 0.001   # spatial step
delta_t = 0.000000932  # time step

x0=0.4
sigma=0.001
k=500
# Define potential function
def V(x):
    return 0.0

# Define initial wave function
def psi_0(x):
    return np.exp(-(x - x0)**2 / (sigma)) * np.exp(1j *k*x)

# Define function to calculate time derivative of wave function
def dpsi_dt(psi, x, t):
    # Calculate second derivative of wave function using central difference method
    d2psi_dx2 = (psi[2:] - 2 * psi[1:-1] + psi[:-2]) / delta_x**2
    
    # Calculate potential energy
    V_vec = np.vectorize(V)
    V_x = V_vec(x[1:-1])
    
    # Calculate time derivative of wave function
    dpsi_dt = -1j * hbar / (2 * m) * d2psi_dx2 + 1j * V_x / hbar * psi[1:-1]
    
    # Add boundary conditions
    dpsi_dt = np.concatenate(([0], dpsi_dt, [0]))
    
    return dpsi_dt

# Set up grid
x = np.linspace(0, 1,100)

# Set up initial wave function
psi = psi_0(x)

# Initialize time
t = 0.0

# Set up plotting
fig, ax = plt.subplots()

# Propagate wave function in time
for i in range(5):
    # Calculate new wave function using fourth-order Runge-Kutta method
    k1 = delta_t * dpsi_dt(psi, x, t)
    k2 = delta_t * dpsi_dt(psi + k1/2, x, t + delta_t/2)
    k3 = delta_t * dpsi_dt(psi + k2/2, x, t + delta_t/2)
    k4 = delta_t * dpsi_dt(psi + k3, x, t + delta_t)
    psi += 1/6 * (k1 + 2*k2 + 2*k3 + k4)
    
    # Update time
    t += delta_t
    
    # Plot wave function
    ax.plot(x, np.abs(psi)**2, label=f't = {t}')
    
# Set plot
ax.set_xlabel('x')
ax.set_ylabel('$|\psi(x)|^2$')
plt.show()