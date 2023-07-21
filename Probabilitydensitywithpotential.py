import numpy as np
import matplotlib.pyplot as plt

# Define constants and initial parameters
k = 500
xo = 0.4
sigma2 = 0.01
delta_x = 0.005
delta_t = 5e-7
V = np.zeros(200)
V[120:] = 1e6


# Initialize wave function and arrays
x = np.linspace(0, 1, 200)
psi = np.exp(-(x-xo)**2/(2*sigma2)) * np.exp(1j*k*x)
psi_new = psi.copy()
psi_next = psi.copy()
psi_sq = np.abs(psi)**2

# Define functions for updating the wave function
def update_psi(psi, psi_new, psi_next, V):
    for i in range(1, len(x)-1):
        psi_next[i] = 2*psi[i] - psi_new[i] + (delta_t/delta_x**2)*(psi[i+1] - 2*psi[i] + psi[i-1]) - delta_t*V[i]*psi[i]
    psi_new = psi.copy()
    psi = psi_next.copy()
    return psi, psi_new, psi_next

# Simulate wave function propagation
fig, axs = plt.subplots(4, 1, figsize=(8, 12))
for i in range(5000):
    psi, psi_new, psi_next = update_psi(psi, psi_new, psi_next, V)
    psi_sq = np.abs(psi)**2
    
    # Plot results at different times
    if i == 0:
        axs[0].plot(x, psi_sq)
        axs[0].set_title('t = 0 s')
        axs[0].set_ylim(0, np.max(psi_sq)*1.1)
        axs[0].axvline(x=0.6, linestyle='--', color='gray')
    elif i == 1500:
        axs[1].plot(x, psi_sq)
        axs[1].set_title('t = 7.5e-4 s')
        axs[1].set_ylim(0, np.max(psi_sq)*1.1)
        axs[1].axvline(x=0.6, linestyle='--', color='gray')
    elif i == 2500:
        axs[2].plot(x, psi_sq)
        axs[2].set_title('t = 1.25e-3 s')
        axs[2].set_ylim(0, np.max(psi_sq)*1.1)
        axs[2].axvline(x=0.6, linestyle='--', color='gray')
    elif i == 3500:
        axs[3].plot(x, psi_sq)
        axs[3].set_title('t = 1.75e-3 s')
        axs[3].set_ylim(0, np.max(psi_sq)*1.1)
        axs[3].axvline(x=0.6, linestyle='--', color='gray')
        break

plt.show()
