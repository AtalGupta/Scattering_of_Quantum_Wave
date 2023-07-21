import numpy as np
import matplotlib.pyplot as plt

# Define the parameters of the Gaussian wave packet
x0 = 0.4   # Mean position
k0 = 500   # Mean wave vector
sigma = 10**(-3/2)   # Width of the wave packet
C = 1

# Define the function for the Gaussian wave packet
def gaussian_wave_packet(x):
    return C*np.exp(-0.5*((x-x0)/sigma)**2 + 1j*k0*x)

# Create an array of x values to plot
x = np.linspace(0.2, 0.6, 1000)

# Evaluate the wave packet function at each value of x
psi = gaussian_wave_packet(x)
psi_squared = [abs(p)**2 for p in psi]

plt.subplot(3, 1, 1) # divide as 2x2, plot top left
plt.plot(x, np.imag(psi),color="red")
plt.xticks([]) 
plt.title('Gaussian quantum wave packet at t=0')
plt.ylabel('Img[$\psi(x)$]')

plt.subplot(3, 1, 2) # divide as 2x2, plot top right
plt.plot(x, np.real(psi),color="blue")
plt.xticks([]) 
plt.ylabel('Re[$\psi(x)$]')

plt.subplot(3, 1, 3) # divide as 2x1, plot bottom


# Plot the real part of the wave packet function

plt.plot(x, (psi_squared),color="green")


# Label the axes and add a title

plt.xlabel('x')
plt.ylabel('psi(x)^2')

plt.show()