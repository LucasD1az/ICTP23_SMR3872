import matplotlib
matplotlib.use('TkAgg')  # Use the TkAgg backend for compatibility with Tkinter

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import os

# Define the number of iterations and the interval between frames
w = 4400
dt = 0.000000025
num_iterations = int(10/w/dt)#10000
frame_interval = int(num_iterations/100)  # Every 10 iterations

# Create a figure and axis for plotting
fig, ax = plt.subplots()

# Function to read data from a file
def read_data(iteration):
    filename = f"laplace{iteration}.dat"
    if os.path.isfile(filename):
        data = np.loadtxt(filename)
        return data
    else:
        return None

# Initialize the plot with data from the first iteration
data = read_data(0)
if data is not None:
    im = ax.imshow(data, cmap='viridis', vmin=-500, vmax=500)

# Function to update the plot for each frame
def update(iteration):
    data = read_data(iteration)
    if data is not None:
        im.set_array(data)
    return im,

# Create the animation
ani = FuncAnimation(fig, update, frames=range(0, num_iterations, frame_interval), blit=True)

# Set up the colorbar with max and min values 100
cbar = plt.colorbar(im)
cbar.set_label("Value")


# Display the animation
plt.show()
