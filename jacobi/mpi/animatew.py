import matplotlib
matplotlib.use('TkAgg')  # Use the TkAgg backend for compatibility with Tkinter

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import os

# Define the number of iterations and the interval between frames
#w = 6000
#dt = 0.000000025
num_iterations = 10000
frame_interval = 100  # Every 10 iterations

# Create a figure and axis for plotting
fig, ax = plt.subplots()

# Function to read data from a file
def read_data(iteration):
    # w=int(iteration*np.pi*2)
    # formatted_w = "{:.6f}".format(w)
    print(iteration)
    filename = f"solution_{iteration}.dat"
    if os.path.isfile(filename):
        data = np.loadtxt(filename)
        return data
    else:
        return None

# Initialize the plot with data from the first iteration
data = read_data(4000)
if data is not None:
    im = ax.imshow(data, cmap='viridis', vmin=0, vmax=150)
    text=ax.text(0.1, 0.9, '', transform=ax.transAxes, color='black', fontsize=12)

    # title = ax.set_title('')
    #ax.set_title(f'Frequency (w): {4000}')

# Function to update the plot for each frame
def update(iteration):
    data = read_data(iteration)
    if data is not None:
        im.set_array(data)
        # ax.text(-0.1, 0.9, f'{iteration}', transform=ax.transAxes, color='black', fontsize=12)
        text.set_text(f'Frequency (w): {iteration}')
        # title.set_text(f'Frequency (w): {iteration}')
        # ax.set_title(f'Frequency (w): {iteration}')
    return im, text,

# Create the animation
ani = FuncAnimation(fig, update, frames=range(4000, num_iterations, frame_interval), blit=True)

# Set up the colorbar with max and min values 100
cbar = plt.colorbar(im)
cbar.set_label("Value")


# Display the animation
plt.show()
