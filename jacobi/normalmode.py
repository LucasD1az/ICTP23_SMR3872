import matplotlib.pyplot as plt
import numpy as np

# Load the data from the file
data = np.loadtxt('modonormal.dat')

# Create a heatmap
plt.imshow(data, cmap='viridis', origin='upper', vmin=0,vmax=500, extent=[0, data.shape[1], 0, data.shape[0]])
plt.colorbar(label='Value')

# Set axis labels (if needed)
plt.xlabel('X Axis')
plt.ylabel('Y Axis')

# Show the heatmap
plt.show()
