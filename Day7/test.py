import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import random
import math
import time

# Constants
TWO_PI = 2 * math.pi

# Initialize global variables
t = 0.0
nring, nactive = 50, 200
v, p, dt, s = 1.0, 0.18, 10.0, TWO_PI / 6.0
active_r, active_a, cell = [0] * 10000, [0] * 10000, np.zeros((51, 301), dtype=int)

# Initialize the graphics
plt.figure(figsize=(8, 8))
plt.xlim(-nring, nring)
plt.ylim(-nring, nring)
plt.gca().set_aspect('equal', adjustable='box')
scatter = plt.scatter([], [], s=1)

# Function to create a star cluster
def create(r, a):
    global nactive
    if a < 1:
        a += 6 * r
    if random.random() < p and cell[r][a] != 15:
        nactive += 1
        active_a[nactive] = a
        active_r[nactive] = r
        cell[r][a] = 15

# Function to evolve the system
def evolve():
    global nactive
    newactive_r, newactive_a = [0] * 10000, [0] * 10000
    newactive = 0

    for i in range(1, nactive + 1):
        r, a = active_r[i], active_a[i]
        create(r, a + 1)
        create(r, a - 1)
        angle = ((a * s + v * t) / r) % TWO_PI

        if r < nring:
            wt = (v * t / (r + 1)) % TWO_PI
            ap = int(((angle - wt) * (r + 1) / s) % (6 * (r + 1)))
            create(r + 1, ap)
            create(r + 1, ap + 1)

        if r > 1:
            wt = (v * t / (r - 1)) % TWO_PI
            am = int(((angle - wt) * (r - 1) / s) % (6 * (r - 1)))
            create(r - 1, am)
            create(r - 1, am + 1)

    nactive = newactive
    active_r[:nactive + 1] = newactive_r[:nactive + 1]
    active_a[:nactive + 1] = newactive_a[:nactive + 1]

# Function to update the plot
def update(frame):
    global t
    plt.title(f"Time = {t:.1f}")
    x, y = [], []

    for r in range(1, nring + 1):
        for a in range(1, 6 * r + 1):
            if cell[r][a] > 0:
                theta = (a * s + v * t) / r
                x.append(r * math.cos(theta))
                y.append(r * math.sin(theta))
                cell[r][a] -= 1

    t += dt
    scatter.set_offsets(np.column_stack((x, y)))
    return scatter,

# Create the animation
ani = FuncAnimation(plt.gcf(), update, frames=1000, interval=100, blit=True)

plt.show()
