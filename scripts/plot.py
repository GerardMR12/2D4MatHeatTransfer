import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load data
j = json.load(open('src\\data.json'))

# Parameters
P1 = j["Geometry"]["p1"]
P2 = j["Geometry"]["p2"]
P3 = j["Geometry"]["p3"]
L = P3[0]
H = P3[1]
N = j["Control Volumes"]["N"]
M = j["Control Volumes"]["M"]
DT = j["Time Parameters"]["dt"]
TI = j["Time Parameters"]["ti"]
TF = j["Time Parameters"]["tf"]
N_T = int((TF-TI)/DT)
CMAP = 'coolwarm'
LEVELS = 15
PLOTDISTR = True

# Load data
path = 'output\\GS N=' + str(N) + ' M=' + str(M) + ' nt=' + str(N_T) + ' results.txt'
data = pd.read_csv(path, sep=r'\s+', skiprows=0, header=None)

# Convert to numpy array
data = data.values

# X axis
X = np.linspace(0, L, N+1)
X[1:] = X[1:] - 0.5*X[1]
X = np.append(X, L)

# Y axis
Y = np.linspace(H, 0, M+1)
Y[1:] = Y[1:] + 0.5*Y[M-1]
Y = np.append(Y, 0)

# Meshgrid
X, Y = np.meshgrid(X, Y)

# Surface plot
fig = plt.figure(figsize=(8, 6))
ax3D = fig.add_subplot(projection='3d')
ax3D.plot_surface(X, Y, data, cmap=CMAP)

# Set the labels
title = 'Resulting field at t = ' + str(int(TF)) + 's'
ax3D.set_title(title, fontsize=18, fontweight='bold', color='red')
ax3D.set_xlabel('X position (m)', fontsize=12, fontweight='bold', color='blue')
ax3D.set_ylabel('Y position (m)', fontsize=12, fontweight='bold', color='blue')
ax3D.set_zlabel('Temperature (ºC)', fontsize=12, fontweight='bold', color='blue')

# Change the numbers of the axes with less values
ax3D.locator_params(axis='x', nbins=5)
ax3D.locator_params(axis='y', nbins=5)
ax3D.locator_params(axis='z', nbins=5)

# Contour plot
fig = plt.figure(figsize=(8, 6))
axC = fig.add_subplot()
contour = axC.contour(X, Y, data, cmap=CMAP, levels=LEVELS)
plt.clabel(contour, colors = 'k', fmt = '%2.1fºC', fontsize=10, inline=True)
contour_filled = axC.contourf(X, Y, data, cmap=CMAP, levels=LEVELS, alpha=0.8)

# Aspect and drawings
axC.set_aspect(1)
plt.plot([P1[0], P1[0]], [0, H], color='black', linestyle='dashed', linewidth=0.5)
plt.plot([P1[0], L], [P2[1], P2[1]], color='black', linestyle='dashed', linewidth=0.5)
plt.plot([0, P1[0]], [P1[1], P1[1]], color='black', linestyle='dashed', linewidth=0.5)

# Set the labels
title = 'Resulting field at t = ' + str(int(TF)) + 's'
axC.set_title(title, fontsize=18, fontweight='bold')
axC.set_xlabel('X position (m)', fontsize=12)
axC.set_ylabel('Y position (m)', fontsize=12)

# Change the numbers of the axes with less values
axC.locator_params(axis='x', nbins=5)
axC.locator_params(axis='y', nbins=5)

# Field plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot()
im = ax.imshow(data, cmap=CMAP)

# Set the labels
ax.set_title(title, fontsize=18, fontweight='bold')
ax.set_xlabel('X position (m)', fontsize=12)
ax.set_ylabel('Y position (m)', fontsize=12)

# Aspect and drawings
ax.set_aspect(H/L)
plt.plot([(N+1)*P1[0]/L, (N + 1)*P1[0]/L], [0, M+1], color='black', linestyle='dashed')
plt.plot([(N+1)*P1[0]/L, N+1], [M+1-(M+2)*P2[1]/H, M+1-(M+2)*P2[1]/H], color='black', linestyle='dashed')
plt.plot([0, (N+1)*P1[0]/L], [M+1-(M+2)*P1[1]/H, M+1-(M+2)*P1[1]/H], color='black', linestyle='dashed')

# Tick labels from 0 to L
ax.set_xticks(np.linspace(0, N+1, 2))
ax.set_xticklabels(np.linspace(0, L, 2))

# Tick labels from 0 to H
ax.set_yticks(np.linspace(0, M+1, 2))
ax.set_yticklabels(np.linspace(H, 0, 2))

# Add a colorbar
fig.colorbar(im, ax=ax, label='Temperature (ºC)', ticks=np.linspace(data.min(), data.max(), 5))

# Plot the mesh
if(PLOTDISTR):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot()
    ones = np.ones((M+2, N+2))
    im = ax.imshow(ones, cmap=CMAP)

    # Set the labels
    ax.set_title('Material distribution', fontsize=18, fontweight='bold')
    ax.set_xlabel('X position (m)', fontsize=12)
    ax.set_ylabel('Y position (m)', fontsize=12)

    # Aspect and drawings
    ax.set_aspect(H/L)
    plt.plot([(N+1)*P1[0]/L, (N + 1)*P1[0]/L], [0, M+1], color='black', linestyle='dashed')
    plt.plot([(N+1)*P1[0]/L, N+1], [M+1-(M+2)*P2[1]/H, M+1-(M+2)*P2[1]/H], color='black', linestyle='dashed')
    plt.plot([0, (N+1)*P1[0]/L], [M+1-(M+2)*P1[1]/H, M+1-(M+2)*P1[1]/H], color='black', linestyle='dashed')

    # Tick labels from 0 to L
    ax.set_xticks(np.linspace(0, N+1, 2))
    ax.set_xticklabels(np.linspace(0, L, 2))

    # Tick labels from 0 to H
    ax.set_yticks(np.linspace(0, M+1, 2))
    ax.set_yticklabels(np.linspace(H, 0, 2))

plt.show()