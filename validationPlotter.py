
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def calculate_plane_points(matrix):
    # Extract x, y, and z coordinates from the matrix
    x = matrix[:, 0]
    y = matrix[:, 1]
    z = matrix[:, 2]

    # Calculate the coefficients of the plane equation
    A = np.linalg.det([[y[0], z[0], 1], [y[1], z[1], 1], [y[2], z[2], 1]])
    B = -np.linalg.det([[x[0], z[0], 1], [x[1], z[1], 1], [x[2], z[2], 1]])
    C = np.linalg.det([[x[0], y[0], 1], [x[1], y[1], 1], [x[2], y[2], 1]])
    D = -np.linalg.det([[x[0], y[0], z[0]], [x[1], y[1], z[1]], [x[2], y[2], z[2]]])

    # Generate points on the plane
    num_points = 5
    x_range = np.linspace(min(x), max(x), num_points)
    y_range = np.linspace(min(y), max(y), num_points)
    X, Y = np.meshgrid(x_range, y_range)
    Z = (-A * X - B * Y - D) / C

    return X, Y, Z



matrix = np.array([[1, 1, 1,4], [1, 2, 3,5], [2, 1, 1,5], [2, 2, 3,6]])

# Extract x, y, and z coordinates from the matrix
x = matrix[:, 0]
y = matrix[:, 1]
z = matrix[:, 2]
z2 = matrix[:, 3]

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the surface
ax.plot_trisurf(x, y, z)
ax.plot_trisurf(x, y, z2)


# plot points on the first plane
plane1 = np.delete(matrix, 3, 1)
plane2 = np.delete(matrix, 2, 1)
X, Y, Z = calculate_plane_points(plane1)
X, Y, Z2 = calculate_plane_points(plane2)

# save the points
points = np.array([X.flatten(), Y.flatten(), Z.flatten()]).T
points2 = np.array([X.flatten(), Y.flatten(), Z2.flatten()]).T

# save the points to a file
np.savetxt('points.txt', points, delimiter=',')
np.savetxt('points2.txt', points2, delimiter=',')

ax.scatter(X, Y, Z, color='r')
ax.scatter(X, Y, Z2, color='g')

# Set labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Surface Plot')

# Show the plot
plt.show()
