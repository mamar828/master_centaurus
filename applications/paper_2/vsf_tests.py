import numpy as np
import graphinglib as gl

from src.tools.statistics.advanced_stats import structure_function


def velocity_structure_function_exact(arr):
    """
    Compute 2D velocity structure function for all pairs of points,
    keeping all exact distances.
    Returns arrays of unique distances and corresponding S(r).
    """
    ny, nx = arr.shape
    distances = []
    sq_diffs = []

    # Loop over all pairs of points
    for y1 in range(ny):
        for x1 in range(nx):
            for y2 in range(ny):
                for x2 in range(nx):
                    dy = y2 - y1
                    dx = x2 - x1
                    r = np.sqrt(dx**2 + dy**2)
                    distances.append(r)
                    # sq_diffs.append((arr[y2, x2] - arr[y1, x1])**2)
                    sq_diffs.append(np.abs(arr[y2, x2] - arr[y1, x1]))

    distances = np.array(distances)
    sq_diffs = np.array(sq_diffs)

    # Find unique distances
    unique_r = np.unique(distances)
    S_r = np.zeros_like(unique_r, dtype=float)

    for i, r_val in enumerate(unique_r):
        mask = distances == r_val
        S_r[i] = np.mean(sq_diffs[mask])
    S_r /= np.var(arr)

    return np.column_stack((unique_r, S_r))

# --- Test array ---
# test_array = np.array([[1, 2],
#                        [7, 8]])
# test_array = np.array([[1, 2, 3],
#                        [4, 5, 6],
#                        [7, 8, 9]])
np.random.seed(0)
test_array = np.random.rand(5, 5) * 10

print(f"σ² = {np.var(test_array)}")

sample_test = velocity_structure_function_exact(test_array)
print(sample_test)

code_test = structure_function(test_array, 1)
print(code_test)
