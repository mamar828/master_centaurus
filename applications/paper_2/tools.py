import numpy as np
from numpy.typing import ArrayLike


def sample_polyline(vertices: ArrayLike, separation: float) -> np.ndarray:
    """
    Samples points along a polyline defined by a sequence of points, with a specified separation distance.

    Parameters
    ----------
    vertices : ArrayLike
        An array of two-tuples representing the vertices of the polyline.
    separation : float
        The desired separation distance between sampled points.

    Returns
    -------
    np.ndarray
        An array of sampled points along the polyline, with the specified separation distance.
    """
    vertices = np.asarray(vertices, dtype=float)

    # Segment vectors and lengths
    seg = np.diff(vertices, axis=0)
    seg_len = np.linalg.norm(seg, axis=1)

    # Cumulative distance along the line
    cumdist = np.concatenate([[0], np.cumsum(seg_len)])
    total_length = cumdist[-1]

    # Distances where we want points
    d = np.arange(0, total_length, separation)
    if (d[-1] != total_length) and (total_length % separation == 0):
        d = np.append(d, total_length)

    # Find which segment each sample belongs to
    idx = np.searchsorted(cumdist, d, side="right") - 1
    idx = np.clip(idx, 0, len(seg_len)-1)

    # Fraction along each segment
    t = (d - cumdist[idx]) / seg_len[idx]

    # Interpolate
    points = vertices[idx] + t[:, None] * seg[idx]

    return points
