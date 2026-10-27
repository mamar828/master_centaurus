import numpy as np
import graphinglib as gl
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io.fits.header import Header
from re import compile as re_compile


def read_ds9_segments(filename: str) -> tuple[list[np.ndarray], list[str]]:
    """
    Read a DS9 region file that contains `segment` elements.

    Parameters
    ----------
    filename : str
        Path to the DS9 region file.

    Returns
    -------
    tuple[list[np.ndarray], list[str]]
        A tuple containing two lists:
        - A list of numpy arrays, where each array contains the coordinates of a segment. Each array has shape (N, 2),
        where N is the number of vertices in the segment, and the two columns correspond to the RA and Dec coordinates
        of the vertices.
        - A list of strings for the color of each segment.
    """
    segment_pattern = re_compile(r"segment\s*\(([^)]*)\)", flags=0)
    color_pattern = re_compile(r"\bcolor\s*=\s*([a-zA-Z]+)", flags=0)

    segments = []
    colors = []

    with open(filename, "r") as f:
        lines = f.readlines()

    for line in lines:
        segment_match = segment_pattern.search(line)

        if segment_match is None:
            continue

        values = np.fromstring(segment_match.group(1), sep=",")
        segments.append(values.reshape(-1, 2))

        color_match = color_pattern.search(line)
        if color_match is None:
            colors.append("black")
        else:
            colors.append(color_match.group(1).lower())

    if not segments:
        raise ValueError(f"No DS9 segment records found in {filename}")

    return segments, colors


def project_segments_to_pixels(segments: list[np.ndarray], header: Header) -> list[np.ndarray]:
    """
    Project DS9 FK5 coordinates to zero-based image pixel coordinates.

    Parameters
    ----------
    segments : list[np.ndarray]
        A list of numpy arrays, where each array contains the coordinates of a segment in FK5 coordinates. Each array
        has shape (N, 2), where N is the number of vertices in the segment, and the two columns correspond to the RA and
        Dec coordinates of the vertices. This should be the first output of the `read_ds9_segments` function.
    header : Header
        The header of the image to which the segments should be projected.

    Returns
    -------
    list[np.ndarray]
        A list of one array of shape (N, 2) per segment, with columns [x_pixel, y_pixel].
    """
    wcs = WCS(header).celestial
    projected_segments = []

    for segment in segments:
        sky_coordinates = SkyCoord(
            ra=segment[:, 0] * u.deg,
            dec=segment[:, 1] * u.deg,
            frame="fk5",
        )

        x_pixel, y_pixel = wcs.world_to_pixel(sky_coordinates)

        projected_segments.append(
            np.column_stack((x_pixel, y_pixel))
        )

    return projected_segments

def plot_segments(segments: list[np.ndarray], colors: list[str]) -> list[gl.Line]:
    """
    Convert line segments to graphinglib Line objects.

    Parameters
    ----------
    segments : list[np.ndarray]
        A list of one array of shape (N, 2) per segment, with columns [x_pixel, y_pixel]. This should be the output of
        the `project_segments_to_pixels` function.
    colors : list[str]
        A list of strings for the color of each segment.

    Returns
    -------
    list[gl.Line]
        A list of graphinglib Line objects corresponding to the input segments.
    """
    gl_lines = []

    for coordinates, color in zip(segments, colors):
        for start, end in zip(coordinates[:-1], coordinates[1:]):
            gl_lines.append(
                gl.Line(
                    start,
                    end,
                    color=color,
                    width=1.5,
                )
            )

    return gl_lines

def sample_polyline(vertices: np.ndarray, separation: float) -> np.ndarray:
    """
    Samples points along a polyline defined by a sequence of points, with a specified separation distance.

    Parameters
    ----------
    vertices : np.ndarray
        An array of two-tuples representing the vertices of the polyline.
    separation : float
        The desired separation distance between sampled points.

    Returns
    -------
    np.ndarray
        An array of sampled points along the polyline, with the specified separation distance.
    """
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
