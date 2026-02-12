import numpy as np
import graphinglib as gl
import cv2
import pvextractor
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.wcs import WCS
from matplotlib.colors import ListedColormap

ROTATION_ANGLE_NIRSPEC = -48  # degrees

def get_smoothed_contour(image: np.ndarray, gaussian_kernel_stddev: float, **kwargs) -> gl.Contour:
    """
    Gives a smoothed Contour object from the input image. inf values in the image will be interpolated using inpainting.
    The contour is created by inverting vertically the image data to match a Heatmap that would be plotted with
    `origin_position="lower"`.

    Parameters
    ----------
    image : np.ndarray
        2D array representing the image data.
    gaussian_kernel_stddev : float
        Standard deviation for the Gaussian kernel used in smoothing.
    **kwargs
        Additional keyword arguments to pass to the Contour constructor.

    Returns
    -------
    gl.Contour
        The smoothed Contour object, with the origin at the lower left.
    """
    # Detect and replace inf values with NaN for convolution
    inf_mask = np.isinf(image)
    image_copy = np.where(inf_mask, np.nan, image)
    smoothed_image = convolve(
        image_copy,
        Gaussian2DKernel(gaussian_kernel_stddev),
        boundary="extend",
        preserve_nan=True,
    )
    # Interpolate NaNs using inpainting
    smoothed_image = cv2.inpaint(
        smoothed_image.astype(np.float32),
        inf_mask.astype(np.uint8),
        inpaintRadius=1,
        flags=cv2.INPAINT_NS,
    )
    contour = gl.Contour(smoothed_image, *np.mgrid[:image_copy.shape[0], :image_copy.shape[1]][::-1], **kwargs)
    return contour

def get_rotated_contour(contour: gl.Contour) -> gl.Contour:
    """
    Returns a rotated copy of the input Contour object by -48 degrees around the origin (0, 0). This allows to plot
    NIRSpec results in a square subplot.

    Parameters
    ----------
    contour : gl.Contour
        The original Contour object to be rotated.

    Returns
    -------
    gl.Contour
        A new Contour object that is the rotated version of the input.
    """
    theta = ROTATION_ANGLE_NIRSPEC * np.pi / 180
    # Add 0.5 to convert from pixel center to edge-based coordinates
    x_mesh_edges, y_mesh_edges = contour.x_mesh + 0.5, contour.y_mesh + 0.5
    cont_rot = contour.copy_with(
        x_mesh=np.cos(theta) * x_mesh_edges - np.sin(theta) * y_mesh_edges,
        y_mesh=np.sin(theta) * x_mesh_edges + np.cos(theta) * y_mesh_edges,
    )
    return cont_rot

def get_rotated_heatmap(heatmap: gl.Heatmap) -> gl.Heatmap:
    """
    Returns a rotated copy of the input Heatmap object by -48 degrees around the origin (0, 0). This allows to plot
    NIRSpec results in a square subplot.

    Parameters
    ----------
    heatmap : gl.Heatmap
        The original Heatmap object to be rotated.

    Returns
    -------
    gl.Heatmap
        A new Heatmap object that is the rotated version of the input.
    """
    theta = ROTATION_ANGLE_NIRSPEC * np.pi / 180
    n, m = heatmap.image.shape
    y, x = np.mgrid[:n+1, :m+1]  # grids of each cell x/y corners

    hm_rot = heatmap.copy_with(
        x_mesh=x * np.cos(theta) - y * np.sin(theta),
        y_mesh=x * np.sin(theta) + y * np.cos(theta),
    )
    return hm_rot

def make_pv_diagram(
    data_cube: np.ndarray,
    wcs: WCS,
    aperture: list[tuple[float, float]] | str,
    width: float = 1.0,
    spacing: float = 1.0,
    contour_levels: list[float] | None = None,
) -> tuple[gl.Polygon, list[gl.Polygon], gl.Arrow, gl.Heatmap, gl.Contour, gl.Contour]:
    """
    Creates a PV diagram from the given aperture path. This function uses the `pvextractor` library to extract the PV
    slice. This library automatically computes the average in each bin along the aperture and weights the pixels by
    their overlap with the aperture.

    Parameters
    ----------
    data_cube : np.ndarray
        The data cube from which to extract the PV diagram.
    wcs : astropy.wcs.WCS
        The WCS associated with the data cube. Note that it may be needed to provide a WCS of the raw cube rather than
        the reduced one, as pvextractor tends to be very picky with WCS objects.
    aperture : list[tuple[float, float]] | str
        This can be either:
        - List of (x, y) tuples defining the aperture path. Currently, only lists with two elements are supported.
        - String path of a DS9 region file of a line. Polygonal regions are not supported. Also note that the line
        region must be saved using either the 'galactic', 'fk5', 'fk4' or 'icrs' coordinate system.
    width : float, default=1.0
        Width of the aperture in pixels.
    spacing : float, default=1.0
        Spacing between samples along the aperture in pixels. For example, a spacing of 2 will create "bins" of 2 pixels
        along the aperture and average the data in each bin.
    contour_levels : list[float], optional
        If provided, the levels for the PV contours. If not provided, the levels are computed automatically.

    Returns
    -------
    tuple[gl.Polygon, list[gl.Polygon], gl.Arrow, gl.Heatmap, gl.Contour, gl.Contour]
        A tuple containing the following elements:
        - Aperture polygon (gl.Polygon)
        - List of bin polygons, showing the region in which each bin is performed (list[gl.Polygon])
        - Arrow indicating the direction of the aperture (gl.Arrow)
        - PV heatmap, constructed from pvextractor output (gl.Heatmap)
        - PV filled contour showing the same data as the PV heatmap (gl.Contour)
        - PV contour showing only the exterior lines of the filled contour (gl.Contour)
    """
    if isinstance(aperture, str):
        path = pvextractor.paths_from_regfile(aperture)[0]
        path.width = width
    else:
        path = pvextractor.Path(aperture, width=width)

    pv_data = pvextractor.extract_pv_slice(data_cube, path, wcs, spacing).data

    # Polygon for the aperture
    path_xy = np.array(path.get_xy(wcs=wcs.celestial))
    angle = np.arctan2(*(path_xy[1] - path_xy[0])[::-1])
    upper_vertices = path_xy + path.width / 2 * np.array([np.sin(angle), -np.cos(angle)])
    lower_vertices = path_xy - path.width / 2 * np.array([np.sin(angle), -np.cos(angle)])
    vertices = np.vstack([upper_vertices, lower_vertices[::-1]])
    aperture_poly = gl.Polygon(vertices, line_width=1.5, fill=False)

    # Polygons for bins
    patches = path.to_patches(spacing, wcs=wcs)
    polygon_vertices = [p.get_xy() for p in patches]
    bin_polygons = [gl.Polygon(v, line_width=0.5, fill=False, edge_color="k") for v in polygon_vertices]

    # Arrow to indicate direction on top of aperture
    arrow_length = 10
    mid_point = (path_xy[0] + path_xy[-1]) / 2
    dir_vector = (path_xy[-1] - path_xy[0])
    dir_vector = dir_vector / np.linalg.norm(dir_vector)
    perp_vector = np.array([-dir_vector[1], dir_vector[0]])
    arrow_start = mid_point - arrow_length / 2 * dir_vector + perp_vector * (path.width / 2 + 1.5)
    arrow_end = arrow_start + dir_vector * arrow_length
    aperture_arrow = gl.Arrow(
        [arrow_start[0], arrow_start[1]],
        [arrow_end[0], arrow_end[1]],
        width=1.5,
        color="black",
        style="->",
    )

    # Heatmap
    pv_hm = gl.Heatmap(pv_data, origin_position="lower")

    # Contours
    meshes = np.meshgrid(np.arange(pv_data.shape[1]), np.arange(pv_data.shape[0]))
    pv_cont = gl.Contour(pv_data, *meshes, levels=contour_levels, color_map="Reds", show_color_bar=False,
                         color_map_range=(contour_levels[0], contour_levels[-1]))
    pv_cont_cont = pv_cont.copy_with(
        color_map=ListedColormap("k"),
        filled=False,
        line_widths=0.5,
    )

    return aperture_poly, bin_polygons, aperture_arrow, pv_hm, pv_cont, pv_cont_cont
