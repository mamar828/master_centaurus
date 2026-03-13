import numpy as np
import graphinglib as gl
import cv2
import pvextractor
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.wcs import WCS
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt

from src.config import ROTATION_ANGLE_NIRSPEC, ROTATION_ANGLE_NIRSPEC_DEG
from src.hdu.header import Header
from src.coordinates.celestial_coords import RA, DEC
from src.coordinates.fits_coords import FitsCoords


def get_smoothed_image(image: np.ndarray, gaussian_kernel_stddev: float) -> np.ndarray:
    """
    Gives a smoothed version of the input image using a Gaussian kernel of the specified standard deviation. inf values
    in the image will be interpolated using inpainting with a radius of 1 pixel.

    Parameters
    ----------
    image : np.ndarray
        2D array representing the image data to smooth.
    gaussian_kernel_stddev : float
        Standard deviation for the Gaussian kernel used in smoothing.

    Returns
    -------
    np.ndarray
        The smoothed image, with inf values replaced by inpainted values.
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
    # smoothed_image = cv2.inpaint(
    #     smoothed_image.astype(np.float32),
    #     inf_mask.astype(np.uint8),
    #     inpaintRadius=1,
    #     flags=cv2.INPAINT_NS,
    # )
    return smoothed_image

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
    smoothed_image = get_smoothed_image(image, gaussian_kernel_stddev)
    contour = gl.Contour(smoothed_image, *np.mgrid[:smoothed_image.shape[0], :smoothed_image.shape[1]][::-1], **kwargs)
    return contour

def rotate_coordinates(
    x: np.ndarray,
    y: np.ndarray,
    theta: float = ROTATION_ANGLE_NIRSPEC,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Apply 2D rotation to coordinate arrays.

    Parameters
    ----------
    x : np.ndarray
        X coordinates to rotate.
    y : np.ndarray
        Y coordinates to rotate.
    theta : float, default=ROTATION_ANGLE_NIRSPEC
        Rotation angle in radians.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Rotated (x, y) coordinates.
    """
    x_rot = x * np.cos(theta) - y * np.sin(theta)
    y_rot = x * np.sin(theta) + y * np.cos(theta)
    return x_rot, y_rot

def rotate(element: gl.Contour | gl.Heatmap | gl.Polygon | gl.Arrow) -> gl.Contour | gl.Heatmap | gl.Polygon | gl.Arrow:
    """
    Rotates the input element by -48 degrees around the origin (0, 0). This allows to plot NIRSpec results in a square
    subplot.

    Parameters
    ----------
    element : gl.Contour | gl.Heatmap | gl.Polygon | gl.Arrow
        The original element to be rotated.

    Returns
    -------
    gl.Contour | gl.Heatmap | gl.Polygon | gl.Arrow
        A new element that is the rotated version of the input.
    """
    match (type(element)):
        case gl.Contour:
            # Add 0.5 to convert from pixel center to edge-based coordinates
            if element.x_mesh is None or element.y_mesh is None:
                x_mesh, y_mesh = np.meshgrid(
                    np.arange(element.z_data.shape[1]),
                    np.arange(element.z_data.shape[0]),
                )
                x_mesh_edges, y_mesh_edges = x_mesh + 0.5, y_mesh + 0.5
            else:
                x_mesh_edges, y_mesh_edges = element.x_mesh + 0.5, element.y_mesh + 0.5
            x_rot, y_rot = rotate_coordinates(x_mesh_edges, y_mesh_edges)
            return element.copy_with(x_mesh=x_rot, y_mesh=y_rot)

        case gl.Heatmap:
            n, m = element.image.shape
            y, x = np.mgrid[:n+1, :m+1]  # grids of each cell x/y corners
            x_rot, y_rot = rotate_coordinates(x, y)
            return element.copy_with(x_mesh=x_rot, y_mesh=y_rot)

        case gl.Polygon:
            vertices = element.vertices + 0.5  # Convert from pixel center to edge-based coordinates
            x_rot, y_rot = rotate_coordinates(vertices[:, 0], vertices[:, 1])
            rot_vertices = np.column_stack([x_rot, y_rot])
            return element.copy_with(vertices=rot_vertices)

        case gl.Arrow:
            vertices = np.array([element.pointA, element.pointB]) + 0.5
            x_rot, y_rot = rotate_coordinates(vertices[:, 0], vertices[:, 1])
            return element.copy_with(pointA=[x_rot[0], y_rot[0]], pointB=[x_rot[1], y_rot[1]])

        case gl.Ellipse:
            center = np.array([element.x_center, element.y_center]) + 0.5
            x_rot, y_rot = rotate_coordinates(center[0], center[1])
            return element.copy_with(x_center=x_rot, y_center=y_rot, angle=(element.angle + ROTATION_ANGLE_NIRSPEC_DEG))

        case _:
            raise TypeError(f"Unsupported element type: {type(element)}")

def make_pv_diagram(
    data_cube: np.ndarray,
    wcs: WCS,
    aperture: list[tuple[float, float]] | str,
    width: float = 1.0,
    spacing: float = 1.0,
    contour_levels: list[float] | None = None,
    arrow_length: float = 10.0,
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
    arrow_length : float, default=10.0
        Length of the arrow indicating the direction of the aperture, in pixels.

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
    mid_point = (path_xy[0] + path_xy[-1]) / 2
    dir_vector = (path_xy[-1] - path_xy[0])
    dir_vector = dir_vector / np.linalg.norm(dir_vector)
    perp_vector = - np.array([-dir_vector[1], dir_vector[0]])
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
    if contour_levels is None:
        contour_levels = np.linspace(np.nanmin(pv_data), np.nanmax(pv_data), 10)[1:-1]
    pv_cont_filled = gl.Contour(pv_data, *meshes, levels=contour_levels, color_map="Reds", show_color_bar=False,
                         color_map_range=(contour_levels[0], contour_levels[-1]))
    pv_cont_lines = pv_cont_filled.copy_with(
        color_map=ListedColormap("k"),
        filled=False,
        line_widths=0.5,
    )

    return aperture_poly, bin_polygons, aperture_arrow, pv_hm, pv_cont_filled, pv_cont_lines

def get_N_E_arrows(
        arrow_length: float = 5.0,
        center: tuple[float, float] = (27, -20),
        theta: float = ROTATION_ANGLE_NIRSPEC,
        arrow_offset: float = 0.2,
        text_offset: float = 0.6,
) -> list[gl.Arrow | gl.Text]:
    """
    Gives plottables for the N and E arrows to be plotted to indicate the cardinal directions.

    Parameters
    ----------
    arrow_length : float, default=5.0
        Length of the arrows in pixels.
    center : tuple[float, float], default=(27, -20)
        Center of the arrows in (x, y) coordinates. This is the point from which the arrows start.
    theta : float, default=ROTATION_ANGLE_NIRSPEC
        Rotation angle in radians to apply to the arrows.
    arrow_offset : float, default=0.2
        Offset in pixels to apply to the starting point of the arrows in order to ensure their starting point overlaps.
    text_offset : float, default=0.6
        Additional offset in pixels to apply to the position of the N and E labels, in order to ensure they do not
        overlap with the tip of the arrows.

    Returns
    -------
    list[gl.Arrow | gl.Text]
        A list of plottable elements, including the N and E arrows and their labels.
    """
    arrow_center = np.array(center)
    north_vector = np.array([-np.sin(theta), np.cos(theta)])
    east_vector = np.array([-np.cos(theta), -np.sin(theta)])
    rotated_arrows = [
        gl.Arrow(arrow_center - arrow_offset*north_vector, arrow_center + arrow_length*north_vector, "k", style="->"),
        gl.Arrow(arrow_center - arrow_offset*east_vector, arrow_center + arrow_length*east_vector, "k", style="->"),
        gl.Text(*(arrow_center + north_vector * (arrow_length + text_offset)), r"\textbf{N}", "k", font_size=15),
        gl.Text(*(arrow_center + east_vector * (arrow_length + text_offset)), r"\textbf{E}", "k", font_size=15),
    ]
    return rotated_arrows

def get_wcs_transformed_contours(
    data: np.ndarray,
    source_wcs: WCS,
    target_wcs: WCS,
    levels: list[float],
    **polygon_kwargs
) -> list[gl.Polygon]:
    """
    Create contours from data in its native resolution, then transform the contour coordinates from the source WCS to
    target WCS pixel coordinates. This preserves the detail in the source data while plotting it in the target
    coordinate system. However, this function returns the contour lines as Polygon objects.

    Parameters
    ----------
    data : np.ndarray
        2D array representing the source image data in its native resolution.
    source_wcs : WCS
        WCS object for the source data coordinate system.
    target_wcs : WCS
        WCS object for the target coordinate system.
    levels : list[float]
        Contour levels to compute.
    **polygon_kwargs
        Additional keyword arguments to pass to each Polygon constructor (e.g., edge_color, line_width).

    Returns
    -------
    list[gl.Polygon]
        List of Polygon objects representing the transformed contour lines.
    """
    # Create a temporary figure to extract contour paths at native resolution
    _, ax = plt.subplots()
    y_pix, x_pix = np.mgrid[0:data.shape[0], 0:data.shape[1]]
    cs = ax.contour(x_pix, y_pix, data, levels=levels)
    plt.close()

    contour_polygons = []
    for level_segments in cs.allsegs:
        for vertices in level_segments:
            if len(vertices) < 3:  # Skip degenerate contours
                continue

            # vertices shape: (N, 2) where vertices[:, 0] is x, vertices[:, 1] is y
            # Transform: source pixel -> source world -> target world -> target pixel
            world_coords = source_wcs.pixel_to_world(vertices[:, 0], vertices[:, 1])
            target_x, target_y = target_wcs.world_to_pixel(world_coords)

            # Create polygon from transformed vertices
            transformed_vertices = np.column_stack([target_x, target_y])
            polygon = gl.Polygon(transformed_vertices, **polygon_kwargs)
            contour_polygons.append(polygon)

    return contour_polygons

def get_AGN_pos(
    header: Header,
    rotated: bool = False,
) -> gl.Point:
    """
    Gets the position of the AGN in pixel coordinates as a gl.Point object. This is done by transforming the AGN world
    coordinates (RA, Dec) to pixel coordinates using the WCS information in the header.

    Parameters
    ----------
    header : Header
        The header containing the WCS information and the AGN world coordinates.
    rotated : bool, default=False
        Whether to apply the NIRSpec rotation to the AGN position. If True, the position will be rotated by -48 degrees
        around the origin (0, 0) to match the orientation of the NIRSpec data.

    Returns
    -------
    gl.Point
        A gl.Point object representing the position of the AGN in pixel coordinates as a red cross.
    """
    agn_world_coords = [RA.from_sexagesimal("12:48:49.2609").degrees, DEC.from_sexagesimal("-41:18:39.417").degrees]
    agn_python_coords = header.celestial.world_to_pixel(agn_world_coords)[0]

    if rotated:
        agn_coords = FitsCoords.from_python(*agn_python_coords)
        theta = ROTATION_ANGLE_NIRSPEC
        agn_edge_coords = np.array(agn_coords.data) + 0.5 - 1  # +0.5 for edge coords, -1 to start plotting at (0, 0)
        rot_matrix = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
        agn_coords = rot_matrix @ agn_edge_coords
    else:
        agn_coords = tuple(reversed(agn_python_coords))

    point = gl.Point(*agn_coords, marker_style="x", face_color="red", marker_size=50)
    return point
