import numpy as np
import graphinglib as gl
import cv2
from astropy.convolution import convolve, Gaussian2DKernel


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
    contour = gl.Contour(*np.mgrid[:image_copy.shape[0], :image_copy.shape[1]][::-1], smoothed_image, **kwargs)
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
    cont_rot = contour.copy()
    # Add 0.5 to convert from pixel center to edge-based coordinates
    x_mesh_edges, y_mesh_edges = contour.x_mesh + 0.5, contour.y_mesh + 0.5
    cont_rot.x_mesh = np.cos(theta) * x_mesh_edges - np.sin(theta) * y_mesh_edges
    cont_rot.y_mesh = np.sin(theta) * x_mesh_edges + np.cos(theta) * y_mesh_edges
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
    hm_rot = heatmap.copy()
    n, m = heatmap.image.shape
    y, x = np.mgrid[:n+1, :m+1]  # grids of each cell x/y corners

    hm_rot._x_coordinates = x * np.cos(theta) - y * np.sin(theta)
    hm_rot._y_coordinates = x * np.sin(theta) + y * np.cos(theta)

    return hm_rot
