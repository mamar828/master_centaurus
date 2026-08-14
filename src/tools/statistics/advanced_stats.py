import numpy as np
import graphinglib as gl
from scipy.optimize import curve_fit
from copy import deepcopy
from uncertainties import ufloat

from src.tools.statistics.stats_library.build.stats_library import str_func_cpp


np_sort = lambda arr: arr[np.argsort(arr[:,0])]

def structure_function(
    data: np.ndarray,
    order: int,
    log_bin_width: float = None,
    bin_start: float = None,
) -> np.ndarray:
    """
    Computes the structure function of a 2D array.

    Parameters
    ----------
    data : np.ndarray
        Data from which to compute the structure function.
    order : int
        Order of the structure function to compute. This corresponds to the exponent applied on the pair differences.
    log_bin_width : float, optional
        Width of the logarithmic bins for regrouping distances, in logspace. If set to None, 0 or negative, no
        logarithmic binning is applied.
    bin_start : float, optional
        Starting point for the logarithmic bins. This is the lower bound of the first bin. If set to None, 0 or
        negative, no logarithmic binning is applied.

    Returns
    -------
    np.ndarray
        Two-dimensional array with every group of three elements representing the lag and its corresponding structure
        function and uncertainty. The returned array is sorted according to the lag value.
    """
    return np_sort(np.array(str_func_cpp(
        deepcopy(data),
        order,
        log_bin_width if log_bin_width is not None else 0.,
        bin_start if bin_start is not None else 0.,
    )))

def get_fitted_structure_function_figure(
    data: np.ndarray,
    fit_bounds: tuple[float, float],
    number_of_iterations: int = 1000,
) -> gl.SmartFigure:
    """
    Gives the figure of a fitted structure function in the given interval, computing the fit using Monte-Carlo
    uncertainties.

    Parameters
    ----------
    data : np.ndarray
        Data from which to compute the structure function. This should be the data outputted by the function
        `structure_function`.
    fit_bounds : tuple[float, float]
        x interval in which to execute the linear fit. This should exclude the first few points and the points until
        decorrelation, i.e. where the curve is not linear anymore.
    number_of_iterations : int, default=1000
        Number of Monte-Carlo iterations to compute the fit uncertainty.

    Returns
    -------
    gl.SmartFigure
        A log-log Figure containing the data points and their uncertainty as well as a linear fit in the given bounds
        with its corresponding fitted slope.
    """
    scatter = gl.Scatter(
        data[:, 0],
        data[:, 1],
        marker_size=3,
        face_color="black",
    )
    scatter.add_errorbars(
        y_error=data[:, 2],
        cap_width=0,
        errorbars_line_width=0.25,
    )

    try:
        # Fit and its uncertainty
        m = (fit_bounds[0] < data[:,0]) & (data[:, 0] < fit_bounds[1])  # generate the fit mask
        x_values_fit = data[m, 0]
        y_values_distributions = np.random.normal(loc=data[m, 1], scale=data[m, 2], size=(number_of_iterations, np.sum(m)))
        parameters = []
        for y_values_fit in y_values_distributions:
            parameters.append(curve_fit(
                f=lambda x, m, b: b * x**m,
                xdata=x_values_fit,
                ydata=y_values_fit,
                p0=[0.5, 20],
                maxfev=100000
            )[0])

        parameters = np.array(parameters)
        m, b = parameters.mean(axis=0)
        dm, db = parameters.std(axis=0)  # uncertainties on the m and b parameters
        slope = ufloat(m, dm)
        fit = gl.Curve.from_function(
            lambda x: b * x**m,
            *fit_bounds,
            color="red",
            label=f"Slope: ${f"{slope:.2u}".replace("+/-", r"\pm")}$",
            line_width=2,
        )
    except Exception as e:
        print(f"Error while fitting the structure function: {e}")
        fit = None

    fig = gl.SmartFigure(elements=[scatter, fit], log_scale_x=True, log_scale_y=True)
    return fig
