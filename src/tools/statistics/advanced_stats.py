import numpy as np
import graphinglib as gl
from scipy.optimize import curve_fit
from copy import deepcopy
from uncertainties import ufloat

from src.tools.statistics.stats_library.stats_library import str_func_cpp
from src.tools.statistics.split_normal import SplitNormal


np_sort = lambda arr: arr[np.argsort(arr[:,0])]

def structure_function(data: np.ndarray, order: int) -> np.ndarray:
    """
    Computes the structure function of a 2D array.

    Parameters
    ----------
    data : np.ndarray
        Data from which to compute the structure function.
    order : int
        Order of the structure function to compute. This corresponds to the exponent applied on the pair differences.

    Returns
    -------
    np.ndarray
        Two-dimensional array with every group of three elements representing the lag and its corresponding structure
        function and uncertainty. The returned array is sorted according to the lag value.
    """
    return np_sort(np.array(str_func_cpp(deepcopy(data), order)))

def get_fitted_structure_function_figure(
    data: np.ndarray,
    fit_bounds: tuple[float, float],
    number_of_iterations: int=10000,
) -> gl.SmartFigure:
    """
    Gives the figure of a fitted structure function in the given interval, computing the fit using Monte-Carlo
    uncertainties. The log10 of the data is taken and a linear fit is computed.

    Parameters
    ----------
    data : np.ndarray
        Data from which to compute the structure function. This should be the data outputted by the function "structure
        function".
    fit_bounds : tuple[float, float]
        x interval in which to execute the linear fit. This should exclude the first few points and the points until
        decorrelation, i.e. where the curve is not linear anymore.
    number_of_iterations : int
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
            p0=[0.1, 0.1],
            maxfev=100000
        )[0])
        # m, b = parameters[-1]
        # fig = gl.SmartFigure(
        #     x_lim=(0.9*21.4, 20*21.4),
        #     y_lim=(100, 360),
        #     elements=[gl.Scatter(x_values_fit, y_values_fit), gl.Curve.from_function(
        #         lambda x: b * x**m,
        #         *fit_bounds,
        #         line_width=2,
        #     )],
        #     title=f"$m={m:.3f}, b={b:.3f}$",
        # ).show()

    parameters = np.array(parameters)
    m, b = parameters.mean(axis=0)
    dm, db = parameters.std(axis=0)  # uncertainties on the m and b parameters
    # print(m, dm, b, db)
    slope = ufloat(m, dm)
    fit = gl.Curve.from_function(
        lambda x: b * x**m,
        *fit_bounds,
        color="red",
        label=f"Slope: {slope:.1u}".replace("+/-", " ± "),
        line_width=2,
    )
    # max_error_curve = gl.Curve.from_function(lambda x: (b - db) * x**(m + dm), *fit_bounds, line_width=5)
    # min_error_curve = gl.Curve.from_function(lambda x: (b - db) * x**(m - dm), *fit_bounds, line_width=0)
    # max_error_curve.fill_between_other_curve = min_error_curve
    # max_error_curve.fill_between_bounds = fit_bounds
    # max_error_curve.fill_between_color = "red"

    fig = gl.SmartFigure(elements=[scatter, fit], log_scale_x=True, log_scale_y=True)
    fig.set_visual_params(use_latex=True, font_family="serif")
    return fig
