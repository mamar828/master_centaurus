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
        with its corresponding equation.
    """
    logged_data = np.log10(data)
    scatter = gl.Scatter(
        logged_data[:,0],
        logged_data[:,1],
        marker_size=3,
        face_color="black",
    )
    # Uncertainties are given in the order left, right or bottom, top
    uncertainties = np.array([
        np.abs(logged_data[:,1] - np.log10(data[:,1] - data[:,2])),
        np.abs(np.log10(data[:,1] + data[:,2]) - logged_data[:,1]),
    ])
    scatter.add_errorbars(
        y_error=uncertainties,
        cap_width=0,
        errorbars_line_width=0.25,
    )

    # Fit and its uncertainty
    m = (fit_bounds[0] < logged_data[:,0]) & (logged_data[:,0] < fit_bounds[1])     # generate the fit mask
    data_distributions = [SplitNormal(loc, *u) for loc, u in zip(logged_data[m,1], uncertainties.T[m])]
    values = np.array([sn.random(number_of_iterations) for sn in data_distributions]).T
    parameters = []
    for val in values:
        parameters.append(curve_fit(
            f=lambda x, m, b: m*x + b,
            xdata=logged_data[m,0],
            ydata=val,
            p0=[0.1,0.1],
            maxfev=100000
        )[0])

    parameters = np.array(parameters)
    m, b = parameters.mean(axis=0)
    dm, db = parameters.std(axis=0)     # uncertainties on the m and b parameters
    slope = ufloat(m, dm)
    fit = gl.Curve.from_function(
        lambda x: m*x + b,
        *fit_bounds,
        color="red",
        label=f"Slope : {slope:.1u}".replace("+/-", " ± "),
        line_width=1,
        number_of_points=2
    )
    y_fit_errors = db + np.array(fit_bounds)*dm
    fit.add_error_curves(
        y_error=y_fit_errors,
        error_curves_line_style=""
    )

    fig = gl.SmartFigure(x_lim=(0,1.35), y_lim=(-0.2,0.5), elements=[scatter, fit])
    fig.set_visual_params(use_latex=True, font_family="serif")
    return fig
