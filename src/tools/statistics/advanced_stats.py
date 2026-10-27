import numpy as np
import graphinglib as gl
import pwlf
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
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
        x interval within which to execute the linear fit. General consensus is that the fit should be done from the
        smallest resolvable scale (e.g. PSF) to half the size of the window within which the structure function is
        computed (e.g. the radius of the circular window).
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
        m = (fit_bounds[0] < data[:, 0]) & (data[:, 0] < fit_bounds[1])  # generate the fit mask
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
        slope = parameters.mean(axis=0)[0]
        slope_uncertainty = parameters.std(axis=0)[0]

        amplitude = curve_fit(
            f=lambda x, b: b * x**slope,
            xdata=x_values_fit,
            ydata=y_values_fit,
            p0=[20],
            maxfev=100000
        )[0]

        fit = gl.Curve.from_function(
            lambda x: amplitude * x**slope,
            *fit_bounds,
            color="red",
            label=f"Slope: ${f"{ufloat(slope, slope_uncertainty):.2u}".replace("+/-", r"\pm")}$",
            line_width=2,
        )

    except Exception as e:
        print(f"Error while fitting the structure function: {e}")
        fit = None

    fig = gl.SmartFigure(elements=[scatter, fit], log_scale_x=True, log_scale_y=True)
    return fig

def piecewise_linear(x: np.ndarray, breakpoint: float, y0: float, m1: float, m2: float) -> np.ndarray:
    """"
    Defines a piecewise linear function with two segments. This is a function that can be used for VSF curve fitting.

    Parameters
    ----------
    x : np.ndarray
        The input values for which to compute the piecewise linear function.
    breakpoint : float
        The x-value at which the function transitions from the first linear segment to the second.
    y0 : float
        The y-value of the function at the breakpoint.
    m1 : float
        The slope of the first linear segment (for x < breakpoint).
    m2 : float
        The slope of the second linear segment (for x >= breakpoint).

    Returns
    -------
    np.ndarray
        The computed y-values of the piecewise linear function corresponding to the input x-values.
    """
    return np.piecewise(
        x,
        [x < breakpoint, x >= breakpoint],
        [lambda x: m1 * (x - breakpoint) + y0, lambda x: m2 * (x - breakpoint) + y0],
    )

def fit_vsf_piecewise(
    data: np.ndarray,
    fit_bounds: tuple[float, float] = None,
    number_of_iterations: int = 1000,
    return_figure: bool = False,
) -> ...:
    """
    Fits a velocity structure function (VSF) to the given data using a piecewise linear model. The function performs a
    Monte-Carlo estimation of the uncertainties on the fitted parameters.

    Parameters
    ----------
    data : np.ndarray
        Data from which to compute the VSF parameters. This should be the data outputted by the function
        `structure_function`.

    Returns
    -------
    np.ndarray
        A two-dimensional array containing the mean and standard deviation (uncertainty) of the fitted parameters. The
        first axis corresponds to the parameters (i.e. breakpoint, y0, m1 and m2 from the `piecewise_linear` function)
        and the second axis corresponds to the mean and standard deviation (uncertainty) of each parameter.
    """
    m = (fit_bounds[0] < data[:, 0]) & (data[:, 0] < fit_bounds[1])  # generate the fit mask

    # Monte-Carlo estimation of the uncertainties on the parameters
    x_data = np.log10(data[m, 0])
    y_data = np.log10(np.random.normal(loc=data[m, 1], scale=data[m, 2], size=(number_of_iterations, np.sum(m))))
    parameters = []
    for y_data_i in y_data:
        parameters.append(curve_fit(
            piecewise_linear,
            x_data,
            y_data_i,
            p0=[0.7, 0.9, 0.9, 0.6],
            maxfev=int(1e6)
        )[0])

    mean_parameters = np.mean(parameters, axis=0)
    std_parameters = np.std(parameters, axis=0)
    stacked_parameters = np.column_stack([mean_parameters, std_parameters])

    if return_figure:
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

        x_predict = np.linspace(x_data[0], x_data[-1])
        y_predict = piecewise_linear(x_predict, *mean_parameters)
        prediction_curve = gl.Curve(
            10**x_predict,
            10**y_predict,
            color="red",
            line_width=2
        )
        fig = gl.SmartFigure(
            log_scale_x=True,
            log_scale_y=True,
            elements=[scatter, prediction_curve],
        )
        return stacked_parameters, fig

    else:
        return stacked_parameters
