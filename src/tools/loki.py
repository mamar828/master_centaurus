import graphinglib as gl
import numpy as np
from astropy.io.fits import open as fits_open
from astropy.constants import c as light_speed

from src.tools.miscellaneous import get_pdf_image_as_array


def get_loki_line_pdfs_figure(
    lines: list[str],
    folder_name: str = "output_NGC4696_G235H_F170LP_full_model",
) -> gl.SmartFigure:
    """
    Gives a figure containing the Loki PDF plots for the specified lines. The figure has 3 columns representing flux,
    velocity offset, and FWHM. These maps are extracted from the "1.flux", "1.voff" and "1.fwhm" files.

    Parameters
    ----------
    lines : list[str]
        The list of spectral lines to include in the figure. These should correspond to Loki output directories in the
        param_maps/lines folder (e.g., "H210_O2").
    folder_name : str, default="output_NGC4696_G235H_F170LP_full_model"
        The folder inside the data/loki directory containing the Loki output.

    Returns
    -------
    gl.SmartFigure
        The figure containing the Loki PDF plots.
    """
    hms = []
    subtitles = []
    for line in lines:
        # for suffix in ["1.flux", "1.voff", "1.fwhm"]:
        for suffix in ["1.flux", "vpeak", "1.fwhm"]:
            hms.append(gl.Heatmap(
                get_pdf_image_as_array(f"data/loki/{folder_name}/param_maps/lines/{line}/{line}.{suffix}.pdf"),
                show_color_bar=False,
            ))
            subtitles.append(f"{line} {suffix.lstrip("1.")}")

    num_rows = len(hms) // 3
    fig = gl.SmartFigure(
        num_rows,
        3,
        elements=hms,
        subtitles=subtitles,
        size=(10, 2.75*num_rows),
        remove_axes=True,
        reference_labels=False,
        aspect_ratio=1,
    )
    return fig

def get_loki_fit_figure(
    model_filename: str,
    spaxel_coordinates: tuple[int, int],
    version: int = 2,
) -> gl.SmartFigure:
    """
    Gives a figure showing the Loki fit for a given spaxel using the specified model file.

    Parameters
    ----------
    model_filename : str
        The path to the Loki model FITS file.
    spaxel_coordinates : tuple[int, int]
        The (x, y) coordinates of the spaxel to visualize.

        .. note::
            To be consistent with Loki's output, the coordinates should start at (1, 1) for the bottom-left spaxel and
            be ordered as (x, y).
    version : int, default=2
        The iteration version of the results used. 2 is for the late october version and 3 is for the december version
        with OQBr tied. 1 is not implemented.

    Returns
    -------
    gl.SmartFigure
        The figure showing the Loki fit for the specified spaxel as well as the data itself.
    """
    if version not in [2, 3]:
        raise ValueError("Only versions 2 and 3 are supported.")
    hdu = fits_open(model_filename)
    data = hdu[1].data
    wavelength_arange = hdu[-1].data[0][0].flatten() / (1 + 0.0099)

    # Building the stellar continuum
    stellar_extinction = hdu[4].data
    if version == 2:
        raw_stellar_continuum = hdu[7].data
        gas_lines_range = range(8, 28)
    else:
        raw_stellar_continuum = hdu[8].data
        polynomials_multiplicative = hdu[7].data
        raw_stellar_continuum *= polynomials_multiplicative
        gas_lines_range = range(9, 29)
    stellar_continuum = raw_stellar_continuum * stellar_extinction

    # Building the gas emission lines
    gas_extinction = hdu[5].data
    silicates_extinction = hdu[6].data
    gas_lines = np.sum([hdu[i].data for i in gas_lines_range], axis=0) * gas_extinction * silicates_extinction

    total_model = stellar_continuum + gas_lines

    # Convert to erg/s/cm²/sr/Hz
    hertz_conversion_factor = light_speed.to("micron/s").value / wavelength_arange
    data *= hertz_conversion_factor[:, None, None]
    total_model *= hertz_conversion_factor[:, None, None]
    stellar_continuum *= hertz_conversion_factor[:, None, None]

    # Building the emission line labels and texts
    lines = [2.2235, 2.1218, 2.0338, 1.9576, 1.8920, 1.8358, 1.7880, 1.7480, 1.7147, 2.4756, 2.5001 , 2.52802, 2.55985,
             2.62688, 2.80251, 3.00387, 1.8745, 2.6259, 2.1661, 1.9451]
    names = ["S(0)", "S(1)", "S(2)", "S(3)", "S(4)", "S(5)", "S(6)", "S(7)", "S(8)", "Q(6)", "Q(7)", "Q(8)", "Q(9)",
             "O(2)", "O(3)", "O(4)", r"Pa$\alpha$", r"Br$\beta$", r"Br$\gamma$", r"Br$\delta$"]
    name_texts = [gl.Text(line, 0.5, name, font_size=8) for line, name in zip(lines, names)]
    line_vlines = gl.Vlines(lines, colors="gray", line_styles="dashed", line_widths=1)

    data_curve = gl.Curve(wavelength_arange, data[:, *spaxel_coordinates],
                          label="data", color="black")
    model_curve = gl.Curve(wavelength_arange, total_model[:, *spaxel_coordinates],
                           label="total model", color="#ff5d00")
    continuum_curve = gl.Curve(wavelength_arange, stellar_continuum[:, *spaxel_coordinates],
                               label="stellar continuum only", color="fuchsia")
    error_curve = data_curve - model_curve
    error_curve.label = "data - model"

    fig = gl.SmartFigure(
        3,
        x_label=r"$\lambda_\mathrm{rest}$ [$\mu$m]",
        sub_y_labels=[None, r"$\nu/_\nu$ [erg s$^{-1}$ cm$^{-2}$ sr$^{-1}$]", None],
        elements=[
            name_texts,
            [data_curve, model_curve, continuum_curve, line_vlines],
            [error_curve, line_vlines],
        ],
        x_lim=(wavelength_arange.min(), wavelength_arange.max()),
        size=(14.4, 9),
        reference_labels=False,
        share_x=True,
        height_ratios=(0.2, 3, 1),
        height_padding=0,
        general_legend=True,
        legend_loc=(0.85, 0.8),
    ).set_visual_params(use_latex=True).set_ticks(minor_x_tick_spacing=0.05)
    return fig
