from graphinglib import SmartFigure, Heatmap
from src.tools.miscellaneous import get_pdf_image_as_array


def get_loki_line_pdfs_figure(
    lines: list[str],
    folder_name: str = "output_NGC4696_G235H_F170LP_full_model",
) -> SmartFigure:
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
    SmartFigure
        The figure containing the Loki PDF plots.
    """
    hms = []
    subtitles = []
    for line in lines:
        for suffix in ["1.flux", "1.voff", "1.fwhm"]:
            hms.append(Heatmap(
                get_pdf_image_as_array(f"data/loki/{folder_name}/param_maps/lines/{line}/{line}.{suffix}.pdf"),
                show_color_bar=False,
            ))
            subtitles.append(f"{line} {suffix.lstrip("1.")}")

    num_rows = len(hms) // 3
    fig = SmartFigure(
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
