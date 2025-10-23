import numpy as np
import os
import warnings
from contextlib import redirect_stdout
from pdf2image import convert_from_path


def silence_function(func):
    """
    Decorates verbose functions to silence their terminal output.
    """
    def inner_func(*args, **kwargs):
        with open(os.devnull, "w") as outer_space, redirect_stdout(outer_space), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return func(*args, **kwargs)

    return inner_func

def get_pdf_image_as_array(filename: str, page_number: int = 0, dpi: int = 300) -> np.ndarray:
    """
    Gives a numpy array representation of a specific page from a PDF file.

    Parameters
    ----------
    filename : str
        Path to the PDF file.
    page_number : int, default=0
        Page number to convert (0-indexed).
    dpi : int, default=300
        Resolution for the conversion, in dots per inch.

    Returns
    -------
    np.ndarray
        Numpy array representation of the specified PDF page.
    """
    pages = convert_from_path(filename, dpi=dpi)
    image = np.array(pages[page_number])
    return image
