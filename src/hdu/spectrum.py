from __future__ import annotations
import numpy as np
from typing import Self
from graphinglib import Curve
from copy import deepcopy

from src.hdu.fits_object import FitsObject
from src.hdu.arrays.array_1d import Array1D
from src.hdu.header import Header
from src.base_objects.silent_none import SilentNone


class Spectrum(FitsObject):
    """
    Encapsulate all the methods of any data cube's spectrum.
    """

    def __init__(self, data: Array1D, header: Header = SilentNone()):
        """
        Initializes a Spectrum object with a certain header, whose spectral information will be taken.

        Parameters
        ----------
        data : Array1D
            Detected intensity at each channel.
        header : Header, default=SilentNone()
            Spectral header of the Spectrum.
        """
        self.data = Array1D(data)
        self.header = header

    def __len__(self) -> int:
        return len(self.data)

    def copy(self) -> Spectrum:
        return deepcopy(self)

    @property
    def isnan(self) -> bool:
        return np.all(np.isnan(self.data))

    @property
    def x_values(self) -> np.ndarray:
        """
        Gives the x values associated with the Spectrum's data.

        Returns
        -------
        np.ndarray
            Range from 1 and has the same length than the data array. The start value is chosen to match with SAOImage
            ds9 and with the headers, whose axes start at 1.
        """
        return np.arange(1, len(self) + 1)

    @property
    def plot(self) -> Curve:
        """
        Gives the plot of the spectrum with a Curve.

        Returns
        -------
        Curve
            Curve representing the spectrum's values at every channel.
        """
        curve = Curve(
            x_data=self.x_values,
            y_data=self.data,
            label="Spectrum"
        )
        return curve

    def bin(self, bin: int, ignore_nans: bool = False) -> Self:
        """
        Bins a Spectrum.

        Parameters
        ----------
        bin : int
            Number of pixels to be binned together. A value of 1 results in the Spetrum not being binned.
        ignore_nans : bool, default=False
            Whether to ignore the nan values in the process of binning. If no nan values are present, this parameter is
            obsolete. If False, the function np.mean is used for binning whereas np.nanmean is used if True. If the nans
            are ignored, the map might increase in size as pixels will take the place of nans. If the nans are not
            ignored, the map might decrease in size as every new pixel that contained a nan will be made a nan also.

        Returns
        -------
        Self
            Binned Spectrum.
        """
        return self.__class__(
            self.data.bin(bin, ignore_nans),
            self.header.bin(bin)
        )
