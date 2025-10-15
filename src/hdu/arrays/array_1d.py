from __future__ import annotations
from graphinglib import Curve
import numpy as np

from src.hdu.arrays.array import Array


class Array1D(Array):
    """
    Encapsulates the methods specific to one-dimensional arrays.
    """

    @property
    def plot(self) -> Curve:
        """
        Gives the plot of the Array1D with a Curve.

        Returns
        -------
        Curve
            Plotted Array1D.
        """
        curve = Curve(
            x_data=np.arange(0, len(self)),
            y_data=self,
        )
        return curve

    def get_nan_cropping_slices(self) -> slice:
        """
        Gives the slices for cropping the border nan values of the array.

        Returns
        -------
        slice
            Slice that must be applied on the Array1D for cropping.
        """
        non_nan_vals = ~np.isnan(self)
        slice_ = slice(np.argmax(non_nan_vals), non_nan_vals.shape[0] - np.argmax(non_nan_vals[::-1]))
        return slice_
