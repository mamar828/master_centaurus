from __future__ import annotations
import numpy as np
from astropy.io import fits
from difflib import get_close_matches

from src.hdu.map import Map
from src.hdu.arrays.array_2d import Array2D
from src.hdu.header import Header
from src.tools.miscellaneous import silence_function


class GroupedMaps:
    """
    This class implements a container for maps that are linked together. This is useful for opening FITS files that
    contain multiple maps as different extensions and converting them to `Map` objects. This class allows to index the
    different maps using either attributes or keys.
    """

    def __init__(self, maps: list[tuple[str, Map]]) -> None:
        """
        Initializes a GroupedMaps object.

        Parameters
        ----------
        maps : list[tuple[str, Map]]
            List of (key, value) pairs where the key is the name of the map and the value is the `Map` object itself.
        """
        self.map_dict = {name: map_ for name, map_ in maps}
        self.names = list(self.map_dict.keys())
        self.maps = list(self.map_dict.values())

    def __len__(self) -> int:
        return len(self.names)

    def __getitem__(self, name: str) -> Map:
        if name in self.names:
            return self.map_dict[name]
        else:
            close_match = get_close_matches(name, self.names, n=1, cutoff=0.6)
            if close_match:
                raise AttributeError(f"GroupedMaps has no attribute '{name}'. Did you mean '{close_match[0]}'?")
            else:
                raise AttributeError(f"GroupedMaps has no attribute '{name}'.")

    def __getattr__(self, name: str) -> Map:
        return self[name]

    def __str__(self) -> str:
        return f"GroupedMaps object with {len(self)} maps."

    @classmethod
    def load(cls, filename: str, mute: bool = False) -> GroupedMaps:
        """
        Loads a collection of maps from a FITS file to a GroupedMaps object

        Parameters
        ----------
        filename : str
            Name of the file to load.
        mute : bool, default=False
            Whether to completely mute the loading process. This is useful when opening FITS files creates warnings.

        Returns
        -------
        GroupedMaps
            A GroupedMaps object containing the maps stored in the file.
        """
        def load_func():
            hdu_list = fits.open(filename)
            maps = [(hdu.name, Map(data=Array2D(hdu.data), header=Header(hdu.header))) for hdu in hdu_list[1:]]
            return maps

        if mute:
            maps = silence_function(load_func)()
        else:
            maps = load_func()

        return cls(maps)
