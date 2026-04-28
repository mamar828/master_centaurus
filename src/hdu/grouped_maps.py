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

    def __init__(self, maps: list[tuple[str, Map]], header: Header | None = None) -> None:
        """
        Initializes a GroupedMaps object.

        Parameters
        ----------
        maps : list[tuple[str, Map]]
            List of (key, value) pairs where the key is the name of the map and the value is the `Map` object itself.
        header : Header, optional
            A header that can be associated with the whole collection of maps. This is typically the first header of a
            hdu list, which isn't associated with any map but contains metadata about the collection of maps.
        """
        self.map_dict = {name: map_ for name, map_ in maps}
        self.names = list(self.map_dict.keys())
        self.maps = list(self.map_dict.values())
        self.header = header

    def __len__(self) -> int:
        return len(self.names)

    def __getitem__(self, key: str | int) -> Map:
        if key in self.names:
            return self.map_dict[key]
        elif isinstance(key, int):
            if 0 <= key < len(self.maps):
                return self.maps[key]
            else:
                raise IndexError(f"Index {key} is out of range for GroupedMaps with {len(self)} maps.")
        else:
            close_match = get_close_matches(key, self.names, n=1, cutoff=0.6)
            if close_match:
                raise AttributeError(f"GroupedMaps has no attribute '{key}'. Did you mean '{close_match[0]}'?")
            else:
                raise AttributeError(f"GroupedMaps has no attribute '{key}'.")

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
            header = Header(hdu_list[0].header)
            maps = [(hdu.name, Map(data=Array2D(hdu.data), header=Header(hdu.header))) for hdu in hdu_list[1:]]
            return maps, header

        if mute:
            maps, header = silence_function(load_func)()
        else:
            maps, header = load_func()

        return cls(maps, header)
