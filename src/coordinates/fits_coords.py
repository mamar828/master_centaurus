from __future__ import annotations


class FitsCoords:
    """
    Encapsulates the methods specific to .fits coordinates, namely 1-based indexing and inverted axes order.
    An object of this class may be used just like a tuple of ints, but the values are automatically converted to
    numpy-compatible coordinates when accessed. You can also unpack this object to get the converted values as ints.

    Examples
    --------
    spectrum = Cube[:,*FitsCoords(5,10)]
    will correctly slice the Cube at the specified coordinates, which returns a Spectrum.

    The following lines show how useful this class is for precisely slicing a FitsObject. To get the value at the
    (x, y, z) = (8, 17, 225) on a .fits file, the following code is equivalent :
    using FitsCoords :   print(Cube[*FitsCoords(8, 17, 225)])
    without FitsCoords : print(Cube[224, 16, 7])
    In the former case, the coordinates can intuitively be entered and in the latter, an awkward double conversion is
    required.
    """

    def __init__(self, *coordinates: tuple[int]):
        """
        Initialize a FitsCoords object.

        Parameters
        ----------
        coordinates : tuple[int]
            Coordinates to initialize the object with. These are given in the same order as in fits files and follow
            one-based indexing.
        """
        self.data = list(coordinates)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, key: int) -> int:
        """
        Gives the value of the specified key by converting the fits coordinate to a numpy one. The conversion is made by
        inverting the coordinates order (e.g. x,y -> y,x) as numpy indexing starts with the "last index", then by
        subtracting 1 because fits indexing starts at 1 and not 0.
        """
        # Roll axis to account for the coordinate switch
        index = len(self) - 1 - key
        if index >= 0:
            # Reduce axis to account for indexing differences
            coord = self.data[index] - 1
            return coord
        else:
            # This clause is required by the unpacking operator
            raise IndexError

    def __str__(self):
        return f"FitsCoords({', '.join(map(str, self.data))})"

    @classmethod
    def from_python(cls, *coordinates: tuple[int]) -> FitsCoords:
        """
        Create a FitsCoords object from Python coordinates.

        Parameters
        ----------
        coordinates : tuple[int]
            Coordinates to initialize the object with. These are given in the same order as in Python.

        Returns
        -------
        FitsCoords
            A new FitsCoords object initialized with the provided coordinates.
        """
        coords = coordinates[::-1]
        coords = [c + 1 for c in coords]
        return cls(*coords)
