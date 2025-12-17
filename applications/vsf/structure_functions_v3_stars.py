import numpy as np
import graphinglib as gl
import pyregion

from src.tools.statistics.advanced_stats import structure_function, get_fitted_structure_function_figure
from src.hdu.grouped_maps import GroupedMaps
from src.hdu.map import Map


# gm = GroupedMaps.load_from_loki("data/loki/output_NGC4696_G235H_F170LP_full_OQBr_tied/NGC4696_G235H_F170LP_full_OQBr_tied_parameter_maps.fits")
# gm["CONTINUUM.STELLAR_KINEMATICS.VEL"].save("data/loki/vsf/v3/STELLAR_KINEMATICS.fits")

velocity = Map.load("data/loki/vsf/v3/STELLAR_KINEMATICS.fits")
print(f"Number of pixels: {np.sum(~np.isnan(velocity.data))}")
# gl.SmartFigure(1, 2, size=(12, 5), elements=[velocity.data.plot, velocity_cut.data.plot]).show()

str_func = structure_function(velocity.data, order=1)

# Pixels are 0.1 arcsec wide, and we have 214pc/arcsec
px_to_pc = 21.4  # pc
str_func[:, 0] *= px_to_pc
fig = get_fitted_structure_function_figure(str_func, (20, 150)).copy_with(
    # x_lim=(0.9*px_to_pc, 600),
    # y_lim=(150, 350),
    x_label=r"$\ell$ [pc]",
    y_label=r"$\left<\left|\delta v\right|\right>$ [km s$^{-1}$]",
)
fig[0][0].label = r"VSF for the stars"
fig.show()
fig.save("figures/tests/structure_functions/stars.pdf")

# ganguly = gl.Scatter(*(np.loadtxt("data/ganguly_2023_inner.csv", delimiter=",", skiprows=1) * np.array([1000, 1])).T,
#                      face_color="red", marker_size=10, marker_style="s", label="Ganguly et al. (2023) data")
# fig.add_elements(ganguly)
# fig.show()
# fig.x_lim = 0.9*px_to_pc, 600
# fig.save("figures/tests/structure_functions/new_data/S1_filament_other_fit_bounds_no_ganguly.pdf")
