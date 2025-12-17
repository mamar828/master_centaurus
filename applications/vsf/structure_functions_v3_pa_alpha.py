import numpy as np
import graphinglib as gl
import pyregion

from src.tools.statistics.advanced_stats import structure_function, get_fitted_structure_function_figure
from src.hdu.grouped_maps import GroupedMaps
from src.hdu.map import Map


# gm = GroupedMaps.load_from_loki("data/loki/output_NGC4696_G235H_F170LP_full_OQBr_tied/NGC4696_G235H_F170LP_full_OQBr_tied_parameter_maps.fits")
# gm["LINES.HI_PA_ALPHA.DELTA_V"].save("data/loki/vsf/v3/HI_Pa_alpha.DELTA_V.fits")
# gm["LINES.HI_PA_ALPHA.TOTAL_SNR"].save("data/loki/vsf/v3/HI_Pa_alpha.TOTAL_SNR.fits")

velocity = Map.load("data/loki/vsf/v3/HI_Pa_alpha.DELTA_V.fits")
velocity.data[velocity.data < -1e3] = np.nan
snr = Map.load("data/loki/vsf/v3/HI_Pa_alpha.TOTAL_SNR.fits")
velocity_cut = velocity.mask(snr.data > 3)
print(f"Number of remaining pixels after SNR cut: {np.sum(~np.isnan(velocity_cut.data))}")
# gl.SmartFigure(1, 2, size=(12, 5), elements=[velocity.data.plot, velocity_cut.data.plot]).show()
masked = velocity_cut.get_masked_region(pyregion.open("filament.reg"))
velocity_cut = masked
# velocity_cut.data[~np.isnan(masked.data)] = np.nan
print(f"Number of remaining pixels after mask: {np.sum(~np.isnan(velocity_cut.data))}")
# gl.SmartFigure(1, 2, size=(12, 5), elements=[velocity.data.plot, velocity_cut.data.plot]).show()

str_func = structure_function(velocity_cut.data, order=1)

# Pixels are 0.1 arcsec wide, and we have 214pc/arcsec
px_to_pc = 21.4  # pc
str_func[:, 0] *= px_to_pc
fig = get_fitted_structure_function_figure(str_func, (20, 150)).copy_with(
    x_lim=(0.9*px_to_pc, 600),
    y_lim=(150, 350),
    x_label=r"$\ell$ [pc]",
    y_label=r"$\left<\left|\delta v\right|\right>$ [km s$^{-1}$]",
)
fig[0][0].label = r"VSF for Pa$\alpha$"
fig.show()
fig.save("figures/tests/structure_functions/Pa_alpha.pdf")

# ganguly = gl.Scatter(*(np.loadtxt("data/ganguly_2023_inner.csv", delimiter=",", skiprows=1) * np.array([1000, 1])).T,
#                      face_color="red", marker_size=10, marker_style="s", label="Ganguly et al. (2023) data")
# fig.add_elements(ganguly)
# fig.show()
# fig.x_lim = 0.9*px_to_pc, 600
# fig.save("figures/tests/structure_functions/new_data/S1_filament_other_fit_bounds_no_ganguly.pdf")
