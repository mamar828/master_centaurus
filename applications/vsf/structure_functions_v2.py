import numpy as np
import graphinglib as gl

from src.tools.statistics.advanced_stats import structure_function, get_fitted_structure_function_figure
from src.hdu.grouped_maps import GroupedMaps
from src.hdu.map import Map


# gm = GroupedMaps.load_from_loki("data/loki/output_NGC4696_G235H_F170LP_full_model_stars_full/"
#                                 "NGC4696_G235H_F170LP_full_model_stars_full_parameter_maps.fits")
# gm["LINES.HI_PA_ALPHA.1.VOFF"].save("data/loki/HI_PA_ALPHA_voff_map.fits")
# gm["LINES.HI_PA_ALPHA.TOTAL_SNR"].save("data/loki/HI_Pa_alpha.TOTAL_SNR.fits")
# gm["LINES.H210_S3.1.VOFF"].save("data/loki/H210_S3_voff_map.fits")
# gm["LINES.H210_S3.TOTAL_SNR"].save("data/loki/H210_S3.TOTAL_SNR.fits")
# gm["LINES.H210_O3.1.VOFF"].save("data/loki/H210_O3_voff_map.fits")
# gm["LINES.H210_O3.TOTAL_SNR"].save("data/loki/H210_O3.TOTAL_SNR.fits")

# snr = Map.load("data/loki/HI_Pa_alpha.TOTAL_SNR.fits")
# velocity = Map.load("data/loki/HI_PA_ALPHA_voff_map.fits")
# snr = Map.load("data/loki/H210_S3.TOTAL_SNR.fits")
# velocity = Map.load("data/loki/H210_S3_voff_map.fits")
snr = Map.load("data/loki/H210_O3.TOTAL_SNR.fits")
velocity = Map.load("data/loki/H210_O3_voff_map.fits")
velocity_cut = velocity.mask(snr.data > 3)  # 796 pixels remain
# gl.SmartFigure(1, 2, size=(12, 5), elements=[velocity.data.plot, velocity_cut.data.plot]).show()

str_func = structure_function(velocity_cut.data, order=1)
# gl.SmartFigure(elements=[gl.Scatter(str_func[:, 0], str_func[:, 1])], log_scale_x=True, log_scale_y=True).show()
# fig = get_fitted_structure_function_figure(str_func, (np.sqrt(5), 10))
# fig.x_lim = 0.9, 20
# fig.y_lim = 1.7e-3, 6e-3
# fig.show()
# fig.save("figures/tests/structure_functions/power.pdf")

# TRY BINNING the str_func in logspace
# bins = np.logspace(np.log10(str_func[0, 0]), np.log10(str_func[-1, 0]), 30)
# centers =
# digitized = np.digitize(str_func[:, 0], bins)
# binned_means = [str_func[digitized == i].mean(axis=0) for i in range(1, len(bins))]
# print(binned_means)
# raise

# Pixels are 0.1 arcsec wide, and we have 214pc/arcsec
str_func_pc = str_func.copy()
str_func_pc[:, 0] *= 21.4
fig = get_fitted_structure_function_figure(str_func_pc, (np.sqrt(5)*21.4, 7*21.4)).copy_with(
    x_lim=(0.9*21.4, 20*21.4),
    y_lim=(10, 300),
    x_label=r"$\ell$ [pc]",
    y_label=r"$\left<\left|\delta v\right|\right>$ [km s$^{-1}$]",
)
fig[0][0].label = "VSF for O(3)"
# fig.add_elements(gl.Vlines(bins*21.4))
fig.show()
# fig.save("figures/tests/structure_functions/S3/parsecs.pdf")

# ganguly = gl.Scatter(*(np.loadtxt("data/ganguly_2023_inner.csv", delimiter=",", skiprows=1) * np.array([1000, 1])).T,
#                      face_color="red", marker_size=10, marker_style="s")
# fig.add_elements(ganguly)
# fig.x_lim = 0.9*21.4, 4000
# fig.show()
# fig.save("figures/tests/structure_functions/O3_bins.pdf")
