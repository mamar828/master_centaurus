import numpy as np
import graphinglib as gl

from src.tools.statistics.advanced_stats import structure_function, get_fitted_structure_function_figure
from src.hdu.grouped_maps import GroupedMaps
from src.hdu.map import Map


# gm = GroupedMaps.load_from_loki("data/loki/output_NGC4696_G235H_F170LP_full_model_stars_full/"
#                                 "NGC4696_G235H_F170LP_full_model_stars_full_parameter_maps.fits")
# gm["LINES.HI_PA_ALPHA.1.VOFF"].save("data/loki/HI_PA_ALPHA_voff_map.fits")
# gm["LINES.HI_PA_ALPHA.TOTAL_SNR"].save("data/loki/HI_Pa_alpha.TOTAL_SNR.fits")

snr = Map.load("data/loki/HI_Pa_alpha.TOTAL_SNR.fits")
velocity = Map.load("data/loki/HI_PA_ALPHA_voff_map.fits")
velocity_cut = velocity.mask(snr.data > 3)  # 796 pixels remain
# gl.SmartFigure(1, 2, size=(12, 5), elements=[velocity.data.plot, velocity_cut.data.plot]).show()

str_func = structure_function(velocity_cut.data, order=1)
# gl.SmartFigure(elements=[gl.Scatter(str_func[:,0], str_func[:,1])], log_scale_x=True, log_scale_y=True).show()
fig = get_fitted_structure_function_figure(str_func, (0.45, 1))
fig1 = fig.copy_with(x_lim=(None, None), y_lim=(-2.8, -2.2))
fig = get_fitted_structure_function_figure(str_func, (0.45, 1.4))
fig2 = fig.copy_with(x_lim=(None, None), y_lim=(-2.8, -2.2))
gl.SmartFigure(1, 2, size=(12, 5), elements=[fig1, fig2]).save("figures/tests/structure_functions/first_test.pdf")
