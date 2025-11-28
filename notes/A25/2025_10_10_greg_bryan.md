# Greg Bryan Meeting (invited speaker)
First talk:
- How do you use machine learning for your predictions of galaxy formation? (U-Net?)
	- Predicts the displacement between the the original and final particles' position. It works better for large scale than for small scales (chaos) and they used the power spectrum as a score metric.
- What do you use to make your simulations/simulation gifs?
- For accretion simulations, how do you deal with feedback/cooling/cooling flows/angular momentum?
	- "AGN HEATING IN SIMULATED COOL-CORE CLUSTERS"
- I didn't really understand physically why the cooling time is shortest at an intermediate temperature.
- Did you ever use structure functions to quantify turbulence/zurflueh filters?
---
Second talk:
- Solving equations of fluid dynamics (continuity equation, poisson's equation, equation of state, ...)
- Hot and cool gas cool slowly and cold, intermediate temperatures cool the most rapidly.
- Examples in astrophysics:
	- Ram pressure: the headwind strips away the gas
	- Accretion of gas onto galaxies:
	- Precipitation of cold gas out of hot gas
	- Galactic winds 
- What controls the growth or destruction of cold clouds in a hot wind?
	- Cold clouds can be rapidly destroyed (mixed): there's a mixing time determined by the size of the cloud and the background flow.
	- Radiative cooling is crucially important since it determines the time to cool. At an intermediate temperature, the cooling time is very short.
	- Very turbulent radiative mixing layer (interface between cold and hot gas) -> it seems to drive turbulence.
	- Theoretically, its the pressure that drives hot gas to cold gas through the interface. With high resolution simulations, they find that the hot/cold interface is fractal, since the area depends on the averaging length. 
	- The resolving theory is that the mixing is not like a waterfall but more like a water wheel.
	- The gas is mixed at the same rate that it cools.
	- The constant inflow of hot gas nurtures this interface.
- Let's test the numerical problem with the galactic wind problem: SN drive out gas, small clouds are shredded, large clouds gros and accelerated out into the CGM -> multiphase winds regulate galaxies cosmic evolution. A single equation can be constructed to model the evolution of the clouds.
*Check cooling flows/cooling flow problem.*
---
AFTER:
- It could be interesting to look at the VSF/power spectrum/fourier transform of the ACR.
- "AGN HEATING IN SIMULATED COOL-CORE CLUSTERS".
- The UNET is used to predict the displacement of particles between the start and the end of the simulation. When looking at large structures, the predictions match the simulations but for small structures, the chaotic nature makes the correspondence less clear.
- They simulated the cooling flow problem and also found that AGN feedback is probably responsible for heating the medium.
