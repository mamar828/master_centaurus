# Michael Reefe meeting (LOKI) to discuss the second fitting results
- Thank you for your time!
- Are you comfortable if we record this meeting since Jorge and Julie aren't able to join us today?
- We want to discuss a few ideas and questions we have after looking at the results from the second fitting run you did for us.
- I'll let Olivia start with the first couple of questions.

Questions
---
<!-- - Is there an easy way with LOKI to extract the complete fitted model for a single spaxel? For example, let's say we want to recreate the spaxel_plots figures for a specific spaxel ourselves and put it in our own figure, what do you recommend for doing that?
- Do you have any insights on how we could get a stellar mass from the fits with multiple components? We want to use this as a measure of "flux".
- The FWHM of many spaxels hits the upper limit of 800 km/s. Physically, we would like to be able to fit broader lines than this limit. Is there also a way to make sure the fit converges to a FWHM lower than the upper limit? We can see this behaviour in many FWHM param maps (e.g. H210_O2 from the "full" recent run). **Prepare figure to show this.**
--- -->
- We notice some absorption features in the spectrum aren’t being fit currently (e.g. near the Paschen alpha line) whereas others are being fit really well. We were wondering if perhaps trying to fit the spectrum only up until the chip gap could be an experiment to see if the fit improves for features like those?
  - We use model atmospheres and not real atmospheres. We could just try to add a gaussian ($e^{-\tau}$). We could add a complete polynomial over the whole spectrum.
  - There are some catalogues of real stellar spectra that could be used instead of model atmospheres.
  - Chip gap: if the calibration before and after the gap is different, it could cause problems to the model. (shouldn't be an issue as it is not mentioned in any paper)
  - We could mask the bad absorption features as a last resort.
- For the multicomponent fits, we had a nomenclature question as to how the “1” and “2” components are labeled. Is it always that “1” is the blueshifted and the “2” is the redshifted component?
  - 1 component is the brightest, and 2 is the second brightest. (sorted by flux)
  - Sort by velocity width or by velocity centroid is possible.
- Because the Br beta and O(2) lines are practically on top of each other, are there ways for us to constrain them? We were thinking that tying the O(2), O(3) and O(4) lines together and then tie all the Bracket lines as well would help break the degeneracy perhaps.
  - Seems like a good idea.
- We noticed some of the spaxels especially near the center are struggling to be fit, for example in the S(3) only run if we look at the result for spaxel (29, 32). Are the multiple components being evaluated for each spaxel individually, and if so are there constraints on how wide they can be that might be limiting the fit for spaxels like this one?
  - It gets stuck in a local minimum.
  - Local optimization is Lev-Mar.
  - When we do a global optimization (simulated annealing), the bad single components fitted where there should be two should be fixed.
- The FWHM of many lines hits an upper limit of 800 km/s. Physically, we would like to be able to fit broader lines than this limit. Is there also a way to make sure the fit converges to a FWHM lower than this upper limit? We can see this behaviour in many FWHM param maps (e.g. H210_O2 from the "full" recent run).
  - This appears for many lines in the full run, e.g. O2, Br$\beta$, Br$\gamma$.
  - It might be because no line has been detected
  - We could decompose it into a broad component and into a narrow component.
  - We could also add a stricter condition on the detection of a line.
  - Let's see if we really want to include the Brackett series in the future.
- Is there an easy way with LOKI to extract the complete fitted model for a single spaxel? For example, let's say we want to recreate the spaxel_plots figures for a specific spaxel ourselves and put it in our own figure, what would you recommend for achieving this?
  - full_model.fits, extension data
  - Each extension is a cube.
- We might want to do some additional runs (always full range so the parameters work for the whole spectrum) (tying = velocity and FWHM):
  - Now that you have a better way of modeling the stellar continuum, we might want to retry fitting with all the H lines tied
  - One where H$_2$ S, O, and Q lines are all tied separately, along with all Brackett lines
  - One where nothing is tied
  - Global optimization turned on only for a few spaxels (pick them and send them to Michael)
<!-- - We might want to do also a run where all the H$_2$ S lines are tied together (velocity and FWHM) to get average gas kinematics? From what I've seen, the S lines have very similar kinematics. -->
- How hard would it be if we wanted to do runs on our own, testing with tying different lines together, changing limits, etc? I don't want to keep bothering you.
---
Notes
- The $1.8\mu m$ absorption feature is in the atmosphere.
- Just send the spectrum. Also provide the size of the aperture because LOKI calculates in terms of intensity
- [stellar library](https://ui.adsabs.harvard.edu/abs/2009ApJS..185..289R/abstract)

Other
---
- Pa alpha and S(3) with two components in the paper, one component for the others.

## Going forward
- To improve the fits:
  - For the absorption features, I could talk to Patrick to see if he has any ideas/resources. If I don't find anything, we could just add a gaussian absorption component and in the worst case simply mask the lines.
  - For the continuum, we could add a polynomial component to the fit.
  - Mask for now, check with atmosphere folks in the meantime.
- Order the components by velocity as this will be the only was I will compare them. It is also easy to change the sorting later. If we do broad/narrow components however, we can sort by FWHM.
- Future runs (keep the double component for Pa$\alpha$ and S(3)):
  - One with all the H$_2$ lines tied.
  - One where all O, Q, and Brackett lines are tied in groups.
  - Pick some spaxels to do global optimization on. For these same spaxels, it would be interesting to try
    - all H$_2$ lines untied
    - all H$_2$ lines tied and Brackett also
    - all O, Q, and Brackett lines tied in groups
  - Did we want to test something for the chip gap? Maybe using the same spaxels than for global optimization. We could do a run where O and Q are tied (not Brackett since there are lines on both sides of the gap).
- After these runs, we will decide if we want to still fit the Brackett lines or not (there's no strong detections).
- Eventually, we will send spectra of each region (and the number of pixels). We might want to wait until we are satisfied with our fitting procedure.
---
- I will make a code to be able to zoom on the fitted spectrum using the full_model.fits file. This will allow us to investigate the FWHM problem (for now it's very hard to see if the fits with broad FWHM are fitting a line or the continuum).
<!-- - spaxel 31_27 -->
<!-- - spaxel 32_26 -->
- Global optimization only for: spaxels 32_27 and 32_37

## Simulated annealing
- Probabilistic technique to approximate the global optimum of a given function.
- For multiple local minimums, this technique is able to find the global minimum.