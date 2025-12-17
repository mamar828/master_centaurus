# JWST Workshop
## Étienne Artigau
### Pre-discussion
- The stellar population we mainly have are M dwarfs. We can see this from the depth of the CO band heads $\implies$ deeper CO band heads $=$ cooler stars.
- He sent me a template for a NIRISS M dwarf spectrum (LHS 1140) that I could use to fit the continuum better.
  - The absorption feature at 1.8 microns seems to be in that spectrum.
  - The feature at 1.845 microns is in a region of the spectrum where there is nothing. We didn't see it in another spectrum he showed me.
- We could get the temperature of the stars from the depth of the CO band heads.
### During the workshop
- Make a plot of the spectral type vs the depth of the CO band head.
- We could do an average spectrum around let's say 10 pc using a weighted average of the different stars.
- Plot: IRTF convolved to NIRSpec resolution, depth vs spectral type. We can also normalize to a given flux.
  - https://irtfweb.ifa.hawaii.edu/~spex/IRTF_Spectral_Library/References_files/K.html
- Supergiants are very less common, but way more bright, so it cancels out ish. Its hard to know if flux weighted, which star type dominates.
- Maybe for a summer internship:
  - Naines IRTF (type V) G to ..., trouver température de l'étoile (SIMBAD), plot du type de l'étoile vs la
  - Type spectral pondéré par le flux en croppant certaines wavelengths du spectre
## Questions
- Loki only gives us the total SNR fits, not per component.
- Why do we take $z=0.0099$ when the proposal uses $z=0.0104$?
  - It seems that its not Joseph who gave the value, but Michael obtained it from somewhere.
- Difference between $v_\mathrm{off}$ and $\Delta v$?
- Weirdly, some S lines have offset velocities really different than the other S lines (skewed).
## Ideas
- Create 3D gaussians and superimpose the stellar dispersion vs the gas dispersion (basically plot the heatmaps as 3D plots instead and add some smoothing, this aims at showing visually how the dispersion varies).
## Conference
- At registering, give choices (4-5)
  - Socks
  - Deck of cards
  - Water bottles
  - Mugs
  - Pen and paper
  - Maybe one UdeM thing
  - Québec tea towels
  - Montreal map (online (google maps list of places) + physical)
- Look at this with Anabel
## VSFs
- Intuition:
  - Single power law on small scales $\implies$ expected from a turbulent flow. Multiple slopes/not linear $\implies$ multiple energy ejections at these scales.
  - Shallow slope: velocities at small and large scales are similar, small scale motions are as strong as large scale motions, motions are dominated by local processes like shocks.
  - Steep slope: not reproduced by simulations with magnetic fields and a little far from kinetic plasma simulations, it most likely comes from SMBH feedback that does not drive "volume-filling turbulence" efficiently
## Meeting with Gary Ferland
- It's hard to get H$_2$ really hot without destroying it.
- The ionization parameter is the ration of star light to gas density
- Next meeting Thursday 8 January 2 PM
## Meeting in January with the whole team
- Joseph summary of data reduction (including a slide on LOKI)
- Send a When2meet for the week of January 19
  - We want to summarize the first results
  - Joseph will send a mailing list
  - Send it to Julie first
### Letter
- Present the data
- Connecting feedback between kpc to sub-kpc scales (completing the feedback loop)
- Filament connecting to the rotating disk
### Kinematics
#### Paper 1 (Mathieu)
- General velocity flow (velocity + vel. dispersion of all detected lines)
- Mention lines that are not detected
- Single components for all the lines (except pa alpha and S3)
- Kinematics of the stars
- P-V diagrams for gas and stars
- Compare with MUSE, ALMA, XRISM
- 3D plot of the velocity dispersion of stars vs gas (overlayed)
- Approaching/receding component for Pa alpha and S(3)
#### Paper 2 (Mathieu)
- VSFs for all the gas vs the filaments vs non filament, add MUSE VSF
- VSFs for the stars (fit two lines and look at the breaking scale)
- Try with multiple components (VSFs narrow vs broad)
- Plots of slope as a function of temperature
- See 1st order vs 2nd order?
- PDFs for relevant gas, also do it for MUSE
#### Paper 3 (another student)
- What are the stars in there
- Maybe codirection with Étienne Artigau?
### Lines
#### Paper 1 (Olivia)
- All lines that are detected and not detected (table that gives the details of each line in each region)
- BPT type diagram with Webb vs MUSE on two axes
- BPT S something vs S(1)
- HST data, most important lines we see in the data
- Excitation diagram
- Some tests with Tremblay's paper
## Email for Michael
- Wrong chip gap taken for the single spaxels
- Chip gap is too large for individual regions also
- Ask where he got the redshift from.
- Ask if he masks the line to fit the continuum or if the continuum and emission lines are fit simultaneously.
- Tie the bracket lines with the narrow component of Pa alpha
- Why is there a spike in the 32_27 spaxel (also in the initial_sum_fit spectrum)
- Why are there grey bands at the low spaxels numbers
- Ask for the high resolution global optimization fits
- What is the difference between v_off and delta_v?
- Also give the smoothed data to Michael (wait for Joseph to reduce the data to the best resolution)
- Send to Michael the region spectra with the new resolution
---
- No global H2, no chip gap specific fits.
## TODO
- [ ] Send email to Michael to ask about the wrong chip gap shown in the fits for global optimization and a single spaxel. Since the cube is basically cropped, it seems that the "flak" file is not indexed correctly which leads to the wrong chip gap data.
- [x] Filament vs non-filament regions for the VSFs
- [x] Implement the Loki fit figure for all versions
- [ ] Convert arcsecs to pc with [NED](https://www.astro.ucla.edu/~wright/CosmoCalc.html)
### For January meeting
- [ ] Approaching/receding component figure for Pa alpha and S(3)
