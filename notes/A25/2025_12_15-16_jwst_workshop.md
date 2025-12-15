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
## VSFs
- I used a normalization by the variance, but this doesn't affect the slope of the VSF. It just changes the y-axis scaling.
## Questions
- Loki only gives us the total SNR fits, not per component.
- Why do we take $z=0.0099$ when the proposal uses $z=0.0104$?
  - It seems that its not Joseph who gave the value, but Michael obtained it from somewhere.
## Ideas
- Create 3D gaussians and superimpose the stellar dispersion vs the gas dispersion (basically plot the heatmaps as 3D plots instead and add some smoothing, this aims at showing visually how the dispersion varies).
## TODO
- Send email to Michael to ask about the wrong chip gap shown in the fits for global optimization and a single spaxel. Since the cube is basically cropped, it seems that the "flak" file is not indexed correctly which leads to the wrong chip gap data.
- Filament vs non-filament regions for the VSFs

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
