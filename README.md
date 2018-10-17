# Xray-Fluorescence-Simulation
Simulates x-ray fluorescence for two capillary tubes containing Yttrium-based nanoparticles in a gel phantom. Written for Python 2.7.

The x-ray source is pencil beam with monoenergetic photons of 17.4 keV. X-ray fluorescence is detected with an x-ray CCD using a slit apeture. The slit is orientated perpendicular to pencil beam. The simulation produces the fluorescence image with simulated detector noise.

Run XF_PB_simulation.py to generate the x-ray fluorescence data. Then run XF_recon.py to perform the pencil beam reconstruction of the imaging slice with the XF data. Finally, run binarymask.py to produce a binary mask of the X-ray fluorescence reconstructions for later use in the x-ray luminescence deconvolution.
