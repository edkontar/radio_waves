## Scattered Image Plot

Python function to plot 2D Gaussian representing an apparent radio source due to anisotropic radio-wave scattering.
Source size and position taken from rays_swind.pro simulation results, and interpolated at input frequency.

**Example:**  
```
from scattering_image_funcs import anis_scat_image  
anis_scat_image(30, 1.0, 0.25, 3000, 30)  
```  

Set the desired parameters above to generate a simulated 2D Gaussian representing a scattered image
based on ansiotropic radio-wave propagation simulation results. The resulting contours show the
50%, 70%, and 90% contour levels. The black cross denotes the initial source location in the sky-plane
and the black cross shows the scattered image centroid.  
  
The parameters are defined as:  
f: Observed frequency in MHz.  
eps0: Scattering strength factor (possible values of 0.5, 0.71, 1.0, 1.41 ,2.0).  
anis0: Anisotropy factor of the density fluctuation spectrum (possible values of 0.19, 0.25, 0.33, 0.42) where 1.0 is isotropic.  
Max plot range: Adjust the maximum x, y range.  
Source heliocentric angle: Change the initial source azimuthal angle from -60 to 60 degrees (0 degrees is the disk centre). 
