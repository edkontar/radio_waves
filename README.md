# radio_waves
radio wave transport in solar wind including anisotropic scattering 

; Ray tracing and scattering based on the following paper:
; http://adsabs.harvard.edu/abs/2019arXiv190900340K
; the routine works in Sun centered cartesian coordinates
; z is along the line of sight
; 
; EXAMPLE: Calling sequence, e.g. for LOFAR observations:
; IDL> ray_new, eps=0.1, anis=0.1, r_init=1.75, asym=1., f_ratio=1.1
; 
; HISTORY:
; Written: April 2018 by eduard@glasgow
; changed output file name to include assymetry and f_ratio July 11, 2019 by eduard@glasgow
; Modified on November 20, 2018 18:19UT
; The code is updated Dec 7, 2018 
; to include 1) new imaging & 2) only isotropic scattering 
; Updated April, 2019 to include anisotropic data 
; Last modified to improve output: September 4, 2019 by Eduard@Glasgow
; changed arrays rx,ry,rz to double for sintheta calculations near angle=0, September 6, 2019 by Eduard@Glasgow
