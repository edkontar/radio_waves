import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter, NullFormatter
from scipy import integrate, interpolate
from astropy.io import fits
from astropy.coordinates import cartesian_to_spherical

class Scattering:
    """
    Ray tracing and scattering based on the following paper:
    http://adsabs.harvard.edu/abs/2019arXiv190900340K
    the routine works in Sun centered cartesian coordinates
    z is along the line of sight
    
    EXAMPLE: Calling sequence, e.g. for LOFAR 30 MHz observations:
    from rays_swind import Scattering
    sim = Scattering(qeps2_scaling = 1.0, anis = 0.3, r_init = 1.75, asym = 1.0, f_ratio = 1.1)
    sim.run()
     
    PARAMETERS:
    qeps2_scaling - scaling of qeps2 function (used 1./4,1/2, 1, 2, 4); default = 1 like in
                    ApJ Paper https://doi.org/10.3847/1538-4357/acf6c1
     
    anis   - wavenumber anisotropy of density turbulence q_par/q_perp   ; default =0.25 
             anis->+infty, quasi-paraller density fluctuations
             anis=1 isotropic density fluctuations 

    r_init - initial location of the source (to check local plasma frequency IDL> omega_pe(r_init)/2/!Pi/1e6 [MHz])
    r_init=1.75D0 ;(around 32 MHz)
    r_init=5.D0 ;(around 2 MHz)
    r_init=2.21D0 ;(around 15 MHz)
    r_init=1.1372D0 ;(around 329 MHz)
    r_init=1.3D0 ;(around 100 MHz)
    r_init=57.D0 ;(around 0.118 MHz)
    r_init=10.D0 ;(around 0.765 MHz)
    r_init=20.D0 ;(around 0.338 MHz)
       
    f_ratio - ratio of emitted frequency to local plasma frequency at the source
              i.e. omega/omega_pe(r_init)
     
    asym  - asymetry parameter density fluctuations assymetry sunward/anisunward 
    asym=1 - symmetric;
    0<asym<1 => more outward density fluctuation (inward radio stronger scattered)
    1<asym<2 => more inward density fluctuations (outward radio stronger scattered)
    <delta n^2>_in= asym *<delta n^2>
    <delta n^2>_out = (2- asym) *<delta n^2>
    z_anis=asym*(k_along_r GT 0.) + (2.-asym)*(k_along_r LE 0.)

    HISTORY:
    Python version updates:
    Translated to Python: July 2020 by DanielClarkson@Glasgow
    Updated python version from ray_new to rays_swind following IDL version: Nov 20, 2023 by DanielClarkson@Glasgow.
    IDL Version updates:
    Written: April 2018 by eduard@glasgow
    changed output file name to include assymetry and f_ratio July 11, 2019 by eduard@glasgow
    Modified on November 20, 2018 18:19UT
    The code is updated Dec 7, 2018
    to include 1) new imaging & 2) only isotropic scattering
    Updated April, 2019 to include anisotropic data
    Last modified to improve output: September 4, 2019 by Eduard@Glasgow
    changed arrays rx,ry,rz to double for sintheta calculations near angle=0, September 6, 2019 by Eduard@Glasgow
    updated: October 2023 by eduard@glasgow
    include qeps2 profile
    """

    # FUNCTIONS
    def __init__(self, qeps2_scaling, anis, r_init, asym, f_ratio):
        self.qeps2_scaling = qeps2_scaling
        self.anis = anis
        self.r_init = r_init
        self.asym = asym
        self.f_ratio = f_ratio

    def density_r(self, rr):
        # returns plasma density in cm^-3

        # density in the low corona
        h0 = 144/6.96e5 # 144 km/solar radius [km]

        h1 = 20/960 # 20 arcsec/solar radius [arcsec]

        #nc=1.17e17*exp(-(rr-1.)/(h0))+1e11*exp(-(rr-1.)/0.02)
        nc = 3e11*np.exp(-(rr-1)/h1)

        # Newkirk model
        #return, 4.2e4*10.^(4.2/rr)

        #Leblanc 1998 model
        #return, 3.3e5*rr^(-2.)+ 4.1e6*rr^(-4.)+8e7*rr^(-6)

        # Model that fits Parker model but simple power-laws
        return 4.8e9/rr**14 + 3e8/rr**6 + 1.39e6/rr**2.3 #+ nc

    def omega_pe(self, r):
        # returns plasma frequency
        return 8.98e3 * np.sqrt(self.density_r(r))*2*np.pi

    def nu_scat_fixed(self, r, omega, qeps2_scaling):
        # scattering power
        c = 2.998e10
        w_pe = self.omega_pe(r)
        q_eps2 = qeps2_scaling * 2e3 / r**0.7 * ((r - 1)/r)**2.7 / 6.96e10
        nu_s = np.pi/8 * q_eps2 * w_pe**4 * c/omega/(omega**2 - w_pe**2)**1.5
        return nu_s

    def nu_scat(self, r, omega, qeps2_scaling):
        c = 2.998e10
        
        nu_s = self.nu_scat_fixed(r, omega, qeps2_scaling)
        # changes here to avoid other changes in the code
        return nu_s # units of rad^2/s

    def domega_dr(self, r):
        # returns domega_pe/dr
        # INPUT: r in distance of R_sun units
        r_tab = 10**(np.linspace(0, 2.99, num = 300))
        d_tab1 = np.gradient(r_tab)
        d_tab2 = np.gradient(self.omega_pe(r_tab))
        d_tab = d_tab2/d_tab1 # gives similar resuls to IDL deriv(x,y)
        interpol = interpolate.interp1d(r_tab, d_tab, kind = 'linear')
        return interpol(r) # this and above line is representation of IDL interpol

    def parker_Bxyz_Earth(self, r, theta, fi):
        # input r, theta, fi so that z is north x to the Earth (HEEQ coordinates)
        # output: BzE, ByE, BxE

        B0 = 5e-9*1e4 # Gauss at 1AU
        r0 = 215
        vsw = 4.2e7 # solar wind speed 1AU
        #vsw = 1e9 # 10000 km/s to check radial
        Omega = 2*np.pi/(27*24*3600) # solar rotation rate

        # parker spiral in spherical coordinates
        Br     = B0*(r0/r)**2
        Btheta = 0
        Bfi    =-B0*(r0/r)*(Omega*r0*6.95e10/vsw)*np.sin(theta)

        # converting to cartesian coordinates
        BxE = Br*np.sin(theta)*np.cos(fi) - Bfi*np.sin(fi)
        ByE = Br*np.sin(theta)*np.sin(fi) + Bfi*np.cos(fi)
        BzE = Br*np.cos(theta)

        return BxE, ByE, BzE

    def run(self):

        qeps2_scaling = self.qeps2_scaling
        anis = self.anis
        r_init = self.r_init
        asym = self.asym
        f_ratio = self.f_ratio
        
        # Constants
        N = 999 # photon number
        R_s = 6.96e10 # solar radius [cm]
        c = 2.998e10 # c [cm s^-1]
        c_r = c/R_s # c [R_s]
        c_s = 2e7/R_s # sound speed [R_s]
        au = 215 # 1 AU [R_s]

        theta_0 = 0
        r = np.linspace(0, N, num = N, dtype = float)*0 + r_init
        rtheta = (theta_0 + np.linspace(0, N, num = N, dtype = float)/N*0) * np.pi/180
        rfi = (0 + np.linspace(0, N, num = N, dtype = float)/N*0) * np.pi/180
        # r = r + (np.random.uniform(0, 1, size = N)) * 1e-1

        Te = 86.0 * 1.5 # temperature, eV
        nu_e = 2.91e-6 * self.density_r(r_init) * 20/(Te**1.5)
        rz = r * np.cos(rtheta)
        rx = r * np.sin(rtheta) * np.cos(rfi)
        ry = r * np.sin(rtheta) * np.sin(rfi)
        Iindex = np.linspace(0, N-1, num = N, dtype = int)

        # large (10 R_sun) line-like source
        #rx = rx + (np.random.uniform(-5, 5, size = N))

        # Initial locations for future use
        r0 = r
        rx0, ry0, rz0 = rx, ry, rz

        f_pe0 = self.omega_pe(r)/2/np.pi
        omega_pe0 = self.omega_pe(r_init) # initial plasma frequency
        f_start = omega_pe0/1e6/2/np.pi
        print('INIT: Source located at R_source = {0:.2f} R_sun, Angle = {1:.2f} deg'.format(r_init, theta_0))
        print('f_pe = {0:.3f} MHz for emission frequency = {1:.3f} MHz'.format(f_start, f_start*f_ratio))
        print('KEY PARAMETERS: qeps2_scaling = {0:.2f}, anisotropy = {1:.2f}, asymmetry = {2:.2f}'.format(qeps2_scaling, anis, asym))

        omega = f_ratio * self.omega_pe(r0)
        # emission at 1.1 of plasma frequency in units of initial omega_pe
        # 2 for harmonic
        omega0 = min(omega)
        omega = omega0

        kc0 = np.sqrt(omega**2 - self.omega_pe(r)**2)
        kc = kc0
        tau = np.zeros(N, dtype = float)

        np.random.seed(1001) # mu and kfi differ the IDL due to random function
        # mu = np.random.uniform(-1, 1, size = N)) # isotropic
        mu = np.random.uniform(0, 1, size = N) # positive mu, only to the observer
        kfi = np.random.uniform(0, 1, size = N)*2*np.pi
        kz = kc * mu
        kx = kc * np.sqrt(1 - mu**2) * np.cos(kfi)
        ky = kc * np.sqrt(1 - mu**2) * np.sin(kfi)

        pp = np.linspace(0, 41, num = 41)/41 * np.pi * 2
        fig, ax = plt.subplots(2, 2, figsize = (8,8))

        ax[0,0].plot(rz0, rx0, 'kx')
        ax[0,0].plot(np.cos(pp), np.sin(pp), linewidth = 2, color = 'gold')
        ax[0,0].plot([0,0], [-3*r_init, 3*r_init], '--', color = 'grey', linewidth = 0.8)
        ax[0,0].plot([-3*r_init, 3*r_init], [0,0], '--', color = 'grey', linewidth = 0.8)
        ax[0,0].set_xlabel('Z [R$_\odot$]')
        ax[0,0].set_ylabel('X [R$_\odot$]')
        ax[0,0].set_xlim(-3*r_init, 3*r_init)
        ax[0,0].set_ylim(-3*r_init, 3*r_init)
        ax[0,0].set_aspect('equal', 'box')

        ax[0,1].plot(rx0, ry0, 'kx')
        ax[0,1].plot(np.cos(pp), np.sin(pp), linewidth = 2, color = 'gold')
        ax[0,1].plot([0,0], [-3*r_init, 3*r_init], '--', color = 'grey', linewidth = 0.8)
        ax[0,1].plot([-3*r_init, 3*r_init], [0,0], '--', color = 'grey', linewidth = 0.8)
        ax[0,1].set_xlabel('X [R$_\odot$]')
        ax[0,1].set_ylabel('Y [R$_\odot$]')
        ax[0,1].set_xlim(-3*r_init, 3*r_init)
        ax[0,1].set_ylim(-3*r_init, 3*r_init)
        ax[0,1].set_aspect('equal', 'box')


        # density fluctuations
        nu_s0 = self.nu_scat(r_init, omega, qeps2_scaling) # IDL version uses min(omega), however omega is single value

        #rint = r_init * (1 + 10**(np.linspace(0, 999, num = 999)/999*5-3.5))
        rint = r_init * (1 + np.linspace(0,3999, num = 3999)/49 + 1e-3)

        N_rint = len(rint)
        dri = np.roll(rint, -1) - rint
        dri[N_rint-1] = 0
        tau_scat = np.zeros(N_rint)

        for j in range(len(rint)-3):
            tau_scat[j] = integrate.simps(
                            self.nu_scat(rint[j:N_rint-1], omega0, qeps2_scaling) / c_r / np.sqrt(1-self.omega_pe(rint[j:N_rint-1])**2 / omega0**2),
                            x = rint[j:N_rint-1]                    
                        )

        tau_scat[len(rint)-1] = tau_scat[len(rint)-4]
        tau_scat[len(rint)-2] = tau_scat[len(rint)-4]
        tau_scat[len(rint)-3] = tau_scat[len(rint)-4] # tau_scat has 3 final values of 0, whereas IDL version has 2, so extra line here required

        tau_scat10 = abs(tau_scat - 1/10) # extra var not present in IDL version, required to get identical mloc10 index
        min_tau10 = min(tau_scat10)
        mloc10 = np.where(tau_scat10 == min_tau10)

        tau_scat1 = abs(tau_scat - 1/1) # extra var not present in IDL version, required to get identical mloc1 index
        min_tau1 = min(tau_scat1)
        mloc1 = np.where(tau_scat1 == min_tau1)

        r_scat = rint[mloc10]
        r_tau1 = rint[mloc1]
        print('Scattering radius = {0}'.format(r_scat[0]))

        # Calculating the scattering rate
        # Scattering rate is by definition D_mu_mu = nu_s * (1 - mu^2)/2 for isotropic scattering
        # nu_s is already multiplied by c/omega_0

        dt_save = 2e7/omega # IDL version uses max(omega), but omega is one value so not required (also doesn't work in Python)
        # dt0 = 0.05/nu_s0
        exp_size = 1.25 * 30 / f_start
        dt0 = 0.01 * exp_size / c_r
        t, tplot = 0, 0

        print('Scattering time = {0} sec, dt save = {1}'.format(1/nu_s0, dt_save))

        # The block below determines the maximum distance
        # omega_pe / omega = 0.1

        #omega_pe_temp = abs(self.omega_pe(rint)/omega-0.01) # extra lines here required to get istop. IDL uses min(omega), but omega is single value
        #min_stop = min(omega_pe_temp)
        #istop = np.where(omega_pe_temp == min_stop)
        #r_stop = rint[istop]
        tau_scat_temp = abs(tau_scat-(0.1*30*np.pi/180/f_start)**2) # extra lines required to get mloc
        min_tau = min(tau_scat_temp) # value differs from IDL version by -0.2
        mloc = np.where(tau_scat_temp == min_tau)

        # angular spread less than (0.1 * 32 / f_start)^2 radians^2
        r_stop = max([rint[mloc], r_scat, r_init*2, 215])
        r_stop = min([r_stop, 215]) # not larger than 1 AU
        dt_save = (r_stop - r_init)/c_r/10

        logfmt = LogFormatter(base=10.0, labelOnlyBase=True)
        ax[1,0].loglog(self.omega_pe(rint)/2/np.pi/1e6, 1/self.nu_scat(rint,omega,qeps2_scaling)) # IDL version uses min(omega), but since omega is single value, not needed (also doesn't work in python)
        ax[1,0].loglog(self.omega_pe(rint)/2/np.pi/1e6, (rint-r_init)/c_r, '--', color = 'grey', linewidth = 0.8)
        ax[1,0].set_xlim(self.omega_pe(r_init)/2/np.pi/1e6, self.omega_pe(r_stop)/2/np.pi/1e6) # this doesn't constrain axes to next order of mag. as in IDL, so using next line manually
        ax[1,0].set_xlim(self.omega_pe(r_init)/2/np.pi/1e6,self.omega_pe(r_stop)/2/np.pi/1e6)
        #ax[1,0].set_ylim(1e-2, 1e8)
        ax[1,0].xaxis.set_major_formatter(logfmt)
        ax[1,0].set_xlabel('Plasma Frequency [MHz]')
        ax[1,0].set_ylabel('Scattering time [seconds]')

        ax[1,1].loglog(rint, tau_scat)
        ax[1,1].loglog([r_stop, r_stop], [1e-1, 10], ':', color = 'grey', linewidth = 0.8)
        ax[1,1].loglog([r_init, r_init], [1e-1, 10], ':', color = 'grey', linewidth = 0.8)
        ax[1,1].loglog([r_tau1, r_tau1], [1e-1, 10], '--', color = 'grey', linewidth = 0.8)
        ax[1,1].set_xlim(1, r_stop)
        ax[1,1].set_ylim(1e-1, 10)
        ax[1,1].yaxis.set_major_formatter(logfmt)
        ax[1,1].xaxis.set_major_formatter(logfmt)
        ax[1,1].xaxis.set_minor_formatter(NullFormatter())
        ax[1,1].set_xlabel('Distance [R$_\odot$]')
        ax[1,1].set_ylabel('Scattering depth $\\tau$')

        print('Stopping distance for calculations = {0} R_sun'.format(r_stop))

        plt.tight_layout()
        #plt.show()

        # Manually changing
        # r_stop = 3
        # while (t <= tmax) and (len(r) >= 10):

        while len(r) >= N/200:

            r = np.sqrt(rx**2 + ry**2 + rz**2)
            kc = np.sqrt(kx**2 + ky**2 + kz**2)
            omega = np.sqrt(self.omega_pe(np.sqrt(rx**2 + ry**2 + rz**2))**2 + kc**2)

            nu_s = self.nu_scat(r, omega, qeps2_scaling)
            nu_s = nu_s0 * (nu_s > nu_s0) + nu_s * (nu_s <= nu_s0)
                
            g0 = np.sqrt(nu_s * kc**2)
            # print(1/max(nu_s))
            dt_ref = min(abs(kc/(self.domega_dr(r) * c_r)/20))
            dt_dr = min(r/(c_r/omega0*kc))/20
            # mean scattering time in seconds
            dt = min([0.1/max(nu_s), dt0, dt_ref, dt_dr])
            
            # distance measured in solar radii. These values differ due to random function
            # used to generate kc, kx etc
            drx_dt = c_r/omega*kc*kx/kc
            dry_dt = c_r/omega*kc*ky/kc
            drz_dt = c_r/omega*kc*kz/kc

            rand = np.random.normal(loc = 0, scale = 1, size = (3,len(rz))) # emulating IDL function randomn with mean = 0, stdev = 1

            wx = rand[2,:] * np.sqrt(dt)#*np.sqrt(3/(2+anis**2))
            wy = rand[1,:] * np.sqrt(dt)#*np.sqrt(3/(2+anis**2))
            wz = rand[0,:] * np.sqrt(dt)#*np.sqrt(3*anis**2/(2+anis**2))

            # calculates positions of the photons in spherical coordinates
            fi0 = np.arctan2(ry, rx)
            sintheta0 = np.sqrt(1 - rz**2/r**2)
            costheta0 = rz/r

            #spherical coordinates with x to the Earth
            r_sph, th_sph, fi_sph = cartesian_to_spherical(rz, rx, ry)
            r_sph = r_sph.value # remove astropy units
            th_sph = th_sph.value
            fi_sph = fi_sph.value

            BxE, ByE, BzE = self.parker_Bxyz_Earth(r_sph, np.pi/2-th_sph, fi_sph)

            BB = np.sqrt(BxE**2 + ByE**2 + BzE**2)
            fi = -np.arctan2(BzE,ByE) # negative sign to match IDL output atan(BzE,ByE)
            sintheta = np.sqrt(1.-BxE**2/BB**2)
            costheta = BxE/BB

            #################################################################
            # Rotate k to r-vector aligned coordinate system

            kcx = -kx*np.sin(fi) + ky*np.cos(fi)
            kcy = -kx*costheta*np.cos(fi) - ky*costheta*np.sin(fi) + kz*sintheta
            kcz = kx*sintheta*np.cos(fi) + ky*sintheta*np.sin(fi) + kz*costheta

            kw = wx*kcx + wy*kcy + wz*kcz*anis
            akc = np.sqrt(kcx**2 + kcy**2 + kcz**2 * anis**2)
            z_asym = (asym * (kcz > 0) + (2 - asym) * (kcz <= 0)) * (kc/akc)**2

            a_perp = nu_s * z_asym * kc / akc**3 * (-(1 + anis**2) * akc**2 + 3 * anis**2 * (anis**2 - 1) * kcz**2) * anis
            a_par = nu_s * z_asym * kc / akc**3 * ((-3 * anis**4 + anis**2) * akc**2 + 3 * anis**4 * (anis**2 - 1) * kcz**2) * anis
            a_g0 = g0 * np.sqrt(z_asym * anis)

            kcx = kcx + a_perp*kcx*dt + a_g0*(wx - kcx*kw/akc**2)
            kcy = kcy + a_perp*kcy*dt + a_g0*(wy - kcy*kw/akc**2)
            kcz = kcz + a_par*kcz*dt + a_g0*(wz - kcz*kw*anis/akc**2) * anis

            #   print('dk^2 = ', min((kcx1)*kcx1+(kcy1)*kcy1+(kcz1)*kcz1), max((kcx1)*kcx1+(kcy1)*kcy1+(kcz1)*kcz1), min(kcx*kcx+kcy*kcy+kcz*kcz), max(kcx*kcx+kcy*kcy+kcz*kcz)

            #   kc_old2new = kc/np.sqrt(kx**2+ky**2+kz**2)
            #   kx = kx * kc_old2new
            #   ky = ky * kc_old2new
            #   kz = kz * kc_old2new
            #   kc = np.sqrt(kx**2+ky**2+kz**2)
            
            #   ax = c_s**2 * omega**2 / kc**2
            #   ay = 0
            #   az = 0
            #   kcx = kcx - nu_s*(1 + ax)*kcx*dt + g0*(wx - kcx*kw/kc**2 * (1 - np.sqrt(2*ax)))#*np.sqrt(np.abs(ay*az))
            #   kcy = kcy - nu_s*(1 + ay)*kcy*dt + g0*(wy - kcy*kw/kc**2 * (1 - np.sqrt(2*ay)))#*np.sqrt(np.abs(ax*az))
            #   kcz = kcz - nu_s*(1 + az)*kcz*dt + g0*(wz - kcz*kw/kc**2 * (1 - np.sqrt(2*az)))#*np.sqrt(np.abs(ax*ay))
            #   kz = kz + g0*np.abs(wz)
            #   kx = kx + g0*(-wx*np.sin(fi) + (wz*sintheta - wy*costheta)*np.cos(fi))
            #   ky = ky + g0*(wx*np.cos(fi) + (wz*sintheta - wy*costheta)*np.sin(fi))
            #   kz = kz + g0*(wz*costheta + wy*sintheta)
            #   Thejappa approach

            kx = -kcx*np.sin(fi) - kcy*costheta*np.cos(fi) + kcz*sintheta*np.cos(fi)
            ky = kcx*np.cos(fi) - kcy*costheta*np.sin(fi) + kcz*sintheta*np.sin(fi)
            kz = kcy*sintheta + kcz*costheta

            # After scattering rotate back from r-aligned to XYZ coordinate system
            #################################################################

            kc_norm = np.sqrt(kx**2 + ky**2 + kz**2)
            kx = kx*kc/kc_norm
            ky = ky*kc/kc_norm
            kz = kz*kc/kc_norm

            dk_dt = (self.omega_pe(r)/omega) * self.domega_dr(r) * c_r # differs from IDL version by 0.1

            kx = kx - dk_dt * dt * (rx/r)
            ky = ky - dk_dt * dt * (ry/r)
            kz = kz - dk_dt * dt * (rz/r)

            if np.isnan(kx).any() == True:
                exit()

            # distance measured in solar radii
            rx = rx + drx_dt * dt
            ry = ry + dry_dt * dt
            rz = rz + drz_dt * dt

            # to conserve the frequency
            kc_new_old = np.sqrt(omega**2 - self.omega_pe(np.sqrt(rx**2 + ry**2 + rz**2))**2) / np.sqrt(kx**2 + ky**2 + kz**2)
            kx = kx * kc_new_old
            ky = ky * kc_new_old
            kz = kz * kc_new_old

            # kc_norm = np.sqrt(kx**2 + ky**2 + kz**2)
            nu_e = 2.91e-6 * self.density_r(r) * 20/Te**1.5 * self.omega_pe(r)**2 / omega**2
            #nu_col = np.pi * (4.8e-10)**4 * self.density_r(r) * 20/np.sqrt(9.1e-28 * (11605*Te*1.38e-16)**3) * self.omega_pe(r)**2 / omega**2

            # optical depth of the photons
            tau = tau + nu_e * dt
            #frac = (1.5e7/3e10)**2
            #kc = kc - nu_s * kc * dt * frac + np.sqrt(g0*frac) * np.random.uniform(0, 1, size = len(rz)) * np.sqrt(dt)

            t = t + dt

            if max(tau) >= 4.61: # exp(-4.61) = 0.01
                # retaining photons that are not absorbed
                not_absorbed = np.where((tau - min(tau)) <= 4.61)[0] # exp(-4.61) = 0.01       
                rx, ry, rz = rx[not_absorbed], ry[not_absorbed], rz[not_absorbed]
                kx, ky, kz = kx[not_absorbed], ky[not_absorbed], kz[not_absorbed]
                tau = tau[not_absorbed]
                #exit()

            #rNone = np.where(r >= r_stop, r, None) # provides array of None where r <= r_stop and r where r >= r_stop
            #if len([x for x in rNone if x is not None]) >= 1: # if number of elements in r >= r_stop is >= 1. Use ths if using above line
            if len(np.where(r >= r_stop)[0]) >= 1: # if number of elements in r >= r_stop is >= 1

                if 'rr' in dir():
                
                    tt = np.append(tt, np.zeros(len(np.where(r >= r_stop)[0])) + t)
                    kxx, kx = np.append(kxx, kx[np.where(r >= r_stop)[0]]), kx[np.where(r < r_stop)[0]]
                    kyy, ky = np.append(kyy, ky[np.where(r >= r_stop)[0]]), ky[np.where(r < r_stop)[0]]
                    kzz, kz = np.append(kzz, kz[np.where(r >= r_stop)[0]]), kz[np.where(r < r_stop)[0]]

                    rxx, rx = np.append(rxx, rx[np.where(r >= r_stop)[0]]), rx[np.where(r < r_stop)[0]]
                    ryy, ry = np.append(ryy, ry[np.where(r >= r_stop)[0]]), ry[np.where(r < r_stop)[0]]
                    rzz, rz = np.append(rzz, rz[np.where(r >= r_stop)[0]]), rz[np.where(r < r_stop)[0]]

                    rxx0, rx0 = np.append(rxx0, rx0[np.where(r >= r_stop)[0]]), rx0[np.where(r < r_stop)[0]]
                    ryy0, ry0 = np.append(ryy0, ry0[np.where(r >= r_stop)[0]]), ry0[np.where(r < r_stop)[0]]
                    rzz0, rz0 = np.append(rzz0, rz0[np.where(r >= r_stop)[0]]), rz0[np.where(r < r_stop)[0]]

                    tau0, tau = np.append(tau0, tau[np.where(r >= r_stop)[0]]), tau[np.where(r < r_stop)[0]]
                    kk, kc = np.append(kk, kc[np.where(r >= r_stop)[0]]), kc[np.where(r < r_stop)[0]]
                    rr, r = np.append(rr, r[np.where(r >= r_stop)[0]]), r[np.where(r < r_stop)[0]]

                else:
                    
                    # tt, rr, kk, rxx, ryy, rzz, kxx, kyy, kzz = 0          
                    tt = np.zeros(len(np.where(r >= r_stop)[0])) + t

                    kxx, kx = kx[np.where(r >= r_stop)[0]], kx[np.where(r < r_stop)[0]]
                    kyy, ky = ky[np.where(r >= r_stop)[0]], ky[np.where(r < r_stop)[0]]
                    kzz, kz = kz[np.where(r >= r_stop)[0]], kz[np.where(r < r_stop)[0]]

                    rxx, rx = rx[np.where(r >= r_stop)[0]], rx[np.where(r < r_stop)[0]]
                    ryy, ry = ry[np.where(r >= r_stop)[0]], ry[np.where(r < r_stop)[0]]
                    rzz, rz = rz[np.where(r >= r_stop)[0]], rz[np.where(r < r_stop)[0]]

                    rxx0, rx0 = rx0[np.where(r >= r_stop)[0]], rx0[np.where(r < r_stop)[0]]
                    ryy0, ry0 = ry0[np.where(r >= r_stop)[0]], ry0[np.where(r < r_stop)[0]]
                    rzz0, rz0 = rz0[np.where(r >= r_stop)[0]], rz0[np.where(r < r_stop)[0]]

                    tau0, tau = tau[np.where(r >= r_stop)[0]], tau[np.where(r < r_stop)[0]]
                    kk, kc = kc[np.where(r >= r_stop)[0]], kc[np.where(r < r_stop)[0]]
                    rr, r = r[np.where(r >= r_stop)[0]], r[np.where(r < r_stop)[0]]

                    # important: rr or rz should last line!
            
            if t >= (tplot + dt_save):

                tplot = tplot + dt_save
                pp = np.linspace(0, 6.28319, num = 41) # findgen(41)/40*!pi*2

                #fig2 = plt.figure(figsize = (5,5))
                #fig2.plot(rz0, rx0, marker = 'x', color = 'k', linewidth = 0.5, markersize = 5)
                #fig2.plot(r_stop*np.cos(pp), r_stop*np.sin(pp), linestyle = (0, (1, 10)), color = 'k', linewidth = 0.5)
                #dt_save = dt

                fig3, ax3 = plt.subplots(2, 1, figsize = (6,10))
                obs_i0 = np.where(kz/kc > 0.9)[0]

                dir_beam = np.average(np.arctan2(kx,kz))

                ax3[0].plot(rx0, ry0, marker = 'x', color = 'k', linewidth = 0.5, markersize = 5)
                ax3[0].plot(np.cos(pp), np.sin(pp), linewidth = 1, color = 'gold')
                ax3[0].scatter([rx], [ry], marker = '.', color = 'green', linewidth = 0.5, s = 1)
                ax3[0].plot(r_stop*np.cos(pp), r_stop*np.sin(pp), linestyle = (0, (1, 10)), color = 'k', linewidth = 0.5)
                ax3[0].plot(rx[obs_i0], ry[obs_i0], linestyle = (0, (1, 10)), color = 'red', linewidth = 0.5)
                ax3[0].plot([0,0], [-r_stop, r_stop], linestyle = (0, (1, 10)), color = 'k', linewidth = 0.5)
                ax3[0].plot([-r_stop, r_stop], [0,0], linestyle = (0, (1, 10)), color = 'k', linewidth = 0.5)
                ax3[0].set_xlim(-round(r_stop+5.1,-1), round(r_stop+5.1,-1))
                ax3[0].set_ylim(-round(r_stop+5.1,-1), round(r_stop+5.1,-1))
                ax3[0].set_aspect('equal', 'box')
                ax3[0].set_xlabel('X [R$_\odot$]')
                ax3[0].set_ylabel('Y [R$_\odot$]')

                #rp = 0.1 * (r_stop - r_init)
                rp = 0.333 * (r_stop - r_init)
                rzc = np.sqrt(r_stop**2 - rp**2)
                ax3[1].plot(rz0, rx0, marker = 'x', color = 'k', linewidth = 0.5, markersize = 5)
                ax3[1].plot(np.cos(pp), np.sin(pp), linewidth = 1, color = 'gold')
                ax3[1].plot(r_stop*np.cos(pp), r_stop*np.sin(pp), linestyle = (0, (1, 10)), color = 'k', linewidth = 0.5)
                ax3[1].plot([0,0], [r_stop, -r_stop], linestyle = (0, (1, 10)), color = 'k', linewidth = 0.5)
                ax3[1].plot([0, r_stop], [0,0], linestyle = (0, (1, 10)), color = 'k', linewidth = 0.5)
                ax3[1].plot([rzc, rzc], [-rp, rp], color = 'c', linewidth = 0.5)
                ax3[1].scatter([rz], [rx], marker = '.', color = 'green', linewidth = 0.5, s = 1)
                ax3[1].plot(rz[obs_i0], rx[obs_i0], linestyle = (0, (1, 10)), color = 'red', linewidth = 0.5)
                
                # oploting emission direction
                #ax3[1].plot([0, np.cos(dir_beam)*r_stop], [0, np.sin(dir_beam)*r_stop], linestyle = 'dashed', color = 'orange')

                rp = np.linspace(1, 215, 100)
                BxE, ByE, BzE = self.parker_Bxyz_Earth(rp, np.pi/2, np.pi/2)
                ps_ang = np.arctan2(ByE,BxE)    
                ax3[1].plot(rp*np.cos(ps_ang), rp*np.sin(ps_ang), color = 'red', linestyle = 'dashed', linewidth = 1)
                BxE, ByE, BzE = self.parker_Bxyz_Earth(rp, np.pi/2, np.pi/4)
                ps_ang = np.arctan2(ByE,BxE)    
                ax3[1].plot(rp*np.cos(ps_ang), rp*np.sin(ps_ang), color = 'red', linestyle = 'dashed', linewidth = 1)
                BxE, ByE, BzE = self.parker_Bxyz_Earth(rp, np.pi/2, 0)
                ps_ang = np.arctan2(ByE,BxE)    
                ax3[1].plot(rp*np.cos(ps_ang), rp*np.sin(ps_ang), color = 'red', linestyle = 'dashed', linewidth = 1)
                BxE, ByE, BzE = self.parker_Bxyz_Earth(rp, np.pi/2, -np.pi/4)
                ps_ang = np.arctan2(ByE,BxE)    
                ax3[1].plot(rp*np.cos(ps_ang), rp*np.sin(ps_ang), color = 'red', linestyle = 'dashed', linewidth = 1)
                
                ax3[1].set_xlim(-r_stop-2, r_stop+2)
                ax3[1].set_ylim(-r_stop-2, r_stop+2)
                ax3[1].set_aspect('equal', 'box')
                ax3[1].set_xlabel('Z [R$_\odot$]')
                ax3[1].set_ylabel('X [R$_\odot$]')

                fig3.tight_layout()   
                
                # uncomment these 3 lines to see scattering plots
                plt.draw()
                plt.pause(1.5)
                plt.close(fig = fig3)

                # if rather stop at each plot and only proceed when you close the window, comment above 3 lines and use:
                #plt.show()

                print('t = {0}, <r> = {1}, N_phot = {2}, <weight> = {3}'.format(t, np.mean(r), len(tau), np.mean(np.exp(-tau))))
                print('F = ', np.mean(np.sqrt(self.omega_pe(r)**2 + kc**2)/omega[0:len(kc)]), min(np.sqrt(self.omega_pe(r)**2 + kc**2)/omega[0:len(kc)]), max(np.sqrt(self.omega_pe(r)**2 + kc**2)/omega[0:len(kc)]))
                # had to add subscript range to omega as python can not broadcast mismatched array sizes, whereas IDL just uses the elements that match the min array size such
                # that the answer is size of smallest array used

                r_tab = r_init + (r_stop - r_init) * np.linspace(0, 900, num = 900)/(900-1)
                d_delay = integrate.simps(1/np.sqrt(1 - self.omega_pe(r_tab)**2 / omega0**2), x = r_tab) / c_r
                tau_exp = integrate.simps(2.91e-6 * self.density_r(r_tab) * 20/Te**1.5 * self.omega_pe(r_tab)**2 / omega0**2 / np.sqrt(1 - self.omega_pe(r_tab)**2 / omega0**2), x = r_tab) / c_r

                time.sleep(0.01)


        # rr, kk, t_free have 995 elements, IDL versions have 996. Values all similar.
        r_free = np.sqrt((rxx - rxx0)**2 + (ryy - ryy0)**2 + (rzz - rzz0)**2)
        rr = np.sqrt(rxx**2 + ryy**2 + rzz**2)
        kk = np.sqrt(kxx**2 + kyy**2 + kzz**2)
        t_free = r_free/c_r
        print('Delta t = ', d_delay - min(t_free)) # differs to IDL version due to min(t_free)

        cosx, cosy, cosz = kxx/kk, kyy/kk, kzz/kk
        fi = np.arctan2(cosy, cosx)

        # cosx0, cosy0, cosz0 = (rxx - rxx0)/r_free, (ryy - ryy0)/r_free, (ryy - ryy0)/r_free
        # xa_free = np.arctan((rxx - rxx0)/(rzz - rzz0))
        # ya_free = np.arctan((ryy - ryy0)/(rzz - rzz0))

        # collecting photons at infinity into angle (far-away source)
        # correct method for far away source
        obs_i = np.where((kzz/kk >= 0.9) & (kzz/kk <= 1.0))
        ww = np.exp(-tau0[obs_i])
        x0, y0 = r_free*kxx/kk, r_free*kyy/kk
        x_im, y_im = rxx[obs_i] - x0[obs_i], ryy[obs_i] - y0[obs_i]

        #x0, y0 = np.sqrt((rzz - rz0)**2 + (rxx - rx0)**2 + (ryy - ry0)**2) * cosx, np.sqrt((rzz - rz0)**2 + (rxx - rx0)**2 + (ryy - ry0)**2) * cosy
        #x0, y0 = (r_stop - r_init)*np.tan(np.arccos(cosr_stop)), (r_stop - r_init)*np.tan(np.arccos(cosr_stop))
        #x0, y0 = (rzz - r_init)*np.sqrt((rzz - r_init)**2 + rxx**2 + ryy**2)*cosx, max(rzz)*cosy

        #x0, y0 = (rzz - rzz0)*kxx/kzz, (rzz - rzz0)*kyy/kzz
        #x0, y0 = kxx/kk*c_r*tt, kyy/kk*c_r*tt

        # time integrated source FWHM size centroids and shifts in R_sun units
        xc, yc = np.mean(x_im * ww)/np.mean(ww), np.mean(y_im * ww)/np.mean(ww)
        sx = np.sqrt(np.mean(x_im**2 * ww - xc**2 * ww)/np.mean(ww)) * 2.355 # differs to IDL by 0.01
        sy = np.sqrt(np.mean(y_im**2 * ww - yc**2 * ww)/np.mean(ww)) * 2.355 # differs to IDL by 0.007

        # collect photons at 1 AU

        z_shift = (au - rzz)/cosz
        x_au = rxx + cosx*z_shift
        y_au = ryy + cosy*z_shift
        z_au = rzz + cosz*z_shift
        rp_1au = np.sqrt(x_au**2 + y_au**2)
        obs_im = np.where(rp_1au <= 0.333*(r_stop - r_init)*au/r_stop)

        xa_im = np.arctan(kxx[obs_im]/kzz[obs_im]) - np.arctan((x_au[obs_im] - rxx0[obs_im])/(z_au[obs_im] - rzz0[obs_im])) + np.arctan(-rxx0[obs_im]/(z_au[obs_im] - rzz0[obs_im]))
        ya_im = np.arctan(kyy[obs_im]/kzz[obs_im]) - np.arctan((y_au[obs_im] - ryy0[obs_im])/(z_au[obs_im] - rzz0[obs_im])) + np.arctan(-ryy0[obs_im]/(z_au[obs_im] - rzz0[obs_im]))
        wa = np.exp(-tau0[obs_im])

        # angle imaging for closely located sources
        xac, yac = np.mean(xa_im*wa)/np.mean(wa), np.mean(ya_im*wa)/np.mean(wa)
        sax = np.sqrt(np.mean(xa_im**2*wa - xac**2 * wa)/np.mean(wa))*2.355
        say = np.sqrt(np.mean(ya_im**2*wa - yac**2 * wa)/np.mean(wa))*2.355 # differs to IDL by 0.005
        print('X/Y angle sizes [degrees]: ', sax*180/np.pi, say*180/np.pi)

        err_xc = sx/np.sqrt(len(x_im) + 1e-8)/2.355
        err_yc = sy/np.sqrt(len(y_im) + 1e-8)/2.355

        err_sx = sx*np.sqrt(2)/np.sqrt(len(x_im) + 1e-8)
        err_sy = sy*np.sqrt(2)/np.sqrt(len(y_im) + 1e-8)

        # without absorption
        xc_noa, yc_noa = np.mean(x_im), np.mean(y_im)
        sx_noa = np.sqrt(np.mean(x_im**2 - xc_noa**2))*2.355
        sy_noa = np.sqrt(np.mean(y_im**2 - yc_noa**2))*2.355

        xlimits = [-1.1*max([xc*2, sx*2, 1]), 1.1*max([xc*2, sx*2, 1])]
        ylimits = [-1.1*max([yc*2, sy*2, 1]), 1.1*max([yc*2, sy*2, 1])]

        # angular distribution for all photons  at r > r_stop
        #dtbin = (np.mean(tt - t_free) - min(tt - t_free))/10
        #t_max = np.mean(tt - t_free)*2
        #pdf_mu, tbins = np.asarray((np.histogram(tt - t_free, weights = np.exp(-tau0*kzz/kk), bins = np.arange(0, t_max, dtbin), density = True)))# / (np.asarray(np.histogram(tt - t_free, weights = np.exp(-tau0), bins = np.arange(0, t_max, dtbin), density = True))+1e-10)
        #pdf_mu1, tbins1 = np.asarray((np.histogram(tt - t_free, weights = kzz/kk, bins = np.arange(0, t_max, dtbin), density = True)))# / (np.asarray(np.histogram(tt - t_free, weights = 1.0, bins = np.arange(0, t_max, dtbin), density = True))+1e-10)

        # Angular distribution of ALL photons
        cosA, w2 = np.concatenate([(kzz/kk - 1), -(kzz/kk - 1)]), np.concatenate([np.exp(-tau0), np.exp(-tau0)])
        muc = np.mean(w2 * cosA) / np.mean(w2)
        smu = np.sqrt(np.mean((cosA - muc)**2 * w2) / np.mean(w2))
        thetac = np.mean(w2 * np.arccos(cosA)) / np.mean(w2)
        s_theta = np.sqrt(np.mean((np.arccos(cosA) - thetac)**2 * w2) / np.average(w2)) * 180/np.pi
        print('Direction = {0} degrees and FWHM angle = {1} degrees'.format(thetac*180/np.pi, s_theta))

        mubin = 0.01
        pdf_mu, mubins = np.histogram((kzz/kk), weights = np.exp(-tau0), bins = np.arange(0, 1, mubin))
        mu_halfw = mubins[np.where(pdf_mu/max(pdf_mu) >= 0.5)].min() + mubin/2
        mubin = 0.05
        pdf_mu, mubins = np.histogram((kzz/kk), weights = np.exp(-tau0), bins = np.arange(0, 1, mubin))

        fig4, ax4 = plt.subplots(2, 3, figsize = (10,6))

        ax4[0,0].step(mubins[:-1] + mubin/2, pdf_mu/max(pdf_mu), color = 'k')
        ax4[0,0].plot([mu_halfw, mu_halfw], [0, 1], 'k--', linewidth = 0.5)
        ax4[0,0].set_xlabel('mu = kzz/kk')
        ax4[0,0].set_xlim(0, 1)
        ax4[0,0].set_ylim(0, 1)
        print('mu_HW = {0}, Theta_HW = {1}'.format(mu_halfw, np.arccos(mu_halfw)*180/np.pi))

        # Only obs_i photons
        deltaT = tt[obs_i] - t_free[obs_i]
        dtbin = (np.mean(deltaT) - min(deltaT))/23
        maxT = np.mean(deltaT)*2.3
        Tpdf, tbins = np.histogram(deltaT, weights = ww, bins = np.arange(0, maxT, dtbin))
        Tpdf_noa, tbins_noa = np.histogram(deltaT, weights = ww, bins = np.arange(0, maxT, dtbin)) # extra line required to get the unweighted histogram
        max_Tpdf_temp = Tpdf/max(Tpdf)
        max_Tpdf = max(max_Tpdf_temp)
        i_maxT = np.where(max_Tpdf_temp == max_Tpdf)
        t_max = tbins[i_maxT] + dtbin/2
        t_HM_before = tbins[np.where(Tpdf/max(Tpdf) >= 0.5)].min() + dtbin/2
        t_HM_after = tbins[np.where(Tpdf/max(Tpdf) >= 0.5)].max() + dtbin/2
        t_1e_after = tbins[np.where(Tpdf/max(Tpdf) >= np.exp(-1))].max() + dtbin/2

        # tbins = 0.5*(tbins[1:]+tbins[:-1])  # use this to plot bar centres rather than edges. Avoids having to use tbins[:-1] in step plot

        ax4[0,1].step(tbins[:-1] + dtbin/2, Tpdf/max(Tpdf), color = 'k')
        ax4[0,1].step(tbins_noa[:-1] + dtbin/2, Tpdf_noa/max(Tpdf_noa),'--', color = 'fuchsia')
        ax4[0,1].plot([t_max, t_max], [min(max_Tpdf_temp), max(max_Tpdf_temp)], 'k--', linewidth = 0.5)
        ax4[0,1].plot([t_HM_before, t_HM_before], [min(max_Tpdf_temp), max(max_Tpdf_temp)], 'k:', linewidth = 0.5)
        ax4[0,1].plot([t_HM_after, t_HM_after], [min(max_Tpdf_temp), max(max_Tpdf_temp)], 'k:', linewidth = 0.5)
        ax4[0,1].set_xlabel('Time [s]')
        ax4[0,1].set_xlim(0, maxT)
        ax4[0,1].set_ylim(0, 1)

        print('Tmax = {0}, Delay = {1}, HM Before = {2}, HM After = {3}'.format(t_max[0], t_max[0], t_HM_before, t_HM_after))
        r_proj = r_init * np.sin(rtheta[0])

        ax4[0,2].scatter(x_im, y_im, marker = '.', color = 'black', linewidth = 0.5, s = 1)
        ax4[0,2].plot(np.cos(pp), np.sin(pp), linewidth = 1, color = 'gold')
        ax4[0,2].plot(r_init*np.cos(pp), r_init*np.sin(pp), 'k--', linewidth = 0.5)
        ax4[0,2].plot(r_proj*np.cos(pp), r_proj*np.sin(pp), 'r--', linewidth = 0.5)
        ax4[0,2].plot(r_scat*np.sin(rtheta[0])*np.cos(pp), r_scat*np.sin(rtheta[0])*np.sin(pp), '--', color = 'fuchsia')
        ax4[0,2].plot([xc, xc], [yc, yc], marker = 'x', color = 'green')
        ax4[0,2].set_xlim(xlimits[0]-0.4, xlimits[1]+0.4)
        ax4[0,2].set_ylim(xlimits[0]-0.4, xlimits[1]+0.4)
        ax4[0,2].set_aspect('equal', 'box')

        print('With absorption:')
        print('X-shift [R_s] = {0} +/- {1}   Y-shift [R_s] = {2} +/- {3}'.format(xc - r_init*np.sin(rtheta[0]), err_xc, yc - 0, err_yc))
        print('Xsize [R_s] = {0} +/- {1}   Ysize [R_s] = {2} +/- {3}'.format(sx, err_sx, sy, err_sy))
        print('Without absorption:')
        print('X-shift [R_s] = {0}   Y-shift [R_s] = {1}'.format(xc_noa - r_init*np.sin(rtheta[0]), yc_noa - 0))
        print('Xsize [R_s] = {0}   Ysize [R_s] = {1}'.format(sx_noa, sy_noa))

        # with absorption
        xc_time, xc_bins_ignore = np.histogram(deltaT, weights = x_im*ww, bins = np.arange(0, maxT, dtbin))
        xc_time = xc_time/np.histogram(deltaT, weights = ww, bins = np.arange(0, maxT, dtbin))[0] + 1e-10
        yc_time, yc_bins_ignore = np.histogram(deltaT, weights = y_im*ww, bins = np.arange(0, maxT, dtbin))
        yc_time = yc_time/np.histogram(deltaT, weights = ww, bins = np.arange(0, maxT, dtbin))[0] + 1e-10
        xs2_time, xc_bins_ignore = np.histogram(deltaT, weights = x_im**2 * ww, bins = np.arange(0, maxT, dtbin))
        xs2_time = xs2_time/np.histogram(deltaT, weights = ww, bins = np.arange(0, maxT, dtbin))[0] + 1e-10
        ys2_time, xc_bins_ignore = np.histogram(deltaT, weights = y_im**2 * ww, bins = np.arange(0, maxT, dtbin))
        ys2_time = ys2_time/np.histogram(deltaT, weights = ww, bins = np.arange(0, maxT, dtbin))[0] + 1e-10

        xs_time = np.sqrt(xs2_time - xc_time**2)*2.355
        ys_time = np.sqrt(ys2_time - yc_time**2)*2.355
        err_xc_time = xs_time/np.sqrt(Tpdf_noa + 1e-8)/2.355
        err_yc_time = ys2_time/np.sqrt(Tpdf_noa + 1e-8)/2.355
        err_xs_time = xs_time/np.sqrt(2*Tpdf_noa + 1e-8)
        err_ys_time = ys_time/np.sqrt(2*Tpdf_noa + 1e-8)

        #xc_time = xc_time/(Tpdf + 1e-10)
        #yc_time = yc_time/(Tpdf + 1e-10)

        ax4[1,0].step(tbins[:-1] + dtbin/2, xc_time, color = 'k', linewidth = 0.5)
        ax4[1,0].errorbar(tbins[:-1] + dtbin/2, xc_time, yerr = err_xc_time, fmt = 'none', capsize = 2, color = 'k', linewidth = 0.5)
        ax4[1,0].step(tbins[:-1] + dtbin/2, yc_time, color = 'green', linewidth = 0.5)
        ax4[1,0].errorbar(tbins[:-1] + dtbin/2, yc_time, yerr = err_yc_time, fmt = 'none', capsize = 2, color = 'green', linewidth = 0.5)

        ax4[1,0].plot([min(tbins + dtbin/2), max(tbins + dtbin/2)], [xc, xc], '--', color = 'k', linewidth = 0.5)
        ax4[1,0].plot([min(tbins + dtbin/2), max(tbins + dtbin/2)], [yc, yc], '--',  color = 'green', linewidth = 0.5)
        ax4[1,0].plot([min(tbins + dtbin/2), max(tbins + dtbin/2)], [min(rxx0), min(rxx0)], color = 'fuchsia')
        ax4[1,0].plot([min(tbins + dtbin/2), max(tbins + dtbin/2)], [min(ryy0), min(ryy0)], color = 'fuchsia')
        ax4[1,0].plot([tbins[i_maxT] + dtbin/2, tbins[i_maxT] + dtbin/2], [min([0, xc, yc]), max([0, xc, yc])], '--', color = 'k', linewidth = 0.5)

        ax4[1,0].set_xlabel('Time [s]')
        ax4[1,0].set_ylabel('X/Y Centroids [R$_s$]')
        ax4[1,0].set_xlim(0, maxT)
        ax4[1,0].set_ylim(min([0, xc, yc])*2, max([0, xc, yc])*2)

        ax4[1,1].step(tbins[:-1] + dtbin/2, xs_time, color = 'k', linewidth = 0.5)
        ax4[1,1].errorbar(tbins[:-1] + dtbin/2, xs_time, yerr = err_xs_time, fmt = 'none', capsize = 2, color = 'k', linewidth = 0.5)
        ax4[1,1].step(tbins[:-1] + dtbin/2, ys_time, color = 'green', linewidth = 0.5)
        ax4[1,1].errorbar(tbins[:-1] + dtbin/2, ys_time, yerr = err_ys_time, fmt = 'none', capsize = 2, color = 'green', linewidth = 0.5)

        ax4[1,1].plot([min(tbins + dtbin/2), max(tbins + dtbin/2)], [sx, sx], '--', color = 'k', linewidth = 0.5)
        ax4[1,1].plot([min(tbins + dtbin/2), max(tbins + dtbin/2)], [sy, sy], '--', color = 'green', linewidth = 0.5)
        ax4[1,1].plot([tbins[i_maxT] + dtbin/2, tbins[i_maxT] + dtbin/2], [min([0, sx, sy]), max([0, sx, sy])], '--', color = 'k', linewidth = 0.5)

        ax4[1,1].set_xlabel('Time [s]')
        ax4[1,1].set_ylabel('FWHM X/Y Size [R$_s$]')
        ax4[1,1].set_xlim(0, maxT)
        ax4[1,1].set_ylim(min([0, sx, sy])*2, max([0, sx, sy])*2)


        # Plotting profiles along x and y
        r_max = max(np.concatenate((x_im, y_im))) * 1.2
        xbin = r_max/30
        xpdf, xloc = np.histogram(x_im, bins = np.arange(-r_max, r_max, xbin))
        ypdf, yloc = np.histogram(y_im, bins = np.arange(-r_max, r_max, xbin))
        ax4[1,2].step(xloc[:-1] + xbin/2, xpdf, color = 'k', linewidth = 0.5)
        ax4[1,2].step(yloc[:-1] + xbin/2, ypdf, color = 'green', linewidth = 0.5)
        ax4[1,2].set_xlim(xlimits[0], xlimits[1])
        ax4[1,2].set_ylim(0, max(np.concatenate((xpdf, ypdf))))

        print('FWHM source sizes (x & y) [degrees] = ', np.arctan(sx*0.5/(au - r_init))*360/np.pi, np.arctan(sy*0.5/(au - r_init))*360/np.pi)
        print('Source area (ellipse) [arcmin^2] = ', np.pi*sx*sy/4*16**2)
        print('Predicted FWHM source sizes [arcmin] =', np.arctan((rint[mloc1] - r_init)*0.5/(au - r_init))*360/np.pi)

        # print('X-size [arcmin] = ', np.sqrt(np.mean(x_im**2 * ww - xc**2 * ww)/np.mean(ww)))
        # print('Y-size [arcmin] = ', np.sqrt(np.mean(y_im**2 * ww - yc**2 * ww)/np.mean(ww)))

        plt.tight_layout()


        # Save data to fits file

        p_hdr = fits.Header()
        t_hdr = fits.Header()

        cols = []

        cols.append(fits.Column(name = 'title', format = '39A', array = ['Ray tracing results in the solar corona']))
        cols.append(fits.Column(name = 'r_init', format = 'E', array = [r_init]))
        cols.append(fits.Column(name = 'r_stop', format = 'E', array = [r_stop]))
        cols.append(fits.Column(name = 'f_pe0', format = 'E', array = [omega_pe0/1e6/2/np.pi]))
        cols.append(fits.Column(name = 'omega0', format = 'E', array = [omega_pe0]))
        cols.append(fits.Column(name = 'fpe_stop', format = 'E', array = omega/1e6/2/np.pi))
        cols.append(fits.Column(name = 'omega', format = 'E', array = omega/1e6/2/np.pi))
        cols.append(fits.Column(name = 't_free', format = 'E', array = t_free))
        cols.append(fits.Column(name = 'kxx', format = 'E', array = kxx))
        cols.append(fits.Column(name = 'kyy', format = 'E', array = kyy))
        cols.append(fits.Column(name = 'kzz', format = 'E', array = kzz))
        cols.append(fits.Column(name = 'kk', format = 'E', array = kk))
        cols.append(fits.Column(name = 'rxx', format = 'E', array = rxx))
        cols.append(fits.Column(name = 'ryy', format = 'E', array = ryy))
        cols.append(fits.Column(name = 'rzz', format = 'E', array = rzz))
        cols.append(fits.Column(name = 'rr', format = 'E', array = rr))
        cols.append(fits.Column(name = 'tt', format = 'E', array = tt))
        cols.append(fits.Column(name = 'tau0', format = 'E', array = tau0))
        cols.append(fits.Column(name = 'qeps2_scaling', format = 'E', array = [qeps2_scaling]))
        cols.append(fits.Column(name = 'asym', format = 'E', array = [asym]))
        cols.append(fits.Column(name = 'anis', format = 'E', array = [anis]))
        cols.append(fits.Column(name = 'rtheta', format = 'E', array = rtheta))
        cols.append(fits.Column(name = 'rfi', format = 'E', array = rfi))
        cols.append(fits.Column(name = 'f_ratio', format = 'E', array = [f_ratio]))
        cols.append(fits.Column(name = 'xcen', format = 'E', array = [xc]))
        cols.append(fits.Column(name = 'ycen', format = 'E', array = [yc]))
        cols.append(fits.Column(name = 'xc_err', format = 'E', array = [err_xc]))
        cols.append(fits.Column(name = 'yc_err', format = 'E', array = [err_yc]))
        cols.append(fits.Column(name = 'xsize', format = 'E', array = [sx]))
        cols.append(fits.Column(name = 'ysize', format = 'E', array = [sy]))
        cols.append(fits.Column(name = 'xs_err', format = 'E', array = [err_sx]))
        cols.append(fits.Column(name = 'ys_err', format = 'E', array = [err_sx]))
        cols.append(fits.Column(name = 'theta_HW', format = 'E', array = [np.arccos(mu_halfw)*180/np.pi]))
        cols.append(fits.Column(name = 't_peak', format = 'E', array = [t_max]))
        cols.append(fits.Column(name = 't_HM_before', format = 'E', array = [t_HM_before]))
        cols.append(fits.Column(name = 't_HM_after', format = 'E', array = [t_HM_after]))
        cols.append(fits.Column(name = 't_1e_after', format = 'E', array = [t_1e_after]))
        cols.append(fits.Column(name = 'terr', format = 'E', array = [dtbin]))

        fname = 'fpe' + (str(int(omega_pe0/1e3/2/np.pi))).zfill(7) + 'kHz_FE' + (str(int(omega0/1e3/2/np.pi))).zfill(7) + 'anis' + str(format(anis, '.2f')) + 'qeps2_s' + str(format(qeps2_scaling, '.3f')) + 'asym' + str(format(asym, '.2f')) + 'fr' + str(format(f_ratio, '.2f'))
        hdulist = fits.HDUList()
        hdulist.append(fits.PrimaryHDU(header = p_hdr))
        hdulist.append(fits.BinTableHDU.from_columns(cols, header=t_hdr, name='DATA'))
        hdulist.writeto(fname + '.fits', overwrite = True)
        hdulist.close()

        print('Data saved to file: ', fname)
        print('All completed .............................................................. OK')

        plt.show()
