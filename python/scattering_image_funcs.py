def anis_scat_image(**params):
    """
    Set the desired parameter kwargs, e.g.
    f=30, eps=1, anis=0.25, xmax=3000, th=30
    to generate a simulated 2D Gaussian representing a scattered image
    based on ansiotropic radio-wave propagation simulation results. The resulting contours show the
    50%, 70%, and 90% contour levels. The black cross denotes the initial source location in the sky-plane
    and the black cross shows the scattered image centroid.
    The parameters are defined as:
    f: Observed frequency in MHz.
    eps: Scattering strength factor (possible values of 0.5, 0.71, 1.0, 1.41 ,2.0).
    anis: Anisotropy factor of the density fluctuation spectrum (possible values of 0.19, 0.25, 0.33, 0.42) where 1.0 is isotropic.
    Max plot range: Adjust the maximum x, y range.
    Source heliocentric angle: Change the initial source azimuthal angle from -60 to 60 degrees (0 degrees is the disk centre). 
    """
    
    # Set parameters
    if len(params) == 0:
        f, eps, anis, xmax, th = 30, 1.0, 0.25, 3000, 30
    elif len(params) == 5:
        f, eps, anis, xmax, th = params.values()
    else:
        f=30 if 'f' not in params else params['f']
        eps=1.0 if 'eps' not in params else params['eps']
        anis=0.25 if 'anis' not in params else params['anis']
        xmax=3000 if 'xmax' not in params else params['xmax']
        th=30 if 'th' not in params else params['th']

    # Check data file exists
    from os.path import exists
    file = 'new_all_results_F_final.sav'
    file_exists = exists(file)
    if file_exists == True:
            
        # check for valid user input
        if (eps in [0.5,0.71,1.0,1.41,2.0]) and (anis in [0.19,0.25,0.33,0.42]) and (0.1 <= f <= 340): 

            import numpy as np
            import matplotlib.pyplot as plt
            import matplotlib.patches as patches
            from matplotlib.ticker import MultipleLocator, AutoMinorLocator
            from scipy.io import readsav
            from scipy.interpolate import interp1d
            from scipy.optimize import curve_fit

            def gauss2d(xy, I, sx, sy, xc, yc, T):
                """
                2D Gaussian
                """
                x, y = xy
                xprime = (x - xc)*np.cos(T) - (y - yc)*np.sin(T)
                yprime = (x - xc)*np.sin(T) + (y - yc)*np.cos(T)
                a = sx/2
                b = sy/2
                U = (xprime / a)**2 + (yprime / b)**2
                return I*np.exp(-U/2).ravel()

            # constants
            r_sun = 950

            # Load data
            data = readsav(file, python_dict=True, verbose=True)

            # parse ds data
            r_init = np.array(data['sims'][0][1])
            fmhz = np.array(data['sims'][0][2])
            r_shift = np.array(data['sims'][0][3])
            s_fwhm = np.array(data['sims'][0][4])
            eps0 = np.array(data['sims'][0][5])
            anis0 = np.array(data['sims'][0][7])
            t_1e_after = np.array(data['sims'][0][10])
            t_peak = np.array(data['sims'][0][11])

            th = np.deg2rad(th)

            # retrieve data for eps, alpha
            idx = np.where((eps0 == eps) & (anis0 == anis))
            fmhz = fmhz[idx]
            r_init = r_init[idx]
            r_shift = r_shift[idx]
            s_fwhm = s_fwhm[idx]
            t_peak = t_peak[idx]
            t_1e_after = t_1e_after[idx]

            # interpolate at f
            s_interp = interp1d(fmhz, s_fwhm)
            tpeak_interp = interp1d(fmhz, t_peak)
            t_1e_after_interp = interp1d(fmhz, t_1e_after)
            r_init_interp = interp1d(fmhz, r_init)
            r_shift_interp = interp1d(fmhz, r_shift)

            size = s_interp(f)*3600
            r_init= r_init_interp(f)
            r_shift = r_shift_interp(f)
            t_peak = tpeak_interp(f)
            decay_1e = t_1e_after_interp(f)
            decay_fwhm = decay_1e*np.sqrt(np.log(2))

            xmax = 3000
            ymax = xmax
            x0 = np.linspace(-xmax, xmax, 240)
            y0 = x0
            z0 = x0*0 + r_shift*r_sun
            xc0 = 0
            yc0 = 0
            zc0 = r_shift*r_sun
            z_init = r_init*r_sun
            xx, yy = np.meshgrid(x0, y0) # get 2D variables instead of 1D
            xy = gauss2d((xx, yy), 1, size, size, 0, 0, 0) # I, sx, sy, xc, yc, T
            xy = (xy - xy.min())/(xy.max() - xy.min()) # normalise 0-1

            # theta (azimuthal) rotation
            x = x0*np.cos(th) + z0*np.sin(th)
            y = y0
            z = -x0*np.sin(th) + z0*np.cos(th)

            xc = xc0*np.cos(th) + zc0*np.sin(th)
            yc = yc0
            zc = -xc0*np.sin(th) + zc0*np.cos(th)

            xc_init = xc0*np.cos(th) + z_init*np.sin(th)
            yc_init = yc0
            zc_init = -xc0*np.sin(th) + z_init*np.cos(th)

            xx, yy = np.meshgrid(x, y) # get 2D variables instead of 1D for rotated arrays

            # difference between r_init and r_obs in rotated sky-plane coordinates
            dr = max([xc_init,xc]) - min([xc_init,xc])

            # fit rotated ellipse
            p0 = (1, size, size, xc, yc, 0) # I, sx, sy, xc, yc, tilt
            popt, pcov = curve_fit(gauss2d, (xx, yy), xy.ravel(), p0=p0)
            fit = gauss2d((xx, yy), *popt)

            # define area and centroid position
            A = np.pi/4 * popt[1] * popt[2] / 3600 # armin^2
            xc = popt[3] # arcsec
            yc = popt[4] # arcsec

            # create 2D xy array for plotting
            xy_plot = xy.reshape(xx.shape)

            # plot image
            fig = plt.figure()
            ax = fig.add_subplot()
            ax.set_aspect('equal')

            cmap = plt.get_cmap('viridis')
            rect = patches.Rectangle((-xmax, -ymax), xmax*2, ymax*2, linewidth=0, facecolor=cmap(0.0), zorder=-1)
            ax.add_patch(rect)

            sun = np.linspace(0, 2*np.pi, 100)

            plt.contourf(x, y, xy_plot, levels=np.linspace(0,1,100), cmap=cmap, extent=[-xmax,xmax,-ymax,ymax], extend='both')
            #plt.contour(x, y, xy_plot, [0.5], linestyles = 'dashed', linewidths = 0.7, colors='white')
            plt.contour(x, y, xy_plot, [0.5,0.7,0.9], colors='w', linewidths=1, extent=[-xmax,xmax,-ymax,ymax])
            plt.scatter(xc, yc, marker='+', color='k', linewidths=0.7)
            plt.scatter(xc_init, 0, marker='x', s=20, color='k', linewidths=0.7)
            plt.plot(r_sun*np.cos(sun), r_sun*np.sin(sun), 'w')

            if A > 1e4:
                plt.text(0.05, 0.95, r'A  = {:.1e} arcmin$^2$'.format(A), color='w', ha='left', va='top', transform=ax.transAxes)
            else:
                plt.text(0.05, 0.95, r'A  = {:.0f} arcmin$^2$'.format(round(A,-1)), color='w', ha='left', va='top', transform=ax.transAxes)

            if xc_init > 1e4:
                plt.text(0.05, 0.90, r'$x_0$ = {:.1e} arcsec'.format(xc_init), color='w', ha='left', va='top', transform=ax.transAxes)
            else:
                plt.text(0.05, 0.90, r'$x_0$ = {:.0f} arcsec'.format(round(xc_init,-1)), color='w', ha='left', va='top', transform=ax.transAxes)

            if xc > 1e4:
                plt.text(0.05, 0.85, r'$x_c$ = {:.1e} arcsec'.format(xc), color='w', ha='left', va='top', transform=ax.transAxes)
            else:
                plt.text(0.05, 0.85, r'$x_c$ = {:.0f} arcsec'.format(round(xc,-1)), color='w', ha='left', va='top', transform=ax.transAxes)

            if dr > 1e4:
                plt.text(0.05, 0.80, r'$\Delta r$ = {:.1e} arcsec ({:.2f} R$_\odot$)'.format(dr, dr/r_sun), color='w', ha='left', va='top', transform=ax.transAxes)
            else:
                plt.text(0.05, 0.80, r'$\Delta r$ = {:.0f} arcsec ({:.2f} R$_\odot$)'.format(round(dr,-1), dr/r_sun), color='w', ha='left', va='top', transform=ax.transAxes)

            plt.text(0.95, 0.95, r'f = {:.2f} MHz'.format(f), color='w', ha='right', va='top', transform=ax.transAxes)
            plt.text(0.95, 0.90, r'$\alpha$ = {:.2f}'.format(anis), color='w', ha='right', va='top', transform=ax.transAxes)
            plt.text(0.95, 0.85, r'$\epsilon$ = {:.2f}'.format(eps), color='w', ha='right', va='top', transform=ax.transAxes)


            plt.xlabel('X [arcsec]')
            plt.ylabel('Y [arcsec]')
            plt.xlim(-xmax,xmax)
            plt.ylim(-ymax,ymax)
            ax.tick_params(axis='both', which='both', direction='in', color='white', top=True, right=True)
            ax.tick_params(axis='both', which='major', length=5)
            ax.tick_params(axis='both', which='minor', length=3)
            minor_locator = AutoMinorLocator()
            ax.xaxis.set_minor_locator(minor_locator)
            ax.yaxis.set_minor_locator(minor_locator)
            # Get the current tick labels
            xticklabels = ax.get_xticklabels()
            yticklabels = ax.get_yticklabels()
            # Set the first and last tick labels to an empty string
            xticklabels[0] = ''
            xticklabels[-1] = ''
            yticklabels[0] = ''
            yticklabels[-1] = ''
            # Set the tick labels
            ax.set_xticklabels(xticklabels)
            ax.set_yticklabels(yticklabels)

            plt.tight_layout()
            plt.show()

        else:

            print("Frequency must be between 0.1-340 MHz.")
            print("Epsilon can only be 0.5, 0.71, 1.0, 1.41, or 2.0.")
            print("Anisotropy can only be 0.19, 0.25, 0.33, or 0.42.")

    else:

        print('Data file "new_all_results_F_final.sav" missing. Place in same directory as scattering_image_funcs.py.')
