import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from flask import Flask, render_template
from flask import request

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from scipy.io import readsav
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

app = Flask(__name__)

@app.route('/')
def run():
    return render_template('index.html')

@app.route('/get_plot', methods = ['GET', 'POST'])
def get_plot():
    if request.method == "POST":
        
        #================================
        # User input
        f = float(request.form['f'])
        eps0 = float(request.form['eps0'])
        anis0 = float(request.form['anis0'])
        xmax = float(request.form['xmax'])
        th = np.deg2rad(float(request.form['th']))
        #================================

        # constants
        r_sun = 950 # solar radius, arcsec

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

        # Load data
        dataDir = './static/'
        file = 'sim_data.sav'
        data = readsav(dataDir + file, python_dict=True, verbose=True)

        # parse ds data
        r_init = np.array(data['sims'][0][1])
        fmhz = np.array(data['sims'][0][2])
        r_shift = np.array(data['sims'][0][3])
        s_fwhm = np.array(data['sims'][0][4])
        eps = np.array(data['sims'][0][5])
        asym = np.array(data['sims'][0][6])
        anis = np.array(data['sims'][0][7])
        f_ratio = np.array(data['sims'][0][8])
        t_1e_fit = np.array(data['sims'][0][9])
        t_1e_after = np.array(data['sims'][0][10])
        t_peak = np.array(data['sims'][0][11])

        # retrieve data for eps, alpha
        idx = np.where((eps == eps0) & (anis == anis0))
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

        ymax = xmax
        x0 = np.linspace(-xmax, xmax, 240)
        y0 = x0
        z0 = x0*0 + r_shift*r_sun
        xc0 = 0
        yc0 = 0
        zc0 = r_shift*r_sun
        z_init = r_init*r_sun

        # Generate 2D Gaussian
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

        plt.text(0.05, 0.95, r'A  = {:.1f} arcmin$^2$'.format(A), color='w', ha='left', va='top', transform=ax.transAxes)
        plt.text(0.05, 0.90, r'$x_0$ = {:.1f} arcsec'.format(xc_init), color='w', ha='left', va='top', transform=ax.transAxes)
        plt.text(0.05, 0.85, r'$x_c$ = {:.1f} arcsec'.format(xc), color='w', ha='left', va='top', transform=ax.transAxes)
        plt.text(0.05, 0.80, r'$\Delta r$ = {:.1f} arcsec ({:.2f} R$_\odot$)'.format(dr, dr/r_sun), color='w', ha='left', va='top', transform=ax.transAxes)

        plt.text(0.95, 0.95, r'f = {:.2f} MHz'.format(f), color='w', ha='right', va='top', transform=ax.transAxes)
        plt.text(0.95, 0.90, r'$\alpha$ = {:.2f}'.format(anis0), color='w', ha='right', va='top', transform=ax.transAxes)
        plt.text(0.95, 0.85, r'$\epsilon$ = {:.2f}'.format(eps0), color='w', ha='right', va='top', transform=ax.transAxes)

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
        xticklabels = ax.get_xticklabels()
        yticklabels = ax.get_yticklabels()
        xticklabels[0] = ''
        xticklabels[-1] = ''
        yticklabels[0] = ''
        yticklabels[-1] = ''
        ax.set_xticklabels(xticklabels)
        ax.set_yticklabels(yticklabels)

        plt.tight_layout()
        plt.savefig('static/my_plot.png', dpi=300)
        plt.savefig('static/my_plot.eps')
        return render_template('index.html', plot_url = "static/my_plot.png", data = data)
    
    else:
    
        return render_template('index.html')
    
app.secret_key = 'some secret that you will never guess'

if __name__ == "__main__":
    app.run('127.0.0.1', 5000, debug = True)