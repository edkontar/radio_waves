pro rays_swind, qeps2_scaling = qeps2_scaling, anis = anis, r_init = R_init, asym = asym, f_ratio = f_ratio
  ; +
  ; Ray tracing and scattering based on the following paper:
  ; http://adsabs.harvard.edu/abs/2019arXiv190900340K
  ; the routine works in Sun centered cartesian coordinates
  ; z is along the line of sight
  ;
  ; EXAMPLE: Calling sequence, e.g. for LOFAR 30 MHz observations:
  ; IDL> rays_swind, qeps2_scaling=0.1, anis=0.25, r_init=1.75, asym=1., f_ratio=2.0
  ; 
  ; PARAMETERS:
  ; qeps2_scaling - scaling of qesp2 function (used 1./4,1/2, 1, 2, 4); default = 1 like in
  ;                  ApJ Paper https://doi.org/10.3847/1538-4357/acf6c1
  ; 
  ; anis   - wavenumber anisotropy of density turbulence q_par/q_perp   ; default =0.25 
  ;          anis->+infty, quasi-paraller density fluctuations
  ;          anis=1 isotropic density fluctuations 

  ; r_init - initial location of the source (to check local plasma frequency IDL> omega_pe(r_init)/2/!Pi/1e6 [MHz])
  ; R_init=1.75D0 ;(around 32 MHz)
  ; R_init=5.D0 ;(around 2 MHz)
  ; R_init=2.21D0 ;(around 15 MHz)
  ; R_init=1.1372D0 ;(around 329 MHz)
  ; R_init=1.3D0 ;(around 100 MHz)
  ; R_init=57.D0 ;(around 0.118 MHz)
  ; R_init=10.D0 ;(around 0.765 MHz)
  ; R_init=20.D0 ;(around 0.338 MHz)
  ;   
  ; f_ration - ratio of emitted frequency to local plasma frequency at the source
  ;            i.e. omega/omega_pe(r_init)
  ; 
  ; asym  - asymetry parameter density fluctuations assymetry sunward/anisunward 
  ; asym=1 - symmetric;
  ; 0<asym<1 => more outward density fluctuation (inward radio stronger scattered)
  ; 1<asym<2 => more inward density fluctuations (outward radio stronger scattered)
  ; <delta n^2>_in= asym *<delta n^2>
  ; <delta n^2>_out = (2- asym) *<delta n^2>
  ; z_anis=asym*(k_along_r GT 0.) + (2.-asym)*(k_along_r LE 0.)
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
  ; updated: October 2023 by eduard@glasgow
  ; include qeps2 profile

  N = 9999l
  ; photon number
  R_S = 6.96e10
  c = 2.998e10
  c_r = c / R_S
  ; speed of light in r_sun units
  c_s = 2e7 / R_S
  ; sound speed in R_sun units
  au = 215.
  ; astronomical unit

  ; if input parameters undefined - we set default parameters
  if (n_elements(qeps2_scaling) ne 1.0) then qeps2_scaling = 1.0
  if (n_elements(R_init) ne 1) then R_init = 20.0 ; r=20 fpe=0.338MHz
  if (n_elements(anis) ne 1) then anis = 0.25
  if (n_elements(f_ratio) ne 1) then f_ratio = 2.0
  if (n_elements(asym) ne 1) then asym = 1. 
 
  theta_0 = 30.d0
  r = findgen(N) * 0.d0 + R_init
  rtheta = (theta_0 + findgen(N) / N * 0.) * !pi / 180.d0
  rfi = (0. + findgen(N) / N * 0.) * !pi / 180.d0
  ; r=r+randomu(seed,N)*1e-1

  Te = 86.0 * 1.5 ; eV temperature in eV
  nu_e = 2.91e-6 * density_r(R_init) * 20. / Te ^ 1.5
  rz = r * cos(rtheta)
  rx = r * sin(rtheta) * cos(rfi)
  ry = r * sin(rtheta) * sin(rfi)
  Iindex = indgen(N)

  ; large (10 R_sun) line-like source
  ; rx=rx+(randomu(seed,N)-0.5)*10.

  r0 = r
  rx0 = rx
  ry0 = ry
  rz0 = rz
  ; initial locations for future use

  f_pe0 = omega_pe(r) / 2. / !pi
  ; initial plasma frequency
  omega_pe0 = omega_pe(R_init)
  f_start = omega_pe0 / 1e6 / 2 / !pi
  print, format = '("INIT: Source located at R_source=",(F8.2)," R_sun, Angle = ",(F6.1)," degrees")', R_init, theta_0
  print, format = '("f_pe=",(F8.3)," MHz for emission frequency=",(F8.3)," MHz")', f_start, f_start * f_ratio
  print, format = '("KEY PARMETERS: qeps2_scaling=",(F5.2),", Anisotropy= ",(F5.2),", Asymmetry= ",(F5.2))', qeps2_scaling, anis, asym

  omega = f_ratio * omega_pe(r0)
  ; emission at 1.1 of plasma frequency in units of initial omega_pe
  ; 2 for harmonic
  omega0 = min(omega)
  omega = omega0

  kc0 = sqrt(omega ^ 2 - omega_pe(r) ^ 2)
  kc = kc0
  tau = fltarr(N) * 0.0

  seed = 1001l
  mu = (randomu(seed, N) - 0.5) * 2.d
  ; isotropic
  ;
  ; mu =randomu(seed,N)
  ; positive mu - only to the observer
  kfi = randomu(seed, N) * 2. * !pi
  kz = kc * mu
  kx = kc * sqrt(1. - mu * mu) * cos(kfi)
  ky = kc * sqrt(1. - mu * mu) * sin(kfi)

  ; stop

  linecolors
  window, 1, xsize = 900, ysize = 800
  !p.multi = [0, 2, 2]

  plot, rz0, rx0, psym = 2, xrange = ([-3 * R_init, 3 * R_init]), yrange = ([-3 * R_init, 3 * R_init]), $
    xtitle = 'Z [R!D!9n!3!N]', ytitle = 'X [R!D!9n!3!N]'
  pp = findgen(41) / 40 * !pi * 2
  oplot, cos(pp), sin(pp), thick = 2, color = 5
  oplot, [0, 0], [-3 * R_init, 3 * R_init], line = 1
  oplot, [-3 * R_init, 3 * R_init], [0, 0], line = 1

  plot, rx0, ry0, psym = 2, xrange = [-3 * R_init, 3 * R_init], yrange = [-3 * R_init, 3 * R_init], $
    xtitle = 'X [R!D!9n!3!N]', ytitle = 'Y [R!D!9n!3!N]'
  pp = findgen(41) / 40 * !pi * 2
  oplot, cos(pp), sin(pp), thick = 2, color = 5
  oplot, [0, 0], [-3 * R_init, 3 * R_init], line = 1
  oplot, [-3 * R_init, 3 * R_init], [0, 0], line = 1

  ; density fluctuations
  nu_s0 = nu_scat(R_init, min(omega), qeps2_scaling)

  ; rint=R_init*(1.+10.^(findgen(999)/999.*5-3.5))
  rint = R_init * (1. + findgen(3999) / 49 + 1e-3)

  N_rint = n_elements(rint)
  dri = shift(rint, -1) - rint
  dri[N_rint - 1] = 0.
  tau_scat = fltarr(N_rint)
  for j = 0, n_elements(rint) - 3 do tau_scat[j] = int_tabulated(rint[j : N_rint - 1], $
    nu_scat(rint[j : N_rint - 1], omega0, qeps2_scaling) / c_r / sqrt(1. - omega_pe(rint[j : N_rint - 1]) ^ 2 / omega0 ^ 2))
  ; total(nu_scat(rint[j:N_rint-1],min(omega),qeps2_scaling)$
  ; *dri[j:N_rint-1]/c_r/sqrt(1.-omega_pe(rint[j:N_rint-1])^2/min(omega)^2))

  tau_scat[n_elements(rint) - 1] = tau_scat[n_elements(rint) - 3]
  tau_scat[n_elements(rint) - 2] = tau_scat[n_elements(rint) - 3]

  min_tau10 = min(abs(tau_scat - 1. / 10), mloc10)
  min_tau1 = min(abs(tau_scat - 1. / 1), mloc1)
  r_scat = rint[mloc10]
  r_tau1 = rint[mloc1]
  print, 'Scatterint radius =', r_scat
  ; calculating the scattering rate
  ; scattering rate is by definition D_mu_mu=nu_s*(1-mu^2)/2 for isotropic scattering
  ; nu_s is already multiplied by c/omega_0

  dT_save = 2e7 / max(omega)
  ; dt0=0.05/nu_s0
  Exp_size = 1.25 * 30. / f_start
  dt0 = 0.01 * Exp_size / c_r
  t = 0.d0
  tplot = 0.
  ; stop

  print, 'scattering time=', 1. / nu_s0, ' sec', '  DT save =', dT_save

  ; the block below determines the maximum distance
  ; omega_pe/omega=0.1
  ; min_stop =min(abs(omega_pe(Rint)/min(omega)-0.01),istop)
  ; R_stop=rint[istop]
  min_tau = min(abs(tau_scat - (0.1 * 30 * !pi / 180 / f_start) ^ 2), mloc)
  ; angular spread less than (0.1*32/f_start)^2 radians^2
  r_stop = max([rint[mloc], r_scat, R_init * 2, 215])
  r_stop = min([r_stop, 215.]) ; Not larger than 1 AU
  dT_save = (r_stop - R_init) / c_r / 10.

  plot_oo, omega_pe(rint) / 2 / !pi / 1e6, 1. / nu_scat(rint, min(omega), qeps2_scaling), $
    xrange = [omega_pe(R_init), omega_pe(r_stop)] / 2. / !pi / 1e6, xtitle = 'Plasma frequency [MHz]', ytitle = 'Scattering time [seconds]'
  oplot, omega_pe(rint) / 2 / !pi / 1e6, (rint - R_init) / c_r, line = 2

  plot_oo, rint, tau_scat, xrange = [1, r_stop], $
    xtitle = 'Distance [R!D!9n!3!N]', yrange = [1e-1, 10], $
    ytitle = 'Scattering depth \tau'
  oplot, [r_stop, r_stop], [1e-1, 10], line = 1
  oplot, [R_init, R_init], [1e-1, 10], line = 1
  oplot, [r_tau1, r_tau1], [1e-1, 10], line = 2
  !p.multi = 0
  print, 'Stopping distance for calculations=', r_stop, ' R_sun'

  ; manually changing
  ; r_stop=3.
  ; stop
  ; while ((T LE tmax) AND (N_elements(r) GE 10)) DO begin

  while (n_elements(r) ge 2) do begin
    r = sqrt(rx * rx + ry * ry + rz * rz)
    kc = sqrt(kx * kx + ky * ky + kz * kz)
    omega = sqrt(omega_pe(sqrt(rx * rx + ry * ry + rz * rz)) ^ 2 + kc * kc)

    nu_s = nu_scat(r, omega, qeps2_scaling)
    nu_s = nu_s0 * (nu_s gt nu_s0) + nu_s * (nu_s le nu_s0)
    g0 = sqrt(nu_s * kc ^ 2)
    ; print, 1./max(nu_s)
    dt_ref = min(abs(kc / (domega_dr(r) * c_r) / 20.))
    dt_dr = min(r / (c_r / omega0 * kc)) / 20.
    ; mean scattering time in second
    dt = min([.3d0 / max(nu_s), dt0, dt_ref, dt_dr])
    ; dt=dt0/5

    ; distance measured in solar radii
    drx_dt = c_r / omega * kc * kx / kc
    dry_dt = c_r / omega * kc * ky / kc
    drz_dt = c_r / omega * kc * kz / kc

    rand = randomn(seed, 3, n_elements(rz))
    Wx = reform(rand[2, *]) * sqrt(dt) ; *sqrt(3./(2.+anis^2))
    Wy = reform(rand[1, *]) * sqrt(dt) ; *sqrt(3./(2.+anis^2))
    Wz = reform(rand[0, *]) * sqrt(dt) ; *sqrt(3.*anis^2/(2+anis^2))

    ; calculates positions of the photons in spherical coordinates

    fi0 = atan(ry, rx)
    sintheta0 = sqrt(1. - rz ^ 2 / r ^ 2)
    costheta0 = rz / r

    sph_E = cv_coord(from_rect = transpose([[rz], [rx], [ry]]), /to_sphere, /double)
    ; costheta1=cos(!Pi/2.D0-reform(sp_[1, *]))
    ; spherical coordinates with x to the Earth

    parker_Bxyz_Earth, reform(sph_E[2, *]), !pi / 2.d0 - reform(sph_E[1, *]), reform(sph_E[0, *]), bxE, byE, bzE

    BB = sqrt(bxE ^ 2 + byE ^ 2 + bzE ^ 2)
    fi = atan(bzE, byE)
    sintheta = sqrt(1. - bxE ^ 2 / BB ^ 2)
    costheta = bxE / BB

    ; stop
    ; ***********************************************
    ; Rotate k to r-vector aligned coordinate system

    kcx = -kx * sin(fi) + ky * cos(fi)
    kcy = -kx * costheta * cos(fi) - ky * costheta * sin(fi) + kz * sintheta
    kcz = kx * sintheta * cos(fi) + ky * sintheta * sin(fi) + kz * costheta

    kw = Wx * kcx + Wy * kcy + Wz * kcz * anis
    Akc = sqrt(kcx * kcx + kcy * kcy + kcz * kcz * anis ^ 2)
    z_asym = (asym * (kcz gt 0.) + (2. - asym) * (kcz le 0.)) * (kc / Akc) ^ 2

    A_perp = nu_s * z_asym * kc / Akc ^ 3 * (-(1 + anis ^ 2) * Akc ^ 2 + 3 * anis ^ 2 * (anis ^ 2 - 1) * kcz ^ 2) * anis
    A_par = nu_s * z_asym * kc / Akc ^ 3 * ((-3 * anis ^ 4 + anis ^ 2) * Akc ^ 2 + 3 * anis ^ 4 * (anis ^ 2 - 1) * kcz ^ 2) * anis
    A_g0 = g0 * sqrt(z_asym * anis)

    kcx = kcx + A_perp * kcx * dt + A_g0 * (Wx - kcx * kw / Akc ^ 2)
    kcy = kcy + A_perp * kcy * dt + A_g0 * (Wy - kcy * kw / Akc ^ 2)
    kcz = kcz + A_par * kcz * dt + A_g0 * (Wz - kcz * kw * anis / Akc ^ 2) * anis

    ; print, 'dk^2 =',minmax((kcx1)*kcx1+(kcy1)*kcy1+(kcz1)*kcz1), minmax(kcx*kcx+kcy*kcy+kcz*kcz)
    ; stop
    ; kc_old2new=kc/sqrt(kx*kx+ky*ky+kz*kz)
    ; kx=kx*kc_old2new
    ; ky=ky*kc_old2new
    ; kz=kz*kc_old2new
    ; kc=sqrt(kx*kx+ky*ky+kz*kz)
    ;
    ;
    ; ax=c_s^2*omega^2/kc^2 & ay=0. & az=0.
    ; kcx=kcx - nu_s*(1.+ax)*kcx*dt + g0*(wx-kcx*kw/kc^2*(1.-sqrt(2*ax)));*sqrt(abs(ay*az))
    ; kcy=kcy - nu_s*(1.+ay)*kcy*dt + g0*(wy-kcy*kw/kc^2*(1.-sqrt(2*ay)));*sqrt(abs(ax*az))
    ; kcz=kcz - nu_s*(1.+az)*kcz*dt + g0*(wz-kcz*kw/kc^2*(1.-sqrt(2*az)));*sqrt(abs(ax*ay))
    ; kz=kz+ g0*abs(Wz)
    ; kx=kx  + g0*(-wx*sin(fi)+ (wz*sintheta-wy*costheta)*cos(fi))
    ; ky=ky  + g0*( wx*cos(fi)+ (wz*sintheta-wy*costheta)*sin(fi))
    ; kz=kz  + g0*(wz*costheta + wy*sintheta)
    ; Thejappa approach

    kx = -kcx * sin(fi) - kcy * costheta * cos(fi) + kcz * sintheta * cos(fi)
    ky = kcx * cos(fi) - kcy * costheta * sin(fi) + kcz * sintheta * sin(fi)
    kz = kcy * sintheta + kcz * costheta

    ; after scattering rotate back from r-aligned to XYZ coordinate system
    ; ***********************************************

    kc_norm = sqrt(kx * kx + ky * ky + kz * kz)
    kx = kx * kc / kc_norm
    ky = ky * kc / kc_norm
    kz = kz * kc / kc_norm

    dk_dt = (omega_pe(r) / omega) * domega_dr(r) * c_r

    kx = kx - dk_dt * dt * (rx / r)
    ky = ky - dk_dt * dt * (ry / r)
    kz = kz - dk_dt * dt * (rz / r)

    if (total(finite(kx, /nan)) ge 1.0) then stop
    ; distance measured in solar radii
    rx = rx + drx_dt * dt
    ry = ry + dry_dt * dt
    rz = rz + drz_dt * dt
    ; stop

    ; To conserve the frequency
    kc_new_old = sqrt(omega ^ 2 - omega_pe(sqrt(rx * rx + ry * ry + rz * rz)) ^ 2) / sqrt(kx * kx + ky * ky + kz * kz)
    kx = kx * kc_new_old
    ky = ky * kc_new_old
    kz = kz * kc_new_old
    ;
    ; kc_norm=sqrt(kx*kx+ky*ky+kz*kz)
    nu_e = 2.91e-6 * density_r(r) * 20. / Te ^ 1.5 * omega_pe(r) ^ 2 / omega ^ 2
    ; nu_col=!PI*4.8d-10^4*density_r(r)*20./sqrt(9.1d-28*(11605.*Te*1.38d-16)^3.)*omega_pe(r)^2/omega^2
    ; stop
    tau = tau + nu_e * dt
    ; optical depth of the photons
    ; stop
    ; frac=(1.5e7/3e10)^2
    ; kc=kc - nu_s*kc*dt*frac + sqrt(g0*frac)*RANDOMN(seed,N_elements(rz))*sqrt(dt)

    ;

    t = t + dt

    if (max(tau) ge 4.61) then begin ; exp(-4.61)=0.01
      ; stop
      ; retaining photons that are not absorbed
      not_absorbed = where((tau - min(tau)) le 4.61) ; exp(-4.61)=0.01
      rx = rx[not_absorbed]
      ry = ry[not_absorbed]
      rz = rz[not_absorbed]
      kx = kx[not_absorbed]
      ky = ky[not_absorbed]
      kz = kz[not_absorbed]
      tau = tau[not_absorbed]
      ; stop
    endif

    if n_elements(where(r ge r_stop, /null)) ge 1 then begin
      if n_elements(rr) eq 0 then begin
        ; tt=t & rr=0. & kk=0. & Rxx=0. &Ryy=0. &Rzz=0. &Kxx=0. & Kyy=0. & Kzz=0.

        tt = fltarr(n_elements(where(r ge r_stop))) + t
        kxx = kx[where(r ge r_stop, /null)]
        kx = kx[where(r lt r_stop)]
        kyy = ky[where(r ge r_stop, /null)]
        ky = ky[where(r lt r_stop)]
        kzz = kz[where(r ge r_stop, /null)]
        kz = kz[where(r lt r_stop)]

        rxx = rx[where(r ge r_stop, /null)]
        rx = rx[where(r lt r_stop)]
        ryy = ry[where(r ge r_stop, /null)]
        ry = ry[where(r lt r_stop)]
        rzz = rz[where(r ge r_stop, /null)]
        rz = rz[where(r lt r_stop)]
        rxx0 = rx0[where(r ge r_stop, /null)]
        rx0 = rx0[where(r lt r_stop)]
        ryy0 = ry0[where(r ge r_stop, /null)]
        ry0 = ry0[where(r lt r_stop)]
        rzz0 = rz0[where(r ge r_stop, /null)]
        rz0 = rz0[where(r lt r_stop)]

        tau0 = tau[where(r ge r_stop, /null)]
        tau = tau[where(r lt r_stop)]
        kk = kc[where(r ge r_stop, /null)]
        kc = kc[where(r lt r_stop)]
        rr = r[where(r ge r_stop, /null)]
        r = r[where(r lt r_stop)]

        ; Important: rr or rz should be last line!
      endif else begin
        tt = [tt, fltarr(n_elements(where(r ge r_stop))) + t]

        kxx = [kxx, kx[where(r ge r_stop, /null)]]
        kx = kx[where(r lt r_stop)]
        kyy = [kyy, ky[where(r ge r_stop, /null)]]
        ky = ky[where(r lt r_stop)]
        kzz = [kzz, kz[where(r ge r_stop, /null)]]
        kz = kz[where(r lt r_stop)]

        rxx = [rxx, rx[where(r ge r_stop, /null)]]
        rx = rx[where(r lt r_stop)]
        ryy = [ryy, ry[where(r ge r_stop, /null)]]
        ry = ry[where(r lt r_stop)]
        rzz = [rzz, rz[where(r ge r_stop, /null)]]
        rz = rz[where(r lt r_stop)]
        rxx0 = [rxx0, rx0[where(r ge r_stop, /null)]]
        rx0 = rx0[where(r lt r_stop)]
        ryy0 = [ryy0, ry0[where(r ge r_stop, /null)]]
        ry0 = ry0[where(r lt r_stop)]
        rzz0 = [rzz0, rz0[where(r ge r_stop, /null)]]
        rz0 = rz0[where(r lt r_stop)]

        tau0 = [tau0, tau[where(r ge r_stop, /null)]]
        tau = tau[where(r lt r_stop)]
        kk = [kk, kc[where(r ge r_stop, /null)]]
        kc = kc[where(r lt r_stop)]
        rr = [rr, r[where(r ge r_stop, /null)]]
        r = r[where(r lt r_stop)]

        ; important: rr or rz should be last line!
      endelse
    endif

    if (t ge tplot + dT_save) then begin
      tplot = tplot + dT_save
      ; window,1,xsize=600,ysize=600
      ; plot,rz0,rx0,psym=2, xrange=[0,8],yrange=[-4,4], xtitle='X [R!D!9n!3!N]',ytitle='Y [R!D!9n!3!N]'
      ; pp = findgen(41)/40*!pi*2
      ; oplot,cos(pp), sin(pp),thick=2
      ; dT_save=dt
      wset, 1
      window, 4, xsize = 600, ysize = 900
      !p.multi = [0, 1, 2]
      obs_i0 = where(kz / kc gt 0.9)

      dir_beam = average(atan(kx, kz))

      plot, [rx0], [ry0], psym = 2, xrange = [-r_stop, r_stop], yrange = [-r_stop, r_stop], xtitle = 'X [R!D!9n!3!N]', ytitle = 'Y [R!D!9n!3!N]'
      pp = findgen(41) / 40 * !pi * 2
      oplot, cos(pp), sin(pp), thick = 2, color = 5
      oplot, [rx], [ry], psym = 3, color = 7
      oplot, r_stop * cos(pp), r_stop * sin(pp), line = 1
      oplot, [rx[obs_i0]], [ry[obs_i0]], psym = 3, color = 3
      oplot, [0, 0], [-r_stop, r_stop], line = 1
      oplot, [-r_stop, r_stop], [0, 0], line = 1

      ; wset, 0
      plot, [rz0], [rx0], psym = 2, xrange = [-r_stop, r_stop], yrange = [-1., 1.] * r_stop, xtitle = 'Z [R!D!9n!3!N]', ytitle = 'X [R!D!9n!3!N]', $
        title = 'Emission angle =' + string(format = '(F6.1)', dir_beam * 180 / !pi) + '  t=' + string(format = '(F7.2)', t)
      pp = findgen(41) / 40 * !pi * 2
      oplot, cos(pp), sin(pp), thick = 2, color = 5
      oplot, r_stop * cos(pp), r_stop * sin(pp), line = 1
      oplot, [0, 0], [r_stop, -r_stop], line = 1
      oplot, [0, r_stop], [0, 0], line = 1
      ; rp =0.1*(r_stop-r_init)
      rp = 0.333 * (r_stop - R_init)
      rzc = sqrt(r_stop ^ 2 - rp ^ 2)
      oplot, [rzc, rzc], [-rp, rp], color = 9
      oplot, [rz], [rx], psym = 3, color = 7
      oplot, [rz[obs_i0]], [rx[obs_i0]], psym = 3, color = 3

      oplot, [0, cos(dir_beam)] * r_stop, [0, sin(dir_beam)] * r_stop, line = 2, color = 4
      ; oploting emission direction

      rp = findgen(100) / 100 * r_stop + 1
      parker_Bxyz_Earth, rp, !pi / 2, !pi / 2., bxE, byE, bzE
      oplot, rp, atan(byE, bxE), color = 2, /polar, line = 2
      parker_Bxyz_Earth, rp, !pi / 2, !pi / 4., bxE, byE, bzE
      oplot, rp, atan(byE, bxE), color = 2, /polar, line = 2
      parker_Bxyz_Earth, rp, !pi / 2, 0., bxE, byE, bzE
      oplot, rp, atan(byE, bxE), color = 2, /polar, line = 2
      parker_Bxyz_Earth, rp, !pi / 2, -!pi / 4., bxE, byE, bzE
      oplot, rp, atan(byE, bxE), color = 2, /polar, line = 2
      ; stop

      tfname = 't' + string(format = '(F07.2)', t) + 'fpe' + string(format = '(I7.7)', omega_pe0 / 1e3 / 2 / !pi) $
        + 'kHz_FE' + string(format = '(I7.7)', omega0 / 1e3 / 2 / !pi) $
        + 'anis' + string(format = '(F4.2)', anis) + 'qeps2_scaling' + string(format = '(F5.3)', qeps2_scaling) $
        + 'asym' + string(format = '(F4.2)', asym) + 'fr' + string(format = '(F4.2)', f_ratio) + '.png'

      write_png, tfname, tvrd(/true)
      wait, 0.1
      !p.multi = [0, 1, 2]

      print, 't = ', t, '<r> =', average(r), ' N_phot=', n_elements(tau), ' <weight>=', average(exp(-tau))
      print, 'F =', average(sqrt(omega_pe(r) ^ 2 + kc ^ 2) / omega), minmax(sqrt(omega_pe(r) ^ 2 + kc ^ 2) / omega)
      ; stop
      r_tab = R_init + (r_stop - R_init) * findgen(900) / (900. - 1)
      d_delay = int_tabulated(r_tab, 1. / sqrt(1. - omega_pe(r_tab) ^ 2 / omega0 ^ 2)) / c_r
      tau_exp = int_tabulated(r_tab, 2.91d-6 * density_r(r_tab) * 20. / Te ^ 1.5 * omega_pe(r_tab) ^ 2 / omega0 ^ 2 / sqrt(1. - omega_pe(r_tab) ^ 2 / omega0 ^ 2)) / c_r

      ; stop

      wait, 0.01
      ; stop
    endif
  endwhile

  r_free = sqrt((rxx - rxx0) ^ 2 + (ryy - ryy0) ^ 2 + (rzz - rzz0) ^ 2)
  rr = sqrt(rxx ^ 2 + ryy ^ 2 + rzz ^ 2)
  kk = sqrt(kxx ^ 2 + kyy ^ 2 + kzz ^ 2)
  t_free = r_free / c_r
  print, 'Delta t= ', d_delay - min(t_free)

  ; plot_rays_new, file=fname
  ; stop
  ; direction cosine
  cosx = kxx / kk
  cosy = kyy / kk
  cosz = kzz / kk
  fi = atan(cosy, cosx)

  ; cosx0=(rxx-rxx0)/r_free & cosy0=(ryy-ryy0)/r_free & cosz0=(ryy-ryy0)/r_free
  ; xa_free=atan((rxx-rxx0)/(rzz-rzz0))
  ; ya_free=atan((ryy-ryy0)/(rzz-rzz0))

  ; collecting photons at infinity into angle (far-away source)
  ; correct method for a far away source
  obs_i = where((kzz / kk ge 0.9) and (kzz / kk le 1.0))
  ww = exp(-tau0[obs_i])
  x0 = r_free * kxx / kk
  y0 = r_free * kyy / kk
  x_im = rxx[obs_i] - x0[obs_i]
  y_im = ryy[obs_i] - y0[obs_i]

  ; x0=sqrt((rzz-Rz0)^2+(rxx-rx0)^2+(ryy-ry0)^2)*cosx & y0=sqrt((rzz-rz0)^2+(rxx-rx0)^2+(ryy-ry0)^2)*cosy
  ; x0=(r_stop-R_init)*tan(acos(cosr_stop)) & y0=(r_stop-R_init)*tan(acos(cosr_stop))
  ; x0=(rzz-R_init)*sqrt((rzz-R_init)^2+rxx^2+ryy^2)/*cosx & y0=max(rzz)*cosy

  ; x0=(rzz-Rzz0)*kxx/kzz & y0=(rzz-Rzz0)*kyy/kzz
  ; x0=kxx/kk*c_r*tt & y0=kyy/kk*c_r*tt

  ; time integrated sourse FWHM sizes centroids and shifts in R_sun units
  xc = average(x_im * ww) / average(ww)
  yc = average(y_im * ww) / average(ww)
  sx = sqrt(average(x_im ^ 2 * ww - xc ^ 2 * ww) / average(ww)) * 2.355
  sy = sqrt(average(y_im ^ 2 * ww - yc ^ 2 * ww) / average(ww)) * 2.355

  ; collecting photons at 1AU

  z_shift = (au - rzz) / cosz
  x_au = rxx + cosx * z_shift
  y_au = ryy + cosy * z_shift
  z_au = rzz + cosz * z_shift
  rp_1au = sqrt(x_au ^ 2 + y_au ^ 2)
  obs_im = where(rp_1au le 0.333 * (r_stop - R_init) * au / r_stop)
  xa_im = atan(kxx[obs_im] / kzz[obs_im]) - $
    atan((x_au[obs_im] - rxx0[obs_im]) / (z_au[obs_im] - rzz0[obs_im])) + atan(-rxx0[obs_im] / (z_au[obs_im] - rzz0[obs_im]))
  ya_im = atan(kyy[obs_im] / kzz[obs_im]) - $
    atan((y_au[obs_im] - ryy0[obs_im]) / (z_au[obs_im] - rzz0[obs_im])) + atan(-ryy0[obs_im] / (z_au[obs_im] - rzz0[obs_im]))
  wa = exp(-tau0[obs_im])
  ; angle imaging for closely located sources
  xac = average(xa_im * wa) / average(wa)
  yac = average(ya_im * wa) / average(wa)
  sax = sqrt(average(xa_im ^ 2 * wa - xac ^ 2 * wa) / average(wa)) * 2.355
  say = sqrt(average(ya_im ^ 2 * wa - yac ^ 2 * wa) / average(wa)) * 2.355
  print, 'X/Y angle sizes [degrees]', sax * 180. / !pi, say * 180. / !pi

  Err_xc = sx / sqrt(n_elements(x_im) + 1e-8) / 2.355
  Err_yc = sy / sqrt(n_elements(y_im) + 1e-8) / 2.355

  Err_sx = sx * sqrt(2.) / sqrt(n_elements(x_im) + 1e-8)
  Err_sy = sy * sqrt(2.) / sqrt(n_elements(y_im) + 1e-8)

  ; without absorption
  xc_noa = average(x_im)
  yc_noa = average(y_im)
  sx_noa = sqrt(average(x_im ^ 2 - xc_noa ^ 2)) * 2.355
  sy_noa = sqrt(average(y_im ^ 2 - yc_noa ^ 2)) * 2.355

  xlimits = [-1.1, 1.1] * max([xc * 2, sx * 2, 1.])
  ylimits = [-1.1, 1.1] * max([yc * 2, sy * 2, 1.])

  window, 7, xsize = 1300, ysize = 800
  !p.multi = [0, 3, 2]
  !p.charsize = 1.4

  ; angular distribution for all photons at r>r_stop
  ; dtbin=(average(tt-t_free)-min(tt-t_free))/10.
  ; t_max=average(tt-t_free)*2
  ; pdf_mu=RAYS_HISTOGRAM_WEIGHT(tt-t_free,WEIGHT=exp(-tau0)*kzz/kk,locations=tbins,bin=dtbin,min=0.,max=t_max)/$
  ; (RAYS_HISTOGRAM_WEIGHT(tt-t_free,WEIGHT=exp(-tau0),bin=dtbin,min=0.,max=t_max)+1e-10)
  ; pdf_mu1=RAYS_HISTOGRAM_WEIGHT(tt-t_free,WEIGHT=kzz/kk,bin=dtbin,min=0.,max=t_max)/$
  ; (RAYS_HISTOGRAM_WEIGHT(tt-t_free,WEIGHT=1.0,bin=dtbin,min=0.,max=t_max)+1e-10)

  ; angular distribution of ALL photons
  cosA = [(kzz / kk - 1), -(kzz / kk - 1)]
  w2 = [exp(-tau0), exp(-tau0)]
  muc = average(w2 * cosA) / average(w2)
  smu = sqrt(average((cosA - muc) ^ 2 * w2) / average(w2))
  thetac = average(w2 * acos(cosA)) / average(w2)
  s_theta = sqrt(average((acos(cosA) - thetac) ^ 2 * w2) / average(w2)) * 180. / !pi
  print, 'direction =', thetac * 180. / !pi, ' degrees and FWHM angle =', s_theta, 'degrees'

  mubin = 0.01
  pdf_mu = rays_histogram_weight((kzz / kk), weight = exp(-tau0), locations = mubins, bin = mubin, min = 0, max = 1)
  mu_HalfW = mubins[min(where(pdf_mu / max(pdf_mu) ge 0.5))] + mubin / 2.
  mubin = 0.05
  pdf_mu = rays_histogram_weight((kzz / kk), weight = exp(-tau0), locations = mubins, bin = mubin, min = 0, max = 1)
  plot, mubins + mubin / 2., pdf_mu / max(pdf_mu), xrange = [0, 1], psym = 10, xtitle = 'mu=kzz/kk', yrange = [0, 1]
  oplot, [1., 1.] * mu_HalfW, [0, 1], line = 1
  print, 'mu_HW = ', mu_HalfW, '  Theta_HW = ', acos(mu_HalfW) * 180. / !pi

  ; ONLY obs_i photons
  ; time variations
  DeltaT = tt[obs_i] - t_free[obs_i]
  dtbin = (average(DeltaT) - min(DeltaT)) / 23.
  maxT = average(DeltaT) * 2.3
  Tpdf = rays_histogram_weight(DeltaT, weight = ww, unweighted = Tpdf_noa, locations = tbins, $
    bin = dtbin, min = 0, max = maxT)

  plot, tbins + dtbin / 2., Tpdf / max(Tpdf), xrange = [0, maxT], psym = 10, xtitle = 'Time [s]'
  oplot, tbins + dtbin / 2., Tpdf_noa / float(max(Tpdf_noa)), psym = 10, line = 2, color = 12
  max_Tpdf = max(Tpdf / max(Tpdf), i_maxT)
  t_max = tbins[i_maxT] + dtbin / 2.

  time_1st = average(DeltaT * ww) / average(ww)
  time_2nd = sqrt(average(DeltaT ^ 2 * ww - time_1st ^ 2 * ww) / average(ww))
  ; 1st and 2nd moments for arrival time

  t_HM_before = tbins[min(where(Tpdf / max(Tpdf) ge 0.5))] + dtbin / 2.
  t_HM_after = tbins[max(where(Tpdf / max(Tpdf) ge 0.5))] + dtbin / 2.
  t_1e_after = tbins[max(where(Tpdf / max(Tpdf) ge exp(-1.)))] + dtbin / 2.

  oplot, [1., 1.] * t_max, minmax(Tpdf / max(Tpdf)), line = 2
  oplot, [1., 1.] * t_HM_before, minmax(Tpdf / max(Tpdf)), line = 1
  oplot, [1., 1.] * t_HM_after, minmax(Tpdf / max(Tpdf)), line = 1
  oplot, [1., 1.] * t_1e_after, minmax(Tpdf / max(Tpdf)), line = 1, color = 7

  print, 'Tmax=', t_max, ' Delay =', t_max, '  HM before =', t_HM_before, 'HM_after   =', t_HM_after

  plot, x_im, y_im, psym = 3, xrange = xlimits, yrange = xlimits
  pp = findgen(41) / 40 * !pi * 2
  oplot, cos(pp), sin(pp), thick = 2, color = 5
  oplot, R_init * cos(pp), R_init * sin(pp), line = 2
  R_proj = R_init * sin(rtheta[0])
  oplot, R_proj * cos(pp), R_proj * sin(pp), line = 2, color = 2
  oplot, r_scat * sin(rtheta[0]) * cos(pp), r_scat * sin(rtheta[0]) * sin(pp), line = 2, color = 12

  oplot, [xc, xc], [yc, yc], psym = 7, color = 7
  print, 'With absorption :'
  print, 'X-shift [R_s]=', xc - R_init * sin(rtheta[0]), '+/-', Err_xc, '  Y-shift [R_s]=', yc - 0., '+/-', Err_yc
  print, 'Xsize [R_s] =', sx, '+/-', Err_sx, ' Ysize [R_s]=', sy, '+/-', Err_sy
  print, 'Without absorption :'
  print, 'X-shift [R_s]=', xc_noa - R_init * sin(rtheta[0]), '  Y-shift [R_s]=', yc_noa - 0.
  print, 'Xsize [R_s]=', sx_noa, ' Ysize [R_s]=', sy_noa

  ; with absorption
  xc_time = rays_histogram_weight(DeltaT, weight = x_im * ww, bin = dtbin, min = 0., max = maxT) / $
    (rays_histogram_weight(DeltaT, weight = ww, bin = dtbin, min = 0., max = maxT) + 1e-10)
  yc_time = rays_histogram_weight(DeltaT, weight = y_im * ww, bin = dtbin, min = 0., max = maxT) / $
    (rays_histogram_weight(DeltaT, weight = ww, bin = dtbin, min = 0., max = maxT) + 1e-10)
  xs2_time = rays_histogram_weight(DeltaT, weight = x_im ^ 2 * ww, bin = dtbin, min = 0., max = maxT) / $
    (rays_histogram_weight(DeltaT, weight = ww, bin = dtbin, min = 0., max = maxT) + 1e-10)
  ys2_time = rays_histogram_weight(DeltaT, weight = y_im ^ 2 * ww, bin = dtbin, min = 0., max = maxT) / $
    (rays_histogram_weight(DeltaT, weight = ww, bin = dtbin, min = 0., max = maxT) + 1e-10)
  xs_time = sqrt(xs2_time - xc_time ^ 2) * 2.355
  ys_time = sqrt(ys2_time - yc_time ^ 2) * 2.355
  Err_xc_time = xs_time / sqrt(Tpdf_noa + 1e-8) / 2.355
  Err_yc_time = ys2_time / sqrt(Tpdf_noa + 1e-8) / 2.355

  Err_xs_time = xs_time / sqrt(2. * Tpdf_noa + 1e-8)
  Err_ys_time = ys_time / sqrt(2. * Tpdf_noa + 1e-8)

  ; xc_time=xc_time/(Tpdf+1e-10)
  ; yc_time=yc_time/(Tpdf+1e-10)

  plot, tbins + dtbin / 2., xc_time, xrange = [0, maxT], psym = 10, xtitle = 'Time [s]', $
    ytitle = 'X/Y centroids [R_s]', yrange = minmax([0, xc, yc]) * 2
  errplot, tbins + dtbin / 2., xc_time - Err_xc_time, xc_time + Err_xc_time

  oplot, tbins + dtbin / 2., yc_time, psym = 10, color = 7
  errplot, tbins + dtbin / 2., yc_time - Err_yc_time, yc_time + Err_yc_time, color = 7

  oplot, minmax(tbins + dtbin / 2.), minmax(xc), line = 1
  oplot, minmax(tbins + dtbin / 2.), minmax(yc), line = 1, color = 7
  oplot, minmax(tbins + dtbin / 2.), minmax(min(rxx0)), line = 2, color = 13
  oplot, minmax(tbins + dtbin / 2.), minmax(min(ryy0)), line = 2, color = 13
  oplot, [1., 1.] * (tbins[i_maxT] + dtbin / 2.), minmax([0, xc, yc]) * 2, line = 1

  plot, tbins + dtbin / 2., xs_time, xrange = [0, maxT], psym = 10, xtitle = 'Time [s]', ytitle = 'FWHM X/Y size [R_s]', yrange = minmax([0, sx, sy]) * 2
  errplot, tbins + dtbin / 2., xs_time - Err_xs_time, xs_time + Err_xs_time

  oplot, tbins + dtbin / 2., ys_time, psym = 10, color = 7
  errplot, tbins + dtbin / 2., ys_time - Err_ys_time, ys_time + Err_ys_time, color = 7

  oplot, minmax(tbins + dtbin / 2.), minmax(sx), line = 1
  oplot, minmax(tbins + dtbin / 2.), minmax(sy), line = 1, color = 7
  oplot, [1., 1.] * (tbins[i_maxT] + dtbin / 2.), minmax([0, sx, sy]) * 2, line = 1

  ; plotting profiles along x and y
  r_max = max([[x_im], [y_im]]) * 1.2
  xbin = r_max / 30
  xpdf = histogram(x_im, locations = xloc, bin = xbin, min = -r_max, max = r_max)
  ypdf = histogram(y_im, locations = yloc, bin = xbin, min = -r_max, max = r_max)
  plot, xloc + xbin / 2., xpdf, xrange = xlimits, psym = 10
  oplot, yloc + xbin / 2., ypdf, psym = 10, color = 7

  !p.charsize = 1
  !p.multi = 0
  print, 'FWHM Source sizes  (x & y) [degrees] =', atan(sx * 0.5 / (au - R_init)) * 360. / !pi, $
    atan(sy * 0.5 / (au - R_init)) * 360. / !pi
  print, 'Source Area (ellipse) [arcmin^2]=', !pi * sx * sy / 4. * 16. ^ 2
  ; print, 'Predicted FWHM source size [degrees] =',atan((rint[mloc1]-R_init)*0.5/(au-r_init))*360./!PI
  ; print, 'X-size [arcmin]=', sqrt(average(x_im^2*ww-xc^2*ww)/average(ww)), $
  ; ' Y-size [arcmin]=', sqrt(average(y_im^2*ww-yc^2*ww)/average(ww))
  print, '1/e decay time =', t_1e_after - t_max, 'T1st=', time_1st, 'T2nd', time_2nd

  out = {title: 'ray tracing results in the solar corona', r_init: R_init, r_stop: r_stop, f_pe0: omega_pe0 / 1e6 / 2 / !pi, $
    omega0: omega_pe0, fpe_stop: omega / 1e6 / 2 / !pi, omega: omega / 1e6 / 2 / !pi, t_free: t_free, $
    kxx: kxx, kyy: kyy, kzz: kzz, kk: kk, rxx: rxx, ryy: ryy, rzz: rzz, rr: rr, tt: tt, tau0: tau0, qeps2_scaling: qeps2_scaling, $
    asym: asym, anis: anis, rtheta: rtheta, rfi: rfi, f_ratio: f_ratio, $
    xcen: xc, ycen: yc, xc_err: Err_xc, yc_err: Err_yc, $
    xsize: sx, ysize: sy, xs_err: Err_sx, ys_err: Err_sy, theta_hw: acos(mu_HalfW) * 180. / !pi, $
    t_peak: t_max, t_hm_before: t_HM_before, t_hm_after: t_HM_after, t_1E_after: t_1e_after, terr: dtbin}

  fname = 'fpe' + string(format = '(I7.7)', omega_pe0 / 1e3 / 2 / !pi) + 'kHz_FE' + string(format = '(I7.7)', omega0 / 1e3 / 2 / !pi) $
    + 'anis' + string(format = '(F4.2)', anis) + 'qeps2_s' + string(format = '(F5.3)', qeps2_scaling) $
    + 'asym' + string(format = '(F4.2)', asym) + 'fr' + string(format = '(F4.2)', f_ratio) + '.sav'

  save, out, filename = fname
  print, 'data saved to file: ', fname
  print, 'All completed .............................................................. OK'

  ; stop
end

function nu_scat, r, omega, qeps2_scaling
  compile_opt idl2
  c = 2.998e10
  ; function to run a specific function for nu_scat 
  ; return,nu_s
  ; units of rad^2/s

  nu_s = nu_scat_fixed(r, omega, qeps2_scaling)
  ; changes here to avoid other changes in the code
  return, nu_s
end

function nu_scat_fixed, r, omega, qesp2_scaling
  compile_opt idl2
  ; scattering qeps2 expression as Equation 14 in ApJ 2003 paper
  ; https://doi.org/10.3847/1538-4357/acf6c1
  ; qesp2_scaling was used 1/4,1/2,1,2,4
  c = 2.998d10
  ; l_i =684.*1e5/sqrt(density_r(r))
  w_pe = omega_pe(r)
  q_eps2 = qesp2_scaling* 2e3 / r ^ 0.7 * ((r - 1.) / r) ^ 2.7 / 6.96d10
  nu_s = !pi / 8 * q_eps2 * w_pe ^ 4 * c / omega / (omega ^ 2 - w_pe ^ 2) ^ 1.5
  return, nu_s
  ; units of rad^2/s
end


function omega_pe, r
  compile_opt idl2
  ; returns plasma frequency
  return, 8.98e3 * sqrt(density_r(r)) * 2 * !pi
end

function domega_dr, r
  compile_opt idl2
  ; returns dOmega_pe/dr
  ; INPUT: r distance in R_sun units
  r_tab = 10. ^ (findgen(300) / 100.)
  d_tab = deriv(r_tab, omega_pe(r_tab))
  return, interpol(d_tab, r_tab, r)
end

function density_r, rr
  compile_opt idl2
  ; returns plasma density in cm^-3
  rr = rr * 1d0
  ; input: rr is the solar radius

  ; density in the low corona
  h0 = 144.d / 6.96e5
  ; 144km/solar radius in km

  h1 = 20.d / 960.
  ; 20 arcsec/solar radius

  ; nc=1.17e17*exp(-(rr-1.)/(h0))+1e11*exp(-(rr-1.)/0.02)
  nc = 3d11 * exp(-(rr - 1.d0) / h1)

  ; Newkirk model below
  ; return, 4.2e4*10.^(4.2/rr)

  ; Leblanc 1998 model below
  ; return, 3.3e5*RR^(-2.)+ 4.1e6*RR^(-4.)+8e7*RR^(-6)
  ;
  ; Model that fits Parker model but simple power-laws

  return, 4.8e9 / rr ^ 14. + 3e8 / rr ^ 6. + 1.39d6 / rr ^ 2.3
  ; +nc
end

function bfield_r, r
  compile_opt idl2
  ; input: R [solar radius]
  ; output B [Gauss]
  ; From 2 to 15 solar radii, Pätzold et al. (1987)
  ; obtained from Faraday rotation the global magnetic field variation of
  ; B  =(6/R^3 + 1.18/R^2)G
  if min(r) lt 1.0 then print, 'WARNING: Minimum r should be 1 solar radius'
  r = r * 1.0
  B_high = 6. / r ^ 3 + 1.18 / r ^ 2
  ; For higher corona we take the function from
  ; https://solarscience.msfc.nasa.gov/papers/garyga/PlasmaBeta.pdf
  ; Equation 2 B(r)=
  r_s = 6.996e10
  ; solar radius
  r1 = 0.5 * 1e8 / r_s
  r2 = 75. * 1e8 / r_s
  r3 = 696. * 1e8 / r_s
  B_low = 2500. / (1. + (r - 1.) / r1) ^ 3 + 50. / (1. + (r - 1.) / r2) ^ 3 + 1. / (1. + (r - 1.) / r3) ^ 3
  return, B_high + B_low
end

;+
; NAME:
;    rays_histogram_weight
;    changed from HISTOGRAM_WEIGHT
;
; PURPOSE:
;    Wrapper to the built-in HISTOGRAM function that calculates a
;    weighted histogram.
;
; CATEGORY:
;    Math
;
; CALLING SEQUENCE:
;    Result = HISTOGRAM_WEIGHT(Data)
;
; INPUTS:
;    Data:   Values whose histogram is to be taken.
;
; KEYWORD PARAMETERS:
;    WEIGHT:       A vector of the same length as Data with the weights
;                  for each data value.
;
;    BIN:          Bin width. Passed through to HISTOGRAM.
;
;    UNWEIGHTED:       Outputs the unweighted histogram.
;
;    REVERSE_INDICES:  Outputs the reverse index array.
;
;    _REF_EXTRA:   All extra keywords are passed through to HISTOGRAM.
;
; OUTPUTS:
;    The function returns a histogram where each Data value has been
;    weighted by its WEIGHT value. The return type is the same type
;    as WEIGHT.
;
; EXAMPLE:
;    IDL> values = 0.1*FINDGEN(40)
;    IDL> PRINT, HISTOGRAM_WEIGHT(values, WEIGHT=VALUES, UNWEIGHTED=plainhist)
;          4.50000      14.5000      24.5000      34.5000
;    IDL> PRINT, plainhist
;              10          10          10          10
;
; MODIFICATION HISTORY:
;    Written by:     Jeremy Bailin
;    12 June 2008    Public release in JBIU
;    11 April 2009   Bug fix
;    8 November 2009 Bug fux for bins with no entries
;-
function rays_histogram_weight, DATA, bin = bin, weight = weight, $
  reverse_indices = ri, unweighted = prehist, _ref_extra = histkeywords
  compile_opt idl2

  prehist = histogram(DATA, bin = bin, _strict_extra = histkeywords, reverse_indices = ri)

  histsize = size(prehist, /dimen)

  outhist = replicate(weight[0], histsize)

  for i = 0l, n_elements(prehist) - 1 do if prehist[i] gt 0 then begin
    q = ri[ri[i] : ri[i + 1] - 1]
    outhist[i] = total(weight[q])
  endif else outhist[i] = 0.

  return, outhist
end

pro parker_Bxyz_Earth, r, theta, fi, bxE, byE, bzE
  compile_opt idl2
  ; input r, theta, fi so that z is north x to the Earth (HEEQ coordinates)
  ; output: BzE, ByE, BxE

  B0 = 5e-9 * 1e4 ; Gauss at 1AU
  r0 = 215.
  vsw = 4.2e7 ; solar wind speed 1AU
  ; vsw=1e9 ; 10000 km/s to check radial
  Omega = 2 * !pi / (27.0 * 24 * 3600) ; solar rotation rate

  ; parker spiral in spherical coordinates
  Br = B0 * (r0 / r) ^ 2
  Btheta = 0.
  Bfi = -B0 * (r0 / r) * (Omega * r0 * 6.95e10 / vsw) * sin(theta)

  ; converting to cartesian coordinates
  bxE = Br * sin(theta) * cos(fi) - Bfi * sin(fi)
  byE = Br * sin(theta) * sin(fi) + Bfi * cos(fi)
  bzE = Br * cos(theta)
end