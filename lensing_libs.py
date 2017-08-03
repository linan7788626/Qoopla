import numpy as np

#--------------------------------------------------------------------
def gauss_2d(x, y, par):

    xo = par[0]
    yo = par[1]
    amplitude = par[3]
    sigma_x = par[2]/par[4]#*0.693
    sigma_y = par[2]*par[4]#*0.693
    theta = np.deg2rad(par[5])
    offset = 0.0

    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    res = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                            + c*((y-yo)**2)))
    return res

def lensing_signals_sie(x1, x2, lpar):
    # lpar = np.asarray([l_xcen,l_ycen,l_re,l_rc,l_axrat,l_pa])
    # x coordinate of the center of lens (in units of Einstein radius).
    xc1 = lpar[0]
    # y coordinate of the center of lens (in units of Einstein radius).
    xc2 = lpar[1]
    re = lpar[2]   # Einstein radius of lens.
    rc = lpar[3]   # Core size of lens (in units of Einstein radius).
    q = lpar[4]   # Axis ratio of lens.
    phl = lpar[5]   # Orintation of lens.

    phirad = np.deg2rad(phl)
    cosa = np.cos(phirad)
    sina = np.sin(phirad)

    xt1 = (x1 - xc1) * cosa + (x2 - xc2) * sina
    xt2 = (x2 - xc2) * cosa - (x1 - xc1) * sina

    phi = np.sqrt(xt2 * xt2 + xt1 * q * xt1 * q + rc * rc)
    sq = np.sqrt(1.0 - q * q)
    pd1 = phi + rc / q
    pd2 = phi + rc * q
    fx1 = sq * xt1 / pd1
    fx2 = sq * xt2 / pd2
    qs = np.sqrt(q)

    a1 = qs / sq * np.arctan(fx1)
    a2 = qs / sq * np.arctanh(fx2)

    xt11 = cosa
    xt22 = cosa
    xt12 = sina
    xt21 = -sina

    fx11 = xt11 / pd1 - xt1 * \
        (xt1 * q * q * xt11 + xt2 * xt21) / (phi * pd1 * pd1)
    fx22 = xt22 / pd2 - xt2 * \
        (xt1 * q * q * xt12 + xt2 * xt22) / (phi * pd2 * pd2)
    fx12 = xt12 / pd1 - xt1 * \
        (xt1 * q * q * xt12 + xt2 * xt22) / (phi * pd1 * pd1)
    fx21 = xt21 / pd2 - xt2 * \
        (xt1 * q * q * xt11 + xt2 * xt21) / (phi * pd2 * pd2)

    a11 = qs / (1.0 + fx1 * fx1) * fx11
    a22 = qs / (1.0 - fx2 * fx2) * fx22
    a12 = qs / (1.0 + fx1 * fx1) * fx12
    a21 = qs / (1.0 - fx2 * fx2) * fx21

    rea11 = (a11 * cosa - a21 * sina) * re
    rea22 = (a22 * cosa + a12 * sina) * re
    rea12 = (a12 * cosa - a22 * sina) * re
    rea21 = (a21 * cosa + a11 * sina) * re

    #kappa = 0.5 * (rea11 + rea22)
    #shear1 = 0.5 * (rea12 + rea21)
    #shear2 = 0.5 * (rea11 - rea22)

    y11 = 1.0 - rea11
    y22 = 1.0 - rea22
    y12 = 0.0 - rea12
    y21 = 0.0 - rea21

    jacobian = y11 * y22 - y12 * y21
    mu = 1.0 / jacobian

    alpha1 = (a1 * cosa - a2 * sina) * re
    alpha2 = (a2 * cosa + a1 * sina) * re

    return alpha1, alpha2, mu#, kappa, shear1, shear2, mu

def lensed_images(x1, x2, lpar, gpar):
    al1, al2, mu = lensing_signals_sie(x1,x2,lpar)
    yi1 = x1-al1
    yi2 = x2-al2

    g_srcsimage = 0.0*yi1
    g_srcsimage = gauss_2d(x1,x2,gpar)

    g_lensimage = 0.0*yi1
    g_lensimage = gauss_2d(yi1,yi2,gpar)

    return g_lensimage, g_srcsimage, mu, yi1, yi2
##--------------------------------------------------------------------
#if __name__ == '__main__':
    #boxsize = 8.0 # (arcsec)
    #nnn = 128
    #dsx = boxsize/nnn

    #xi1 = np.linspace(-boxsize/2.0,boxsize/2.0,nnn)+0.5*dsx
    #xi2 = np.linspace(-boxsize/2.0,boxsize/2.0,nnn)+0.5*dsx
    #xi1,xi2 = np.meshgrid(xi1,xi2)
    ##----------------------------------------------------------------------
    #l_xcen = 0.0    # x position of center (also try (0.0,0.14)
    #l_ycen = 0.0    # y position of center
    #l_re = 2.0   # Einstein radius of lens.
    #l_rc = 0.0   # Core size of lens (in units of Einstein radius).
    #l_axrat = 0.7   # Axis ratio of lens.
    #l_pa = 49.0   # Orintation of lens.

    #lpar = np.asarray([l_xcen,l_ycen,l_re,l_rc,l_axrat,l_pa])
    ##----------------------------------------------------------------------
    #g_amp = 1.0     # peak brightness value
    #g_sig = 0.1     # Gaussian "sigma" (i.e., size)
    #g_xcen = 0.0    # x position of center (also try (0.0,0.14)
    #g_ycen = 0.0    # y position of center
    #g_axrat = 0.7   # minor-to-major axis ratio
    #g_pa = 45.0     # major-axis position angle (degrees) c.c.w. from x axis

    #gpar = np.asarray([g_xcen,g_ycen,g_sig,g_amp,g_axrat,g_pa])
    ##----------------------------------------------------------------------
    #g_source = 0.0*xi1
    #g_limage = 0.0*xi1
    #mua = 0.0*xi1
    #g_limage,g_source,mua,yi1,yi2 = lensed_images(xi1,xi2,lpar,gpar)

    ##--------------------------lens images contour------------------------
    #levels = [0.5,]
    #pl.figure(num=None,figsize=(10,5),dpi=80, facecolor='w', edgecolor='k')

    #a = pl.axes([0.05,0.1,0.4,0.8])
    #a.set_xlim(-boxsize/2.0,boxsize/2.0)
    #a.set_ylim(-boxsize/2.0,boxsize/2.0)
    #a.contour(xi1,xi2,g_source,levels,colors=('b'))
    #a.contour(yi1,yi2,mua,0,colors=('g'),linewidths = 2.0)

    #b = pl.axes([0.55,0.1,0.4,0.8])
    #b.set_xlim(-boxsize/2.0,boxsize/2.0)
    #b.set_ylim(-boxsize/2.0,boxsize/2.0)
    #b.contour(xi1,xi2,g_limage,levels,colors=('b'))
    #b.contour(xi1,xi2,mua,colors=('r'),linewidths = 2.0)
    #pl.show()
