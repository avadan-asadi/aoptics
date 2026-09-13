"""
opticspy.beams - Laser beam generation
"""
import numpy as np
from scipy.special import genlaguerre, hermite, jv, airy as scipy_airy

def make_grid(N=512, L=1e-3):
    """Return x,y,r,phi 2-D grids."""
    x1 = np.linspace(-L, L, N)
    x, y = np.meshgrid(x1, x1)
    r = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return x, y, r, phi

def plane_wave(N=512, L=1e-3, wavelength=633e-9, kx=0.0, ky=0.0):
    x, y, _, _ = make_grid(N, L)
    return np.exp(1j*(kx*x + ky*y))

def gaussian_beam(N=512, L=1e-3, w0=100e-6, wavelength=633e-9, z=0.0):
    x, y, r, _ = make_grid(N, L)
    k = 2*np.pi/wavelength
    zR = np.pi*w0**2/wavelength
    if z == 0:
        return np.exp(-r**2/w0**2).astype(complex)
    wz = w0*np.sqrt(1+(z/zR)**2)
    Rz = z*(1+(zR/z)**2)
    gouy = np.arctan(z/zR)
    return (w0/wz)*np.exp(-r**2/wz**2)*np.exp(1j*(k*z + k*r**2/(2*Rz) - gouy))

def hermite_gaussian(m, n, N=512, L=1e-3, w0=100e-6, wavelength=633e-9, z=0.0):
    x, y, r, _ = make_grid(N, L)
    k = 2*np.pi/wavelength
    zR = np.pi*w0**2/wavelength
    wz = w0*np.sqrt(1+(z/zR)**2) if z!=0 else w0
    Hm = hermite(m)(np.sqrt(2)*x/wz)
    Hn = hermite(n)(np.sqrt(2)*y/wz)
    field = Hm*Hn*np.exp(-(x**2+y**2)/wz**2)
    if z != 0:
        Rz = z*(1+(zR/z)**2)
        gouy = (m+n+1)*np.arctan(z/zR)
        field = field*np.exp(1j*(k*z + k*r**2/(2*Rz) - gouy))
    mx = np.max(np.abs(field))
    return (field/(mx+1e-30)).astype(complex)

def laguerre_gaussian(p, l, N=512, L=1e-3, w0=100e-6, wavelength=633e-9, z=0.0):
    x, y, r, phi = make_grid(N, L)
    k = 2*np.pi/wavelength
    zR = np.pi*w0**2/wavelength
    wz = w0*np.sqrt(1+(z/zR)**2) if z!=0 else w0
    rho = np.sqrt(2)*r/wz
    Lpl = genlaguerre(p, abs(l))(rho**2)
    field = (rho**abs(l))*np.exp(-rho**2/2)*Lpl*np.exp(1j*l*phi)
    if z != 0:
        Rz = z*(1+(zR/z)**2)
        gouy = (2*p+abs(l)+1)*np.arctan(z/zR)
        field = field*np.exp(1j*(k*z + k*r**2/(2*Rz) - gouy))
    mx = np.max(np.abs(field))
    return (field/(mx+1e-30)).astype(complex)

def bessel_beam(m=0, N=512, L=1e-3, wavelength=633e-9, kr=1e4):
    _, _, r, phi = make_grid(N, L)
    return (jv(m, kr*r)*np.exp(1j*m*phi)).astype(complex)

def airy_beam(N=512, L=1e-3, a=0.1, scale=1e-4):
    x, y, _, _ = make_grid(N, L)
    Ai_x, _, _, _ = scipy_airy(x/scale)
    Ai_y, _, _, _ = scipy_airy(y/scale)
    field = Ai_x*Ai_y*np.exp(a*(x+y)/scale)
    mx = np.max(np.abs(field))
    return (field/(mx+1e-30)).astype(complex)

def tophat_beam(N=512, L=1e-3, radius=100e-6):
    _, _, r, _ = make_grid(N, L)
    f = np.zeros((N,N), dtype=complex)
    f[r<=radius] = 1.0
    return f

def intensity(field): return np.abs(field)**2
def phase(field): return np.angle(field)
