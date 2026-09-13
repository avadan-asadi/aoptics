"""
opticspy.utils
==============
Visualization helpers and common utility functions.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

PHASE_CMAP = LinearSegmentedColormap.from_list('phase',
    ['#3d0066','#0000ff','#00ffff','#00ff00','#ffff00','#ff0000','#3d0066'])

def show_field(field, title='Field', L=1e-3, figsize=(12,4), filename=None):
    I = np.abs(field)**2; ph = np.angle(field)
    extent = [-L*1e3, L*1e3, -L*1e3, L*1e3]
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    im0 = axes[0].imshow(I, extent=extent, origin='lower', cmap='inferno', interpolation='bilinear')
    axes[0].set_title(f'{title} – Intensity'); plt.colorbar(im0, ax=axes[0])
    im1 = axes[1].imshow(ph, extent=extent, origin='lower', cmap=PHASE_CMAP,
                          vmin=-np.pi, vmax=np.pi, interpolation='bilinear')
    axes[1].set_title(f'{title} – Phase'); cbar=plt.colorbar(im1, ax=axes[1])
    cbar.set_ticks([-np.pi,0,np.pi]); cbar.set_ticklabels(['-π','0','π'])
    plt.tight_layout()
    if filename: fig.savefig(filename, dpi=150, bbox_inches='tight')
    return fig, axes

def show_interference(intensity, phase_map, title='Interference', L=1e-3, filename=None):
    extent = [-L*1e3, L*1e3, -L*1e3, L*1e3]
    fig, axes = plt.subplots(1, 2, figsize=(12,5))
    axes[0].imshow(intensity, extent=extent, origin='lower', cmap='gray', interpolation='bilinear')
    axes[0].set_title(f'{title} – Intensity')
    im1 = axes[1].imshow(phase_map, extent=extent, origin='lower', cmap=PHASE_CMAP,
                          vmin=-np.pi, vmax=np.pi, interpolation='bilinear')
    axes[1].set_title(f'{title} – Phase')
    cbar=plt.colorbar(im1, ax=axes[1]); cbar.set_ticks([-np.pi,0,np.pi]); cbar.set_ticklabels(['-π','0','π'])
    plt.tight_layout()
    if filename: fig.savefig(filename, dpi=150, bbox_inches='tight')
    return fig, axes

def show_talbot_carpet(carpet, z_norm, title='Talbot Carpet', filename=None):
    fig, ax = plt.subplots(figsize=(10,6))
    ax.imshow(carpet, aspect='auto', cmap='inferno', interpolation='bilinear',
              extent=[0, carpet.shape[1], z_norm[-1], z_norm[0]])
    ax.set_ylabel('z / z_Talbot'); ax.set_title(title)
    plt.tight_layout()
    if filename: fig.savefig(filename, dpi=150, bbox_inches='tight')
    return fig, ax

def show_wigner(W, xvec, pvec, title='Wigner Function', filename=None):
    vmax = np.max(np.abs(W))
    fig, ax = plt.subplots(figsize=(6,5))
    im = ax.pcolormesh(xvec, pvec, W, cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='auto')
    ax.set_title(title); ax.set_xlabel('x'); ax.set_ylabel('p'); ax.set_aspect('equal')
    plt.colorbar(im, ax=ax); plt.tight_layout()
    if filename: fig.savefig(filename, dpi=150, bbox_inches='tight')
    return fig, ax

def show_photon_statistics(n_vals, P_n, title='Photon Statistics', filename=None):
    fig, ax = plt.subplots(figsize=(8,4))
    ax.bar(n_vals, P_n, color='steelblue', edgecolor='navy', alpha=0.8)
    ax.set_xlabel('Photon number n'); ax.set_ylabel('P(n)'); ax.set_title(title)
    plt.tight_layout()
    if filename: fig.savefig(filename, dpi=150, bbox_inches='tight')
    return fig, ax

def normalize(arr):
    mn, mx = arr.min(), arr.max()
    return (arr - mn) / (mx - mn + 1e-30)

def radial_profile(image, center=None):
    if center is None: center = (image.shape[0]//2, image.shape[1]//2)
    y, x = np.indices(image.shape)
    r = np.sqrt((x-center[1])**2 + (y-center[0])**2).astype(int)
    r_max = r.max()
    profile = np.array([image[r==ri].mean() if (r==ri).any() else 0 for ri in range(r_max+1)])
    return np.arange(r_max+1), profile

def fwhm(y, x=None):
    above = y >= y.max()/2; idx = np.where(above)[0]
    if len(idx)<2: return 0.0
    return float(x[idx[-1]]-x[idx[0]]) if x is not None else float(idx[-1]-idx[0])

def unwrap_phase(phase_map):
    return np.unwrap(np.unwrap(phase_map, axis=1), axis=0)

def power_in_bucket(field, radius_pixels):
    N = field.shape[0]; y,x = np.indices((N,N))
    r = np.sqrt((x-N//2)**2+(y-N//2)**2); I = np.abs(field)**2
    return float(I[r<=radius_pixels].sum()/(I.sum()+1e-30))
