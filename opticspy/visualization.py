"""
visualization.py - ابزارهای نمایش تصویر
Plotting utilities for optics simulations.
"""
import numpy as np
import matplotlib.pyplot as plt


def show_field(field, title='Field', figsize=(12, 4), cmap_int='hot',
               cmap_phase='hsv', show_phase=True, extent=None, save=None):
    I = np.abs(field)**2
    phi = np.angle(field)
    n_plots = 2 if show_phase else 1
    fig, axes = plt.subplots(1, n_plots, figsize=figsize)
    if n_plots == 1:
        axes = [axes]
    im0 = axes[0].imshow(I, cmap=cmap_int, origin='lower', extent=extent)
    axes[0].set_title(f'{title} - Intensity')
    plt.colorbar(im0, ax=axes[0], label='I [a.u.]')
    if show_phase:
        im1 = axes[1].imshow(phi, cmap=cmap_phase, origin='lower',
                             vmin=-np.pi, vmax=np.pi, extent=extent)
        axes[1].set_title(f'{title} - Phase')
        plt.colorbar(im1, ax=axes[1], label='Phase [rad]')
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches='tight')
    return fig, axes


def show_interference(I, phase, title='Interference', figsize=(12,4), extent=None, save=None):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    im1 = ax1.imshow(I, cmap='gray', origin='lower', extent=extent)
    ax1.set_title(f'{title} - Intensity')
    plt.colorbar(im1, ax=ax1)
    im2 = ax2.imshow(phase, cmap='RdBu', origin='lower', vmin=-np.pi, vmax=np.pi, extent=extent)
    ax2.set_title(f'{title} - Phase Map')
    plt.colorbar(im2, ax=ax2, label='rad')
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches='tight')
    return fig


def show_talbot_carpet(carpet, x, z, z_T, title='Talbot Carpet', save=None):
    fig, ax = plt.subplots(figsize=(10, 6))
    x_mm = x * 1e3; z_mm = z * 1e3
    extent = [x_mm.min(), x_mm.max(), z_mm.min(), z_mm.max()]
    ax.imshow(carpet, aspect='auto', origin='lower', cmap='inferno', extent=extent)
    for n in range(1, int(z.max()/z_T) + 1):
        ax.axhline(n * z_T * 1e3, color='cyan', linewidth=0.8, linestyle='--')
        ax.axhline((n-0.5)*z_T*1e3, color='lime', linewidth=0.5, linestyle=':', alpha=0.7)
    ax.set_xlabel('x [mm]'); ax.set_ylabel('z [mm]')
    ax.set_title(title + f'  (z_T = {z_T*1e3:.3f} mm)')
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches='tight')
    return fig


def show_moire(I_moire, T1, T2, title='Moiré Pattern', figsize=(15,4), save=None):
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    for ax, img, t in zip(axes, [T1, T2, I_moire], ['Grating 1', 'Grating 2', 'Moiré']):
        ax.imshow(img, cmap='gray', origin='lower'); ax.set_title(t); ax.axis('off')
    plt.suptitle(title); plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches='tight')
    return fig


def show_wigner(W, alpha_re, alpha_im, title='Wigner Function', save=None):
    fig, ax = plt.subplots(figsize=(6, 5))
    lim = np.max(np.abs(alpha_re))
    extent = [-lim, lim, -lim, lim]
    vmax = np.max(np.abs(W))
    im = ax.imshow(W, cmap='RdBu', origin='lower', extent=extent, vmin=-vmax, vmax=vmax)
    ax.set_xlabel('Re(α)'); ax.set_ylabel('Im(α)'); ax.set_title(title)
    plt.colorbar(im, ax=ax, label='W(α)')
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches='tight')
    return fig


def show_photon_statistics(n, P, Q=None, g2=None, state_label='', save=None):
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(n, P, color='steelblue', alpha=0.8, edgecolor='k', linewidth=0.5)
    ax.set_xlabel('Photon number n'); ax.set_ylabel('P(n)')
    info = state_label
    if Q is not None: info += f'  Q = {Q:.3f}'
    if g2 is not None: info += f'  g²(0) = {g2:.3f}'
    ax.set_title(f'Photon Statistics  {info}')
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches='tight')
    return fig


def show_beam_profiles_comparison(fields, labels, x=None, figsize=(15,4), save=None):
    N = len(fields)
    fig, axes = plt.subplots(2, N, figsize=figsize)
    for i, (E, label) in enumerate(zip(fields, labels)):
        I = np.abs(E)**2; phi = np.angle(E)
        extent = None
        if x is not None:
            xmm = x*1e3; extent = [xmm.min(), xmm.max(), xmm.min(), xmm.max()]
        axes[0,i].imshow(I, cmap='hot', origin='lower', extent=extent)
        axes[0,i].set_title(label); axes[0,i].axis('off')
        axes[1,i].imshow(phi, cmap='hsv', origin='lower', vmin=-np.pi, vmax=np.pi, extent=extent)
        axes[1,i].axis('off')
    plt.suptitle('Beam Profile Comparison', fontsize=14)
    plt.tight_layout()
    if save:
        plt.savefig(save, dpi=150, bbox_inches='tight')
    return fig
