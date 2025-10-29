# -*- coding: utf-8 -*-
"""
Scriptet measure_wavelength.py køres bedst fra kommandolinjen. Denne version
er ChatGPT's oversættelse, lavet til at virke i IPython og Jupyter notebooks.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pathlib import Path
from scipy.optimize import curve_fit
import sys

# ----------------------------------------------------------------------------
# Hjælpefunktion: Vent på brugerens input (virker både i IPython og terminal)
# ----------------------------------------------------------------------------
def wait_for_user(fig=None, message="Tryk 'q' i figuren og dernæst Enter i terminalen for at fortsætte..."):
    try:
        from IPython import get_ipython
        in_ipython = get_ipython() is not None
    except Exception:
        in_ipython = False

    if fig is not None:
        done = False
        def on_key(event):
            nonlocal done
            if event.key == 'q':
                done = True
                plt.close(fig)

        fig.canvas.mpl_connect('key_press_event', on_key)

    if in_ipython:
        print(message)
        input()  # wait for Enter in Jupyter cell
        if fig is not None:
            plt.close(fig)
    else:
        print(message)
        plt.show()  # Waits until window closed normally


# ---- Hardcodet input ----
fitsfile    = 'UGC914_skysub.fits'  # .fits-fil med 2D-spektrum
arc_per_pix = 0.163                 # buesekunder pr. pixel
kpc_per_arc = 0.2138                # kpc pr. buesekund
ncenter     = 740                   # Pixelindex for galaksens centrum
binsize     = 15                    # Antal pixels, som integreres over

# ---- Åbn inputfilen ----
hdul = fits.open(fitsfile)
data = hdul[0].data
header = hdul[0].header
hdul.close()

if data.ndim != 2:
    raise ValueError(f"Forventer 2D-billede, men fik {data.ndim}D")

ny, nx = data.shape

# ----------------------------------------------------------------------------
# Figur 1: 2D-spektrum og valg af række
# ----------------------------------------------------------------------------
cdelt1 = header['CDELT1']
crval1 = header['CRVAL1']
wavelength = crval1 + cdelt1 * np.arange(nx, dtype=float)
y_index = np.arange(ny, dtype=float)
distance_kpc = (y_index - ncenter) * arc_per_pix * kpc_per_arc

plt.figure(figsize=(9,6))
im = plt.imshow(data, aspect='auto', origin='lower',
                extent=[wavelength[0], wavelength[-1], distance_kpc[0], distance_kpc[-1]],
                interpolation='nearest', vmin=0, vmax=1000)
cbar = plt.colorbar(im, pad=0.05)
cbar.set_label('Flux / tællinger')
instruktion = "Klik på den ønskede position langs spalten (y-aksen)."
plt.xlabel('Bølgelængde / Å')
plt.ylabel('Afstand langs spalte / kpc')
plt.title(instruktion)
plt.tight_layout()
print('\n' + instruktion + '\n')
plt.show(block=False)

tpoints = plt.ginput(n=1, timeout=0)
nfit = tpoints[0][1]
distance_upper = nfit - binsize/2 * arc_per_pix * kpc_per_arc
distance_lower = nfit + binsize/2 * arc_per_pix * kpc_per_arc
plt.axhline(distance_upper, color='g', lw=0.5)
plt.axhline(distance_lower, color='g', lw=0.5)
plt.draw()

wait_for_user(plt.gcf())
plt.close()

# ----------------------------------------------------------------------------
# Figur 2: 1D-spektret
# ----------------------------------------------------------------------------
nlower = int(distance_upper / (arc_per_pix * kpc_per_arc) + ncenter)
nupper = int(distance_lower / (arc_per_pix * kpc_per_arc) + ncenter)
flux = np.zeros(nx)
for i in range(nlower, nupper):
    flux += data[i, :]

instruktion = "Klik i nærheden af en linje for at vælge zoomområde."
print('\n' + instruktion)
fig, ax = plt.subplots()
ax.plot(wavelength, flux, 'r-', lw=0.25)
ax.set_title(instruktion)
ax.set_ylabel('Signal')
ax.set_xlabel('Bølgelængde / Å')
plt.show(block=False)

tpoints = plt.ginput(n=1, timeout=0)
lambda0, _ = tpoints[0]
ax.axvline(lambda0, color='b', lw=0.8)
plt.draw()
print(f"Du valgte området omkring {lambda0:.2f} Å")
wait_for_user(fig)
plt.close(fig)

# ----------------------------------------------------------------------------
# Figur 3: Zoom ind og klik tre gange
# ----------------------------------------------------------------------------
instruktion = ("Klik tre gange: 1) venstre for linjen, 2) på linjen, "
               "3) højre for linjen.\nTryk derefter 'q'  plottet og Enter i terminalen for at fortsætte.")
print('\n' + instruktion + '\n')
fig, ax = plt.subplots()
ax.plot(wavelength, flux, 'r-', lw=0.25)
ax.set_xlim(lambda0 - 250, lambda0 + 250)
ax.set_ylim(bottom=0)
ax.set_title(instruktion)
plt.show(block=False)

slitregions = []
while len(slitregions) < 3:
    tpoints = plt.ginput(n=1, timeout=0)
    if not tpoints:
        print("Ingen klik registreret – prøv igen.")
        continue
    x, _ = tpoints[0]
    slitregions.append(int(x))
    ax.axvline(x, color='b', lw=0.5)
    fig.canvas.draw()

lambdal, lambdac, lambdar = slitregions[:3]
wait_for_user(fig)
plt.close(fig)

# ----------------------------------------------------------------------------
# Figur 4: Fit Gauss-funktion
# ----------------------------------------------------------------------------
def gauss_function(x, I, mu, sigma, background):
    return I/(np.sqrt(2*np.pi)*sigma)*np.exp(-0.5*((x-mu)/sigma)**2)+background

fitrange = np.where((wavelength >= lambdal) & (wavelength <= lambdar))[0]
wlguess = lambdac
parameter_guesses = [np.max(flux[fitrange]), wlguess, 4, 2000]
val, cov = curve_fit(gauss_function, wavelength[fitrange], flux[fitrange],
                     p0=parameter_guesses,
                     bounds=([0, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]))

print('------------------------------------------------')
print(f'Afstand fra galaksens centrum: {nfit:.3f} kpc')
print(f'Linjecentrums bølgelængde:     {val[1]:.6f} ± {np.sqrt(cov[1][1]):.2f} Å')
print(f'Linjebredde:                   {val[2]:.2f} Å')
print('------------------------------------------------')

# Plot resultat
wavelength_hires = np.linspace(wavelength[0], wavelength[-1], 10*len(wavelength))
plt.figure()
plt.plot(wavelength, flux, 'r-', label='Observeret spektrum')
plt.plot(wavelength_hires, gauss_function(wavelength_hires, *val), 'b--', label='Gauss-fit')
plt.ylabel('Signal')
plt.xlabel('Bølgelængde / Å')
plt.legend()
plt.xlim(lambda0 - 250, lambda0 + 250)
plt.ylim(0, val[3] + val[0]/(np.sqrt(2*np.pi)*val[2])*1.2)
plt.show()
