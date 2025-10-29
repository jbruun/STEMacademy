# -*- coding: utf-8 -*-
"""
Denne funktion finder den centrale bølgelængde af en emissionslinje i et 2D-spektrum,
hvis data ligger i en .fits-fil med navnet givet ved parameteren `fitsfile`.
Udover denne parameter er et par stykker mere hardcodet (se i scriptet under `fitsfile`).

Først vises 2D-spektret med:
 - y-akse: afstand i kpc (distance = (y - ncenter) * arc_per_pix * kpc_per_arc)
 - x-akse: bølgelængde i Ångström (np.arange(crval1, x*cdelt1 + crval1, cdelt1))
I dette plot vælges en position langs galaksens storakse (y-aksen).

Dernæst vises hele spektret i dette område. Klik på den linje, du gerne vil undersøge, f.eks. Hα,
hvis hvilebølgelængde ligger ved 6565 Ångström, men som i spektret er rødforskudt med en faktir (1+z).

Der zoomes nu ind på den valgt linje, og linjen markeres mere præcis på følgende måde:
Først klikkes et godt stykke til venstre for linjen, dernæst på selve linjen, og til sidst et stykke til højre herfor.

Herefter fittes en Gauss-kurve til linje, og resultatet plottes og skrives ud i terminalen.
"""

# Load diverse Python-pakker:
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pathlib import Path
from scipy.optimize import curve_fit

# ---- Hardcodet input ----
fitsfile    = 'UGC914_skysub.fits'  # .fits-fil med 2D-spektrum
arc_per_pix = 0.163                 # Antal buesekunder pr. pixel langs positionsaksen
kpc_per_arc = 0.2138                # Antal kiloparsec i den givne afstand, som ét buesekund dækker
ncenter     = 740                   # Pixelindex for galaksens centrum
binsize     = 15                    # Antal pixels i afstands-retningen, som linjen integreres over

# ---- Åbn inputfilen ----
hdul = fits.open(fitsfile)  # Åbn .fits-filen
data = hdul[0].data         # Gem selve spektret i variablen `data`
header = hdul[0].header     # Variablen `header` indeholder diverse meta-data
hdul.close()                # Luk filen igen, ellers bliver der ballade
# Check at spektret er 2D: (ny, nx)
if data.ndim != 2:
    raise ValueError(f"Forventer 2D-billede, men fik {data.ndim}D")
ny, nx = data.shape         # Antal pixels i hhv. y- og x-retningen (hhv. afstand og bølgelængde)


# ----------------------------------------------------------------------------
# Figur 1: Plot 2D-spektrum og vælg række
# ----------------------------------------------------------------------------

# ---- Byg akser til første plot ----
# np.arange inkluderer start og ekskluderer slut hvis det ikke lander præcist;
# for at sikre længden er nx, laver vi eksplicit array:
cdelt1 = header['CDELT1']                                       # Bølgelængdeskala i Ångström pr. pixel
crval1 = header['CRVAL1']                                       # Start-bølgelængde i Ångström
wavelength = crval1 + cdelt1 * np.arange(nx, dtype=float)       # Bølgelængde-akse
y_index = np.arange(ny, dtype=float)                            # 0...ny-1
distance_kpc = (y_index - ncenter) * arc_per_pix * kpc_per_arc  # Afstand langs spalten (kpc)

# ---- Lav plottevinduet ----
plt.figure(figsize=(9,6))
im = plt.imshow(data, aspect='auto',
                origin='lower', # origin='lower' så mindre y-index (typisk top) placeres nederst i kpc-skalaen
                extent=[wavelength[0], wavelength[-1], distance_kpc[0], distance_kpc[-1]],  # Brug extent så akserne bliver i fysiske enheder
                interpolation='nearest', vmin=0, vmax=1000)
cbar = plt.colorbar(im, pad=0.05)
cbar.set_label('Flux / tællinger')
instruktion = "Klik på den ønskede position langs spalten/galaksen (y-aksen).\nTryk 'q' for at fortsætte."
plt.xlabel('Bølgelængde / Å')
plt.ylabel('Afstand langs spalte / kpc')
plt.title(instruktion)
plt.tight_layout()
print('\n'+instruktion+'\n')

# ---- Sæt op til at modtage input fra musen ----
slitregions = list(range(1))
nregions = 0
while nregions <= 0:
    tpoints = plt.ginput(n=1, timeout=30, show_clicks=True, mouse_add=1, mouse_stop=2)
    pix_ref = tpoints[0]
    slitregions[nregions] = pix_ref[1]
    nregions = nregions+1
    distance_upper = slitregions[0]-binsize/2 * arc_per_pix * kpc_per_arc
    distance_lower = slitregions[0]+binsize/2 * arc_per_pix * kpc_per_arc
    if (nregions == 1):
           plt.axhline(distance_upper, xmin=-1., xmax=2000., color='g', lw=0.5)
           plt.axhline(distance_lower, xmin=-1., xmax=2000., color='g', lw=0.5)
    plt.draw()
plt.show()
nfit = slitregions[0]


# ----------------------------------------------------------------------------
# Figur 2: Lav 1D-spektret for den valgte bin
# ----------------------------------------------------------------------------
flux = np.array(nx) 
nlower = int(distance_upper / (arc_per_pix * kpc_per_arc) + ncenter)
nupper = int(distance_lower / (arc_per_pix * kpc_per_arc) + ncenter)
for i in range(nlower, nupper): 
    flux = flux+data[i,:]

# ---- Lav plottevinduet ----
instruktion = "Klik i nærheden af en linje for at vælge et område at zoome ind på."
print('\n'+instruktion)
fig, ax = plt.subplots()
ax.plot(wavelength, flux, 'r-',lw=.25)
ax.set_title(instruktion)
ax.set_ylabel('Signal')
ax.set_xlabel('Bølgelængde / Å')

# ---- Sæt op til at modtage input fra musen ----
tpoints = plt.ginput(n=1, timeout=0)
lambda0, _ = tpoints[0]
ax.axvline(lambda0, color='b', lw=0.8)
plt.draw()
plt.title(instruktion)
instruktion = "Du valgte området omkring {:.2f} Å\nTryk 'q' for at fortsætte.".format(lambda0)
plt.text(0.1, 0.9,instruktion, horizontalalignment='left',
    verticalalignment='center', transform=ax.transAxes)
print(instruktion.format(lambda0))
done = False
def on_key(event):
    global done
    if event.key == 'q':
        done = True
        plt.close()
fig.canvas.mpl_connect('key_press_event', on_key)
plt.show()  # blocks until window closed


# ----------------------------------------------------------------------------
# Figur 3: Zoom ind på linjen og markér hhv. kontinuet til venstre, selve linjen, og kontinuet til højre
# ----------------------------------------------------------------------------

# ---- Lav plottevinduet ----
instruktion = ("Klik tre gange med musen: 1) venstre for linjen (baggrund),\n"
               "2) på linjen, og 3) højre for linjen (baggrund).\n"
               "Tryk derefter 'q' for at fortsætte.")
print('\n'+instruktion + '\n')
fig, ax = plt.subplots()
ax.plot(wavelength, flux, 'r-', mec='k')
ax.set_ylabel('Signal')
ax.set_xlabel('Bølgelængde [Å]')
ax.set_xlim(lambda0 - 250, lambda0 + 250)
ax.set_ylim(bottom=0)
ax.set_title(instruktion)
plt.show(block=False)   # display the window, but don't block Python

# ---- Sæt op til at modtage input fra musen ----
slitregions = list(range(6))
nregions = 0
while nregions <= 2:
    tpoints = plt.ginput(n=1, timeout=0, show_clicks=True)
    if not tpoints:
        print("Ingen klik registreret – prøv igen.")
        continue
    pix_ref, _ = tpoints[0]
    slitregions[nregions] = int(pix_ref)
    ax.axvline(slitregions[nregions], ymin=-1., ymax=2., color='b', lw=0.5)
    fig.canvas.draw()
    nregions += 1
lambdal, lambdac, lambdar = slitregions[:3]
plt.show(block=True) # After interaction, optionally pause to show result or close window:


# ----------------------------------------------------------------------------
# Figur 4: Fit en Gauss-funktion til linjen, plot det, og udskriv resultatet
# ----------------------------------------------------------------------------

# ---- Definér en Gauss-funktion som ikke (nødvendigvis) gåd mod nul ----
def gauss_function(x,I,mu,sigma,background):
    return I/(np.sqrt(2.0*np.pi)*sigma)*np.exp(-0.5*((x-mu)/sigma)**2)+background

# ---- Udfør fittet ----
fitrange = np.where((wavelength >= lambdal) & (wavelength <= lambdar))[0]
wlguess = lambdac
parameter_guesses = [np.max(flux[fitrange]),wlguess,4,2000] #4 og 2000 er hhv. bredden af linjen og baggrundsniveauet
val, cov = curve_fit(gauss_function, wavelength[fitrange], flux[fitrange],
                     p0=parameter_guesses,
                     bounds=([0, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf])) # Denne linje er for at sikre positiv linjebredde

#Udskriv resultatet
print('------------------------------------------------')
print('Afstand fra galaksens centrum: {:.3}'.format(nfit)+' kpc')
print('Linjecentrums bølgelængde:     {:.6} ± {:.2}'.format(val[1], np.sqrt(cov[1][1])) +' Å')
print('Linjebredde:                   {:.2}'.format(val[2])+' Å')
print('------------------------------------------------')

# Plot fittet oven på 1d spektret
wavelength_hires = np.linspace(wavelength[0], wavelength[-1], 10*len(wavelength)) # Spektret er lav opløsning, så vi laver en higher-res version til at plotte fittet på. Det er pænere.
plt.figure()
plt.plot(wavelength, flux,'r-',mec = 'k',label='Observeret spektrum')
plt.plot(wavelength_hires, gauss_function(wavelength_hires, *val),'b--', label='Gauss-fit')
plt.ylabel('Signal')
plt.xlabel('Bølgelængde / Å')
plt.legend(loc='upper left', fancybox=True, shadow=True, framealpha= 1,\
    facecolor='#d8dcd6',edgecolor='black', prop={'size': 8}, markerscale=2)
plt.xlim(6500,6800)
plt.ylim(0,val[3]+val[0]/(np.sqrt(2.0*np.pi)*val[2])*1.2)
plt.show()
