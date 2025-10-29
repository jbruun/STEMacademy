import numpy as np
import math

"""
Dette program gør to ting:
    1. Læser de moderne afstande og heliocentriske hastigheder
       for Hubbles 24 galakser fra filen "hubble_tab1_modern.dat".
       Denne fil skal indeholde to søjler: afstand (Mpc) og v_helio (km/s).
    2. Retter hastighederne for Solens bevægelse i forhold til den kosmiske
       mikrobølgebaggrund (CMB). De korrigerede værdier gemmes i
       "hubble_tab1_modern_solkorr.dat".

Filen "ra_dec.dat" skal indeholde RA og deklination for de 24 galakser
i tekstformat, én pr. linje, samme rækkefølge som datafilen.
"""

# ---------- Hjælpefunktioner ----------

def ra_hms_to_rad(s):
    """Konverterer RA fra 'hhmmss' eller 'hh:mm:ss' format til radianer."""
    s = s.replace('h', ' ').replace('m',' ').replace('s',' ').replace(':',' ')
    parts = s.split()
    h = float(parts[0]); m = float(parts[1]); sec = float(parts[2])
    hours = h + m/60.0 + sec/3600.0
    deg = hours * 15.0
    return math.radians(deg)

def dec_dms_to_rad(s):
    """Konverterer deklination fra '±ddmmss' eller '±dd:mm:ss' til radianer."""
    s = s.replace('d',' ').replace('m',' ').replace('s',' ').replace(':',' ')
    parts = s.split()
    sign = 1
    if parts[0].startswith('-'):
        sign = -1
        parts[0] = parts[0][1:]
    if parts[0].startswith('+'):
        parts[0] = parts[0][1:]
    d = float(parts[0]); m = float(parts[1]); sec = float(parts[2])
    deg = d + m/60.0 + sec/3600.0
    return math.radians(sign*deg)

# ---------- Solens bevægelse i forhold til CMB ----------

# Hastighed og retning (Planck 2018 værdier)
V_sun = 369.82   # km/s
l = math.radians(264.021)   # galaktisk længde
b = math.radians(48.253)    # galaktisk bredde

# Konverter til rektascension/deklination (for brug i samme formel som Hubble)
# Aproksimation: Galaktiske koordinater -> Ækvatorial (J2000) for dipolretning
# Normalt skal der bruges en rotationsmatrix, men vi hardcoder RA,dec for CMB-dipol:
ra_cmb  = math.radians(167.942)  # 11h11m34s
dec_cmb = math.radians(-6.944)   # -6°56' (Planck 2018)

# Enhedsvektor for Solens bevægelse
X = math.cos(ra_cmb)*math.cos(dec_cmb)
Y = math.sin(ra_cmb)*math.cos(dec_cmb)
Z = math.sin(dec_cmb)

# Skaler med hastighed
Vx = V_sun * X
Vy = V_sun * Y
Vz = V_sun * Z

# ---------- Indlæs galakse-data ----------

# Moderne afstande og heliocentriske hastigheder
dist, v_helio = np.loadtxt("hubble_tab1_modern.dat", unpack=True)

# RA, dec fra fil
ra_str, dec_str = np.loadtxt("ra_dec.dat", unpack=True, dtype=str)
ra_rad  = np.array([ra_hms_to_rad(s) for s in ra_str])
dec_rad = np.array([dec_dms_to_rad(s) for s in dec_str])

# ---------- Korrigér hastigheder ----------

v_corr = []
for ra, dec, v in zip(ra_rad, dec_rad, v_helio):
    # Retningsvektor til galaksen
    nx = math.cos(ra)*math.cos(dec)
    ny = math.sin(ra)*math.cos(dec)
    nz = math.sin(dec)
    # Projektion af Solens bevægelse på denne retning
    v_proj = Vx*nx + Vy*ny + Vz*nz
    # Korrigeret hastighed
    v_cmb = v + v_proj
    v_corr.append(v_cmb)

v_corr = np.array(v_corr)

# ---------- Gem resultat ----------
np.savetxt("hubble_tab1_modern_solkorr.dat",
           np.column_stack([dist, v_corr]),
           fmt="%.6f  %.3f",
           header="distance_Mpc  velocity_CMBcorr_km_s")

print("Korrigerede hastigheder gemt i 'hubble_tab1_modern_solkorr.dat'")
