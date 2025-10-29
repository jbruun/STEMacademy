import numpy as np
import math

"""
Denne funktion gør to ting:
    1. Først beregner den Solens hastighed gennem rummet ud fra Hubbles data.
       Disse data skal ligge i en fil kaldet "hubble_tab1.dat" med to søjler,
       afstand i Mpc og hastighed i km/s, og én række for hver galakse.
       Til beregningen skal også bruges galaksernes koordinater. Disse
       koordinater skal gives i en fil kaldet "ra_dec.dat" ligeledess med 24
       rækker og to søjler med galaksernes rektascension of deklination
       (se forklaring nedenfor).
    2. Dernæst trækker den denne hastighed fra Hubbles observerede hastigheder.
       Resultatet gemmes i en fil kaldet "hubble_tab1_solkorr.dat".

Den kosmiske referenceramme
---------------------------
Når vi måler himmellegemernes hastighed, afhænger resultatet af vores egen
hastighed. Denne hastighed må fratrækkes de målte hastigheder for at de kan
give mening i en kosmologisk sammenhæng, som ikke afhænger af vores synspunkt.

Nu om dage kan vi måle vores hastighed ift. den kosmiske
mikrobølgebaggrundsstråling (CMB), men den kendte Hubble intet til.

En anden måde at måle vores hastighed på er at tage gennemsnittet af alle
galakser i Universet. Hubble havde desværre kun 24, så han måtte foretage en
mindste-kvadrat-regression til disse data.

Ækvatorialkoordinater (RA og dek.)
----------------------------------
Solens hastighed gennem rummet beregnes i det koordinatsystem,
som
    1. man normalt udtrykker himmellegemers koordinater i, og som
    2. ligger stille ift. de 24 galakser, som Hubble målte afstande til.

Himmellegemers position på himlen gives typisk i "ækvatorialkoordinater", dvs.
en vinkel langs ækvator, som måles i timer, minutter og sekunder mellem 0 og
24h, samt en vinkel væk fra ækvator, som måles i grader, bueminutter og
buesekunder mellem -90 og +90 grader.

Den første vinkel kaldes rektascension, og den anden deklination.
    * Rektascension skrives også "RA" eller "α" (alpha).
    * Deklination skrives også "dec"/"dek" eller "δ" (delta).

Ækvatorialkoordinater (x, y og z)
---------------------------------
Vi kan i stedet udtrykke koordinaterne som x-, y- og z-koordinater, hvor
   * x-aksen ligger i ækvatorialplanen (δ = 0°) og peger mod α = 0h, δ = 0°.
   * y-aksen ligger i samme plan, men er vinkelret på x og peger mod
     α = 6h, δ = 0° (dvs. 90° øst for Forårspunktet), og
   * z-aksen peger langs nordpolsaksen, dvs. δ = +90°.

Resultat
----------------
Resultatet for Solens hastighed giver Hubble på side 3 i sin 1929-artikel:
    * Hastighed i x-retningen er X =  -65 ± 50 km/s
    * Hastighed i y-retningen er Y = +226 ± 95 km/s
    * Hastighed i z-retningen er Z = -195 ± 40 km/s
    * Den totale hastighed er V0   = √(65² + 226² + 195²) = 306 km/s.

Med dette program fås resultater som er lidt anderledes, men konsistente.
Forskellene skyldes nok afrundinger og evt. lidt andre værdier for galaksernes
koordinater (som Hubble ikke angiver i artiklen).

Resultatet for galaksernes Sol-hastigheds-korrigerede hastigheder angiver
Hubble ikke i artiklen, men viser dem i sin figur. Du kan vise dem med

> python fit_hubble.py hubble_tab1_solkorr.dat

og eventuelt sammenligne med de ikke-korrigerede hastigheder med

> python fit_hubble.py hubble_tab1.dat
"""

# Funktion som omregner fra rektascension til radianer
def ra_hms_to_rad(s):
    """
    Denne funktion bruges til at omregne en rektascension til radianer.

    Funktionen tager som input en tekststreng, f.eks. '13h 37m 00.9s'
    (mellemrum er ikke nødvendigt, og dens output er så et tal mellem 0 og 2π.

    Eksempel:
    >>> ra_hms_to_rad('13h37m00.9s')
    3.564900447045368
    """
    s = s.replace('h', ' ').replace('m',' ').replace('s',' ').replace(':',' ')
    parts = s.split()
    h = float(parts[0]); m = float(parts[1]); sec = float(parts[2])
    hours = h + m/60.0 + sec/3600.0
    deg = hours * 15.0
    return math.radians(deg)

# Funktion som omregner fra deklination til radianer
def dec_dms_to_rad(s):
    """
    Denne funktion bruges til at omregne en deklination til radianer.

    Funktionen tager som input en tekststreng, f.eks. '-29d 51m 56.7s'
    (mellemrum er ikke nødvendigt, og dens output er så et tal mellem -π/2 og +π/2.

    Eksempel:
    >>> dec_dms_to_rad('-29d51m56.7s')
    -0.5212556710774965
    """
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

# Load data for de 24 galakser i Hubbles Tab. 1. Filen "hubble_tab1.dat" skal
# indeholde søjler, én med afstande i Mpc og én med hastigheder i km/s.
dist,vel = np.loadtxt('hubble_tab1.dat',unpack=True)

# Load de 24 galaksers koordinater. Filen skal indeholde rektascension og
# deklination i to søjler, ordnet på samme måde som hubble_tab1.dat.
ra,dec = np.loadtxt('ra_dec.dat',unpack=True,dtype=str)

ra_rad  = np.array([ra_hms_to_rad(s) for s in ra])
dec_rad = np.array([dec_dms_to_rad(s) for s in dec])

# Opstil en matrix for de ubekendte [K, X, Y, Z]:
A = np.column_stack([dist, np.cos(ra_rad)*np.cos(dec_rad),
                           np.sin(ra_rad)*np.cos(dec_rad),
                           np.sin(dec_rad)])

# Løs med mindste-kvadrat regression
params, residuals, rank, s = np.linalg.lstsq(A, vel, rcond=None)
K,X,Y,Z = params # K er dét vi i dag kalder Hubble-konstanten, H0.

# Estimér 1 sigma-usikkerheder fra residuals og kovarians
N = len(vel); p = 4
if residuals.size > 0:
    sigma2 = residuals[0]/(N-p)
    cov = np.linalg.inv(A.T @ A) * sigma2
    err = np.sqrt(np.diag(cov))
else:
    err = np.full(4, np.nan)
sigmaK,sigmaX,sigmaY,sigmaZ = err

# Total fart og usikkerhed
V0      = np.sqrt(X**2 + Y**2 + Z**2)      # Den totale fart er summen af de tre komponenter i kvadratur
sigmaV0 = 1/V0 * np.sqrt(X**2 * sigmaX**2
                       + Y**2 * sigmaY**2
                       + Z**2 * sigmaZ**2) # Denne approksimation antager, at der ikke er nogen kovarians mellem usikkerhederne på X, Y og Z

print('Bedste fit for Solens hastighed er')
print('  X = {:.1f} ± {:.1f} km/s'.format(X,sigmaX))
print('  Y = {:.1f} ± {:.1f} km/s'.format(Y,sigmaY))
print('  Z = {:.1f} ± {:.1f} km/s'.format(Z,sigmaZ))
print('Total hastighed (V0) = {:.1f} ± {:.1f} km/s'.format(V0,sigmaV0))

# Beregn Solens projektion på hver galakse
v_proj = X * np.cos(ra_rad) * np.cos(dec_rad) \
       + Y * np.sin(ra_rad) * np.cos(dec_rad) \
       + Z * np.sin(dec_rad)

# Korrigerede hastigheder
vel_corr = vel - v_proj

# Gem resultater i en ny fil
outdata = np.column_stack([dist, vel_corr])
np.savetxt("hubble_tab1_solkorr.dat", outdata,
           header="Afstand[Mpc]  V_solkorr[km/s]",
           fmt="%12.6f %12.6f")

print("Korrigerede hastigheder gemt i 'hubble_tab1_solkorr.dat'")
