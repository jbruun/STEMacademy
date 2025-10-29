# Load først en række Python-pakker med forskellige redskaber:
import numpy as np                      # Matematik- og array-redskaber
import matplotlib.pyplot as plt         # Plotteredskaber
from scipy.optimize import curve_fit    # Statistiske redskaber, inkl. fitte-funktionen "curve_fit"

# Giv navnet på fil med to søjler: hastigheder (i km/s) og afstande (i Mpc):
inputfil = "hubble_tab1.dat"                  # <-- dette filnavn kan ændres til en anden fil
# inputfil = "hubble_tab1_solkorr.dat"        # <-- dette filnavn kan ændres til en anden fil
# inputfil = "hubble_tab1_modern.dat"         # <-- dette filnavn kan ændres til en anden fil
# inputfil = "hubble_tab1_modern_solkorr.dat" # <-- dette filnavn kan ændres til en anden fil
# inputfil = "cosmicflows.dat"                # <-- dette filnavn kan ændres til en anden fil
# inputfil = "cosmicflows_adj.dat"            # <-- dette filnavn kan ændres til en anden fil
# inputfil = "pantheon_sh0es.dat"             # <-- dette filnavn kan ændres til en anden fil

# Load dataene fra inputfilen:
dist,vel = np.loadtxt(inputfil,unpack=True)


# Definér en lineær funktion til regression:
def linear(x,a): # Ingen `b`, da vi vil tvinge fittet til at skære (0,0)
    return(a*x)


# Fit en lineær funktion til dataene (altså "lineær regression"):
popt,pcov = curve_fit(linear,dist,vel)   # Input er (funktion, x-værdier, y-værdier); output er (hældning, varians på hældning)
H0        = popt[0]                      # curve_fit returnerer arrays af parametre, men i dette tilfælde er der kun ét tal pr. array,
varH0     = pcov[0][0]                   #   fordi den lineære funktion kun har én uafhængig variabel
sigH0     = np.sqrt(varH0)               # Usikkerheden (1 sigma) er kvadratroden af variansen
d_akse    = np.linspace(0,1.1*max(dist)) # Definér en x-akse for afstande til at plotte fittet på
bestfit   = linear(d_akse,H0)            # Linje for bedste fit, baseret på koefficienterne i `coeff`


# Plot data og bedste fit
plt.clf()                                        # Ryd plottevinduet, i tilfælde af at du har plottet noget andet
if len(vel) < 100:                               # Hvis der ikke er så mange galakser...
    plt.scatter(dist,vel,color='b')              # ...så lav en prik for hvert datapunkt ('b' betyder blå)
else:                                            # Men hvis der er mange...
    plt.scatter(dist,vel,color='b',s=1,alpha=.1) # ...så gør prikkerne små og gennemsigtige
plt.plot(d_akse,bestfit,'r')                     # Plot fittet ('r' betyder rød)
#plt.plot(d_akse,d_akse*70,'r--')                # Plot linjen for H0 = 70 km/s/Mpc. Denne linje kan bruges til at vise, at hastigheder ved store afstande afviger fra den lineære sammenhæng
plt.xlabel('Afstand / Mpc')                      # Tekst på x-aksen
plt.ylabel('Hastighed / km/s')                   # Tekst på y-aksen
plt.title(inputfil)                              # Skriv filnavnet øverst i plottet (praktisk hvis du har flere åbne plots)


# Print resultatet. {:.1f} betyder "kommatal med ét decimal"
print('Bedste fit for Hubble-konstanten er H0 = {:.1f} ± {:.1f} km/s/Mpc'.format(H0,sigH0))


# Hvis du vil have plottet til at ligne Hubbles figur helt, med samme akselængder osv.,
# så skriv hubble1929=True i stedet for False nedenfor.
# (men husk, at hvis du plotter andre data med større afstande eller hastigheder, så falder de udenfor plottevinduet)
hubble1929 = False
if hubble1929:
    plt.xlim([-.333,2.333])
    plt.ylim([-333,1333])
    plt.xticks([0,1,2])
    plt.yticks([0,500,1000])
    plt.grid()
    plt.gca().set_axisbelow(True)

plt.show()  # Vis plottevinduet (skal stå nederst, da programmet ikke kører videre)

# Hvis du vil gemme plottet, så sæt gemplot=True nedefor.
# Så får du en .png-fil med samme navn som datafilen, bortset fra at ".dat" er skiftet ud med ".png"
gemplot = False
if gemplot:
    plt.savefig(inputfil[0:-3]+'png')
