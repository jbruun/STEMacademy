import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
import scipy.integrate as spi

# Constants
H0 = 70.0  # Hubble constant in km/s/Mpc
H0_SI = H0 * 1000 / (3.0857e22)  # Convert to SI units (1/s)

def compute_scale_factor(Omega_m, Omega_Lambda):
    """For open universes numerically solve for the scale factor a(t) over cosmic time, including the collapse. For a closed universe use parametric solution."""
    Omega_m = round(Omega_m,2)
    Omega_Lambda = round(Omega_Lambda,2)
    if Omega_m>1 and Omega_Lambda==0:
       """Parametric solution for a closed universe with Omega_m > 1 and Omega_Lambda = 0."""
       # Define the range of the parameter theta (from 0 to 2*pi)
       theta = np.linspace(0, 2 * np.pi, 5000)

       # Parametric solution for a(t) and t(theta)
       a = (Omega_m / (2 * (Omega_m - 1))) * (1 - np.cos(theta))
       t = (Omega_m / (2 * (Omega_m - 1)**(3/2))) * (theta - np.sin(theta))
       return t, a  # Return time in Gyr and scale factor a(t)
    else:
       def dtda(a):
          Omega_k = 1. - Omega_m - Omega_Lambda
          return 1./np.sqrt(Omega_m/a+Omega_Lambda*a**2+Omega_k)

       a = 10**(np.arange(0,130,1)/40.-2.)
       t = [0 for x in range(130)]
       for n in range(130):
          I = spi.quad(dtda,np.min(a),a[n])
          t[n] = I[0]
       return t, a  # Return time in Gyr and scale factor


def plot_friedmann(ax1, ax2, Omega_m, Omega_Lambda):
    """Update the bottom plot with the scale factor evolution for selected Omega_m and Omega_Lambda."""
    # Compute the scale factor evolution
    t, a = compute_scale_factor(Omega_m, Omega_Lambda)

   # Clear the previous plot and update with new data
    ax2.cla()
    ax2.plot(t, a, label=f"Ω_m={Omega_m:.2f}, Ω_Λ={Omega_Lambda:.2f}")
    ax2.set_title("Scale Factor Evolution")
    ax2.set_xlabel("Cosmic Time (Hubble time)")
    ax2.set_ylabel("Scale Factor a(t)")
    ax2.legend()
    ax2.grid(True)

    # Redraw the figure
    plt.draw()

def on_click(event, ax1, ax2):
    """Handle mouse click event on the top plot to select Omega_m and Omega_Lambda."""
    if event.inaxes == ax1:
        Omega_Lambda = event.xdata
        Omega_m = event.ydata
        print(f"Selected Ω_m={Omega_m:.2f}, Ω_Λ={Omega_Lambda:.2f}")
        
        # Update the bottom plot with the scale factor for the selected Omegas
        plot_friedmann(ax1, ax2, Omega_m, Omega_Lambda)

def main():
    # Create the figure and two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))
    
    # Top plot: Omega_Lambda vs Omega_matter (user clicks here)
    Omega_Lambda_values = np.linspace(0, 2.0, 100)
    Omega_m_values = np.linspace(0, 2.0, 100)  # Extended to 2.0
    X, Y = np.meshgrid(Omega_Lambda_values, Omega_m_values)
    Omega_k = 1.0 - X - Y
    
    ax1.contourf(X, Y, Omega_k, levels=np.linspace(-1, 1, 20), cmap='RdYlBu')
    ax1.set_title("Click to Select Ω_m and Ω_Λ")
    ax1.set_xlabel("Ω_Λ (Omega Lambda)")
    ax1.set_ylabel("Ω_m (Omega Matter)")
    ax1.grid(True)
    
    # Bottom plot: Scale factor evolution (initially blank)
    ax2.set_title("Scale Factor Evolution")
    ax2.set_xlabel("Cosmic Time (Gyr)")
    ax2.set_ylabel("Scale Factor a(t)")
    ax2.grid(True)
    
    # Interactive cursor and click events
    cursor = Cursor(ax1, useblit=True, color='black', linewidth=1)
    fig.canvas.mpl_connect('button_press_event', lambda event: on_click(event, ax1, ax2))

    # Show the plot
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()

