import numpy as np

from scipy import constants as c
from matplotlib import pyplot as plt
from matplotlib import animation

class DrawGraph():
  def initGraph(self, x_array, psi_init, V_array, n_energy, filename, isShow=False):
    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(111)
    ax1.plot(x_array, np.abs(psi_init) ,color='black', label=r'$|\psi(x)|$')
    ax1.plot(x_array, np.real(psi_init),color='blue', alpha=0.5, label=r'$Re{\psi(x)}$')
    ax1.plot(x_array, np.imag(psi_init),color='red', alpha=0.5, label=r'$Im{\psi(x)}$')

    ax2 = ax1.twinx()
    ax2.fill_between(x_array, V_array, color='gray', alpha=0.3)
    ax2.axhline(n_energy, color='black', linestyle='--', label='neutron energy (free)')

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, loc='upper right')

    ax1.set_xlabel('position (m)')
    ax1.set_ylabel('probability density')
    ax1.grid(True)
    ax2.set_ylabel('potential')
    
    fig.savefig(filename)
    if(isShow):
      plt.show()
    plt.close()
    
  def animation(self, x_array, t_array, psi_array, V_array, n_energy, filename):
    x_array = x_array*1e10
    V_array = V_array/c.e*1e3
    n_energy = n_energy/c.e*1e3
    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(111)

    ax2 = ax1.twinx()
    ax2.fill_between(x_array, V_array, color='gray', alpha=0.3)
    ax2.axhline(n_energy, color='black', linestyle='--', label='neutron energy (free)')
    ax1.set_ylim(np.min(np.real(psi_array))*1.1, np.max(np.real(psi_array))*1.1)
    ax2.set_ylim(0, np.max((np.max(V_array), n_energy)) * 1.1)
    
    
    lines = ax1.plot([], [], color='black', label=r'$|\psi(x)|$')
    lines2 = ax1.plot([], [], color='blue', alpha=0.5, label=r'$Re\{\psi(x)\}$')
    lines3 = ax1.plot([], [], color='red', alpha=0.5, label=r'$Im\{\psi(x)\}$')
    text = ax1.text(0.05, 0.95, '', transform=ax1.transAxes, fontsize=12, verticalalignment='top')
    line1 = lines[0]
    line2 = lines2[0]
    line3 = lines3[0]
    
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, loc='upper right')

    ax1.set_xlabel(r'position ($\AA$)')
    ax1.set_ylabel('sqrt of probability density')
    ax2.set_ylabel('potential(meV)')
    ax1.grid(True)
    
    def update(i):
      line1.set_data(x_array, np.abs(psi_array[i]))
      line2.set_data(x_array, np.real(psi_array[i]))
      line3.set_data(x_array, np.imag(psi_array[i]))
      text.set_text(f't={t_array[i]*1e9:.3f}ns')
      return [line1, line2, line3, text]
      
    ani = animation.FuncAnimation(
      fig,
      update,
      frames=range(0, len(t_array), 30),
      interval=1,
      blit=True
      )
    
    ani.save(filename, writer='pillow')
    plt.close()