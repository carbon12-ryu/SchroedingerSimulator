import numpy as np

from scipy import constants as c


class InitPsi():
  def main(self, x_array, setting):
    n_energy = c.h**2/(2*c.m_n*setting["wavelength"]**2)
    if(setting["type"]=="gauss"):
      psi_init_array = self.gauss_wave(
        x_array,
        setting["sigma"],
        setting["center"],
        setting["wavelength"]
        )
    return psi_init_array, n_energy
  
  def gauss_wave(self, x_array, sigma, center, wavelength):
    k = 2*np.pi/wavelength
    psi_init_array = (
      np.exp(1j*k*x_array)
      *np.exp(-(x_array-center)**2/sigma**2)/(np.sqrt(2*np.pi)*sigma)
      )
    return psi_init_array