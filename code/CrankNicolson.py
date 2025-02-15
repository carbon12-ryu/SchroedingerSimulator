import numpy as np

from scipy import constants as c

class CrankNicolson():  
  def krondel(self, nn,mm):
    if nn==mm:
        return 1
    else:
        return 0
  
  def main(self, x_array, t_array, delta_x, delta_t, V_array, psi_init_array):
    Nx = len(x_array)
    I_mat=np.identity(Nx, dtype=float)
    Hmat=np.empty([Nx, Nx])  
    for n in range(Nx):
        for m in range(Nx):
            Hmat[n,m] = (-c.hbar**2/(2*c.m_n*delta_x**2))*(self.krondel(n,m-1)-2*self.krondel(n,m)+self.krondel(n,m+1))+V_array[n]*self.krondel(n,m)
    Ap = I_mat+(1j*delta_t/(2*c.hbar))*Hmat
    Am = I_mat-(1j*delta_t/(2*c.hbar))*Hmat
    Ainv=np.linalg.inv(Ap)
    
    AA = np.dot(Ainv, Am)
    print("end calc Hamiltonian")
    psi_array = np.zeros((len(t_array) + 1, len(psi_init_array)), dtype=np.complex)
    psi_array[0, :] = psi_init_array
    psi_old = psi_init_array.copy()

    for i in range(1, len(t_array) + 1):
        print(f"{i}/{len(t_array)}")
        psi_old[0] = 0
        psi_old[-1] = 0
        psi_new = np.dot(AA, psi_old)
        psi_array[i, :] = psi_new
        psi_old = psi_new
    
    print("end calc")
    return psi_array

# https://qiita.com/sci_Haru/items/e551488eddbc84818cc0