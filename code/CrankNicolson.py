import numpy as np

from scipy import constants as c
from scipy.sparse import linalg
from scipy import sparse

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
  
  def main_sparse(self, x_array, t_array, delta_x, delta_t, V_array, psi_init_array):
    Nx = len(x_array)
    H_row = []
    H_col = []
    H_data = []
    
    I_row = []
    I_col = []
    I_data = []
    for n in range(Nx):
      if(n!=0):
        H_row.append(n)
        H_col.append(n-1)
        H_data.append((-c.hbar**2/(2*c.m_n*delta_x**2)))
      
      H_row.append(n)
      H_col.append(n)
      H_data.append((-c.hbar**2/(2*c.m_n*delta_x**2))*(-2)+V_array[n])
      
      I_row.append(n)
      I_col.append(n)
      I_data.append(1)
      
      if(n!=Nx-1):
        H_row.append(n)
        H_col.append(n+1)
        H_data.append((-c.hbar**2/(2*c.m_n*delta_x**2)))
        
    H = sparse.csc_matrix((H_data, (H_row, H_col)))
    I = sparse.csc_matrix((I_data, (I_row, I_col)))
      
    Ap = I+(1j*delta_t/(2*c.hbar))*H
    Am = I-(1j*delta_t/(2*c.hbar))*H
    
    print("end calc Hamiltonian")
    
    psi_array = np.zeros((len(t_array) + 1, len(psi_init_array)), dtype=np.complex)
    psi_array[0, :] = psi_init_array
    psi_old = psi_init_array.copy()
        
    solve_Ap = linalg.factorized(Ap)
    
    for i in range(1, len(t_array) + 1):
        print(f"{i}/{len(t_array)}")
        psi_new = solve_Ap(Am @ psi_old)
        psi_old[0] = psi_old[-1] = 0
        psi_old = psi_new
        psi_array[i, :] = psi_new
    
    print("end calc")
    return psi_array

# https://qiita.com/sci_Haru/items/e551488eddbc84818cc0
# https://qiita.com/iwasaki620/items/603220d9102e82d4438e