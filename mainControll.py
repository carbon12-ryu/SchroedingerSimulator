import json
import numpy as np

from code.initPsi import InitPsi
from code.potential import Potential
from code.CrankNicolson import CrankNicolson
from code.drawGraph import DrawGraph

class MainControll():
  def __init__(self):
    self.initPsi = InitPsi()
    self.potential = Potential()
    self.CrankNicolson = CrankNicolson()
    self.drawGraph = DrawGraph()
    
  def main(self):
    json_open = open('setting.json', 'r')
    setting = json.load(json_open)
    calc_setting = setting["calc_setting"]
    x_array = np.arange(
      calc_setting["x_start"],
      calc_setting["x_end"],
      calc_setting["delta_x"]
      )
    t_array = np.arange(
      calc_setting["t_start"],
      calc_setting["t_end"],
      calc_setting["delta_t"]
      )
    
    psi_init, n_energy = self.initPsi.main(
      x_array,
      setting["init_psi"]
      )
    V = self.potential.main(
      x_array,
      setting["potential"],
      )
    
    self.drawGraph.initGraph(
      x_array, 
      psi_init, 
      V,
      n_energy,
      setting["graph"]["init_name"],
      isShow=False
      )
    
    psi = self.CrankNicolson.main(
      x_array=x_array,
      t_array=t_array,
      delta_x=calc_setting["delta_x"],
      delta_t=calc_setting["delta_t"],
      V_array=V,
      psi_init_array=psi_init
    )
    
    self.drawGraph.animation(
      x_array=x_array,
      t_array=t_array,
      psi_array=psi,  
      V_array=V,
      n_energy=n_energy,
      filename=setting["graph"]["ani_name"],
    )
    
MainControll().main()

