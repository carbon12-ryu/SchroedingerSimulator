import numpy as np


class Potential():
  def main(self, x_array, settings):
    V = np.zeros_like(x_array)
    for setting in settings:
      if(setting["type"] == "rect"):
        V+=self.rect(
          x_array, 
          setting["x_start"],
          setting["x_end"],
          setting["v_val"]
          )
    return V
    
  def rect(self, x_array, x_start, x_end, v_val):
    V = np.zeros_like(x_array)
    mask = (x_start <= x_array) & (x_array <= x_end)
    V[mask] = v_val
    return V