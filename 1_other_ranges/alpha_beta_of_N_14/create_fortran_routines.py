import numpy as np
import matplotlib.pyplot as pl
from matplotlib import patches
from scipy.optimize import curve_fit, root, fsolve
from scipy.signal import argrelextrema
from numpy import dot, outer
from os import listdir
from os.path import isfile, join

def main():
    
    # Set parameters
    n_minimax = 14                     # Number of minimax points
    R_minimax = int(146)              # Range of the minimax approximation
    n_x       = 12000                   # total number of points on the x-axis for optimization
    eps_diff  = 10**(-10)

    ranges = [50,70,100,150,300,600,1000,1500,3000,10000,30000]

    f = open('minimax_grid.F', 'w')
    print('SUBROUTINE get_coeff_'+str(n_minimax)+'(k,E_range, aw)', file=f) 
    print('INTEGER, INTENT(IN)                                :: k', file=f) 
    print('REAL(KIND=dp), INTENT(IN)                          :: E_range', file=f) 
    print('REAL(KIND=dp), DIMENSION(2*k), INTENT(OUT)         :: aw', file=f) 
    print('', file=f) 

    for R_minimax in ranges:

      alphas_betas = find_file(n_minimax, R_minimax)


      if R_minimax == ranges[0]:
         print('IF(E_range < '+str(R_minimax)+'.0_dp) THEN', file=f) 
      elif R_minimax == ranges[-1]:
         print('ELSE', file=f) 
      else:
         print('ELSE IF(E_range < '+str(R_minimax)+'.0_dp) THEN', file=f) 

      print('aw(:) = (/ &', file=f)

      for index, alpha_beta in enumerate(alphas_betas):
         if index == np.size(alphas_betas)-1:
            print(str(alpha_beta)+"_dp /)", file=f) 
         else:
            print(str(alpha_beta)+"_dp, &", file=f) 

    print('END IF', file=f) 
    print('END SUBROUTINE', file=f) 

    f.close

def find_file(n_minimax, R_minimax):

    files = [f for f in listdir(".") if isfile(join(".", f))]

    for f in files:
        if(not f.startswith("alpha")): continue
        R_file = int(f[21:34])
        if( R_file == R_minimax ): 
            with open(f) as filetoread:
                alphas_betas = [np.float64(x) for x in filetoread]

    return alphas_betas


if __name__ == "__main__":
    main()

