import numpy as np
import matplotlib.pyplot as pl
from matplotlib import patches
from scipy.optimize import curve_fit, root, fsolve
from scipy.signal import argrelextrema
from numpy import dot, outer

def main():
    
    # Set parameters
    n_minimax = 22                     # Number of minimax points
    R_minimax = 10**10                 # Range of the minimax approximation
    n_x       = 1000                   # total number of points on the x-axis for optimization
    eps_diff  = 10**(-10)

    n_previous = (n_minimax-10)//2

    minimax_range = (np.linspace(10, n_minimax-2, n_previous)).astype(int)

    print("range", minimax_range)

    min_alpha, max_alpha, min_beta, max_beta = [], [], [], []
    for i_minimax in minimax_range:
       alphas_betas = np.loadtxt("alpha_beta_of_N_"+str(i_minimax))
       min_alpha.append(alphas_betas[0])
       max_alpha.append(alphas_betas[i_minimax-1])
       min_beta.append(alphas_betas[i_minimax])
       max_beta.append(alphas_betas[2*i_minimax-1])

       print("")
       print("================")
       print("i_minimax", i_minimax)
       print("================")
       print("")
       print("alphas_fac", np.array(alphas_betas[1:i_minimax])/np.array(alphas_betas[0:i_minimax-1]))
       print("")
       print("betas_fac", np.array(alphas_betas[i_minimax+1:])/np.array(alphas_betas[i_minimax:2*i_minimax-1]))


    min_alpha_factor = np.array(min_alpha[1:])/np.array(min_alpha[0:n_previous-1])

    print("")
    print("min_alpha", min_alpha)
    print("max_alpha", max_alpha)
    print("min_beta", min_beta)
    print("max_beta", max_beta)
    print("")
    print("min_alpha_factor", min_alpha_factor)
    print("max_alpha_factor", np.array(max_alpha[1:])/np.array(max_alpha[0:n_previous-1]))
    print("min_beta_factor", np.array(min_beta[1:])/np.array(min_beta[0:n_previous-1]))
    print("max_beta_factor", np.array(max_beta[1:])/np.array(max_beta[0:n_previous-1]))

    min_alpha_extra_factor = 2*min_alpha_factor[n_previous-2] - min_alpha_factor[n_previous-3]

    print("")
    print("min_alpha_extra_factor", min_alpha_extra_factor)
    print("")
    print("")




if __name__ == "__main__":
    main()

