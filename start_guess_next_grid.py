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

    min_alpha, max_alpha, min_beta, max_beta, alphas_fac, betas_fac = [], [], [], [], [], []
    for i_minimax in minimax_range:
       alphas_betas = np.loadtxt("alpha_beta_of_N_"+str(i_minimax))
       min_alpha.append(alphas_betas[0])
       max_alpha.append(alphas_betas[i_minimax-1])
       min_beta.append(alphas_betas[i_minimax])
       max_beta.append(alphas_betas[2*i_minimax-1])

       alphas_fac.append(np.array(alphas_betas[1:i_minimax])/np.array(alphas_betas[0:i_minimax-1]))
       betas_fac.append(np.array(alphas_betas[i_minimax+1:])/np.array(alphas_betas[i_minimax:2*i_minimax-1]))
       print("")
       print("================")
       print("i_minimax", i_minimax)
       print("================")
       print("")
       print("alphas_fac", np.array(alphas_betas[1:i_minimax])/np.array(alphas_betas[0:i_minimax-1]))
       print("")
       print("alpha_max/alpha_min", alphas_betas[i_minimax-1]/alphas_betas[0])
       print("")
#       print("betas_fac", np.array(alphas_betas[i_minimax+1:])/np.array(alphas_betas[i_minimax:2*i_minimax-1]))
#       print("")
#       print("beta_max/beta_min", alphas_betas[2*i_minimax-1]/alphas_betas[i_minimax])
#       print("")




    min_alpha_factor = np.array(min_alpha[1:])/np.array(min_alpha[0:n_previous-1])
    max_alpha_factor = np.array(max_alpha[1:])/np.array(max_alpha[0:n_previous-1])
    min_beta_factor = np.array(min_beta[1:])/np.array(min_beta[0:n_previous-1])
    max_beta_factor = np.array(max_beta[1:])/np.array(max_beta[0:n_previous-1])


    print("")
    print("min_alpha", min_alpha)
    print("max_alpha", max_alpha)
#    print("min_beta", min_beta)
#    print("max_beta", max_beta)
    print("")
    print("min_alpha_factor", min_alpha_factor)
    print("max_alpha_factor", max_alpha_factor)
#    print("min_beta_factor", min_beta_factor)
#    print("max_beta_factor", max_beta_factor)

    min_alpha_extra_factor = 2*min_alpha_factor[n_previous-2] - min_alpha_factor[n_previous-3]
    max_alpha_extra_factor = 2*max_alpha_factor[n_previous-2] - max_alpha_factor[n_previous-3]
    min_beta_extra_factor  = 2*min_beta_factor[n_previous-2] - min_beta_factor[n_previous-3]
    max_beta_extra_factor   = 2*max_beta_factor[n_previous-2] - max_beta_factor[n_previous-3]

    print("")
    print("min_alpha_extra_factor", min_alpha_extra_factor)
    print("max_alpha_extra_factor", max_alpha_extra_factor)
#    print("min_beta_extra_factor", min_beta_extra_factor)
#    print("max_beta_extra_factor", max_beta_extra_factor)
    print("")
    print("")



#####################################################
# 1. Extrapolate the nodes alpha                    #
#####################################################

    upper_alphas_fac_extrapol = np.array(alphas_betas[5:i_minimax])/np.array(alphas_betas[4:i_minimax-1])
    n_lower = 5
    lower_alphas_fac_extrapol = 2*alphas_fac[n_previous-1][0:n_lower+1] - alphas_fac[n_previous-2][0:n_lower+1]
    n_rest = n_minimax - np.size(upper_alphas_fac_extrapol) - np.size(lower_alphas_fac_extrapol) - 1
    middle_alphas_fac_extrapol = np.linspace(lower_alphas_fac_extrapol[n_lower], upper_alphas_fac_extrapol[0], n_rest+2)

    alphas_fac_extrapol = np.append(lower_alphas_fac_extrapol, np.append(middle_alphas_fac_extrapol[1:n_rest+1], upper_alphas_fac_extrapol))

    print("")
    print("")
    print("alpha_fac_extrapol =", alphas_fac_extrapol)
    print("")
    print("")

    alphas_betas_extrapol = []

#####################################################
# 2. Extrapolate the weights beta                   #
#####################################################

    upper_betas_fac_extrapol = np.array(alphas_betas[i_minimax+5:])/np.array(alphas_betas[i_minimax+4:2*i_minimax-1])
    n_lower = 2
    lower_betas_fac_extrapol = 2*betas_fac[n_previous-1][0:n_lower+1] - betas_fac[n_previous-2][0:n_lower+1]
    n_rest = n_minimax - np.size(upper_betas_fac_extrapol) - np.size(lower_betas_fac_extrapol) - 1
    middle_betas_fac_extrapol = np.linspace(lower_betas_fac_extrapol[n_lower], upper_betas_fac_extrapol[0], n_rest+2)

    betas_fac_extrapol = np.append(lower_betas_fac_extrapol, np.append(middle_betas_fac_extrapol[1:n_rest+1], upper_betas_fac_extrapol))

    missing_factor = np.prod(betas_fac[n_previous-1])/np.prod(betas_fac_extrapol)*max_beta_extra_factor/min_beta_extra_factor

    multiplication_factor = missing_factor**(2/(n_minimax/2*(n_minimax/2+1)))

    factor = []
    for n in range(n_minimax//2):
      factor.append(multiplication_factor**(n_minimax//2-n))

    factor = np.append(np.array(factor), np.ones(n_minimax//2-1))

    betas_fac_extrapol = betas_fac_extrapol * factor

    missing_factor_new = np.prod(betas_fac[n_previous-1])/np.prod(betas_fac_extrapol)*max_beta_extra_factor/min_beta_extra_factor

    alphas_betas_extrapol.append(min_alpha_extra_factor*min_beta[n_previous-1])
    for i_beta in range(1,n_minimax):
      alphas_betas_extrapol.append(alphas_betas_extrapol[i_beta-1]*betas_fac_extrapol[i_beta-1])

    print("")
    print("")
    print("betas =", alphas_betas_extrapol)
    print("")
    print("")

if __name__ == "__main__":
    main()

