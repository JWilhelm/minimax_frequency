import numpy as np
import matplotlib.pyplot as pl
from matplotlib import patches
from scipy.optimize import curve_fit, root, fsolve
from scipy.signal import argrelextrema
from numpy import dot, outer

def main():
    
    # Set parameters
    n_minimax = 38                     # Number of minimax points
    R_minimax = 10**13                 # Range of the minimax approximation
    n_x       = 4000                   # total number of points on the x-axis for optimization
    eps_diff  = 10**(-10)

    xdata = 10**(np.logspace(0,np.log10(np.log10(R_minimax*10)),n_x))/10
    ydata = np.zeros(n_x)

    alphas_betas_init = np.loadtxt("alpha_beta_of_N_"+str(n_minimax))

#    print("xdata", xdata)

    alphas_betas_L2_opt, alphas_betas_conv = curve_fit(eta, xdata, ydata, p0=alphas_betas_init)

    alphas_betas_E = np.append(alphas_betas_L2_opt,1)

    E_old = alphas_betas_E[-1]*2

    sort_indices = np.argsort(alphas_betas_E[0:n_minimax])
#    np.savetxt("alpha_beta_of_N_"+str(n_minimax)+"_L2", np.append(alphas_betas_E[sort_indices],alphas_betas_E[sort_indices+n_minimax]) )
    i = 0
    while (alphas_betas_E[-1]/E_old < 1-eps_diff or alphas_betas_E[-1] > E_old):

        E_old = alphas_betas_E[-1]
        extrema_x = np.append(xdata[0], xdata[argrelextrema(eta_plotting(xdata,alphas_betas_E[0:np.size(alphas_betas_E)-1]), np.greater)[0]])
        extrema_x = np.append(extrema_x, xdata[argrelextrema(eta_plotting(xdata,alphas_betas_E[0:np.size(alphas_betas_E)-1]), np.less)[0]])
        print("number of extrema =", np.size(extrema_x))
        alphas_betas_E[np.size(alphas_betas_E)-1] = np.average(np.abs(eta_plotting(extrema_x,alphas_betas_E[0:np.size(alphas_betas_E)-1])))
        i += 1
        alphas_betas_E = fsolve(eta_for_alphas_betas_E_update, x0=alphas_betas_E, args=extrema_x)
        print("iteration =", i, "E =",  alphas_betas_E[-1])

    sort_indices = np.argsort(alphas_betas_E[0:n_minimax])
    np.savetxt("alpha_beta_of_N_"+str(n_minimax), np.append(alphas_betas_E[sort_indices],alphas_betas_E[sort_indices+n_minimax]) )

    fig1, (axis1) = pl.subplots(1,1)
    axis1.set_xlim((0.8,R_minimax))
    axis1.semilogx(xdata,eta_plotting(xdata,alphas_betas_L2_opt))
    axis1.semilogx(xdata,eta_plotting(xdata,alphas_betas_E))
    axis1.semilogx([0.8,R_minimax], [alphas_betas_E[-1],alphas_betas_E[-1]])
    axis1.semilogx([0.8,R_minimax], [-alphas_betas_E[-1],-alphas_betas_E[-1]])

    pl.show()

def build_denominator_matrix(xdata,alphas):
    matrix = []
    for x in xdata:
      matrix.append((2*x/(x**2+np.array(alphas)**2))**2/np.pi)
    return matrix

def eta(x, *params):
#    print("params", params)
    denominator_matrix = build_denominator_matrix(x, np.array(params[0:np.size(params)//2]))
    return 1/x - np.array(denominator_matrix).dot(params[np.size(params)//2:])

def eta_plotting(x, *params):
    params_1d = np.transpose(params)[:,0]
    denominator_matrix = build_denominator_matrix(x, np.array(params_1d[0:np.size(params)//2]))
    return 1/x - np.array(denominator_matrix).dot(params_1d[np.size(params)//2:(np.size(params)//2)*2])

def eta_for_alphas_betas_E_update(x, *params):
    params_1d = np.transpose(params)[:,0]
    size_params = np.size(params_1d)
    size_x = np.size(x)
    E = np.empty(size_x)
    E[0:size_x//2+1] = x[size_params-1]
    E[size_x//2+1:] = -x[size_params-1]
    denominator_matrix = build_denominator_matrix(params_1d,np.array(x[0:np.size(x)//2]))
    return 1/params_1d - np.array(denominator_matrix).dot(x[np.size(x)//2:(np.size(x)//2)*2]) - E

if __name__ == "__main__":
    main()

