

# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 09:51:33 2021

@author: ew126
"""

import dadi
import Selection
import pickle
import numpy
import pandas as pd
from sys import argv

syn_SFS = argv[1]
nonsyn_SFS = argv[2]
h = argv[3]
replicate = int(argv[4])
sfs_size = int(argv[5])
let_prop = argv[6]

###############Demography inference with dadi####################################

#the demographic function we will be using for our analysis. This describes
#a two epoch model with one historical size change.
def two_epoch(params, ns, pts):
    nu, T = params
    #add dadi.Integration.timescale_factor = 0.001
    xx = dadi.Numerics.default_grid(pts)
    
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T, nu)
    
    #calculate spectrum
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

#Demographic model with selection
def two_epoch_sel(params, ns, pts):
    nu, T, gamma = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx, gamma=gamma)
    phi = dadi.Integration.one_pop(phi, xx, T, nu, gamma=gamma)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

# Variables description:
# pneu = neutral paramter
# alpha + beta = gamma parameters
# let = lethal parameter

def neugammalet(mgamma, pneu, plet, alpha, beta):
    mgamma = -mgamma
    
    #neutral
    if ( 0 <= mgamma ) and ( mgamma < 1e-4 ) :
        return pneu/(1e-4) + (1-pneu-plet)*Selection.gamma_dist(-mgamma, alpha, beta)
    
    #gamma
    else:
        return Selection.gamma_dist(-mgamma, alpha, beta) * (1 - pneu-plet)
  
#load synonymous SFS
syn_sfs = numpy.genfromtxt(syn_SFS,delimiter="\n")
syn_sfs = dadi.Spectrum(syn_sfs)

#initial parameter guesses, convert to coalescent units
Nanc = 10000. 
N1 = 1 
T1 = 0.01 


pts_l = [1200, 1600, 2000] 
func = two_epoch
func_ex = dadi.Numerics.make_extrap_log_func(func)
params = [N1, T1] 
lower_bound = [1e-2, 0.0001]
upper_bound = [10,50.]
fixed_params = [None, None]

# randomize starting point
p0 = dadi.Misc.perturb_params(params, upper_bound=upper_bound)

print("This is a initial random point")
print(p0)

# optimize
popt = dadi.Inference.optimize_log(p0, syn_sfs, func_ex, pts_l, 
								   lower_bound=lower_bound, 
								   upper_bound=upper_bound,
								   fixed_params=fixed_params,
								   verbose=len(p0), maxiter=30)

print("This is the demographic model")

print(popt)
demo = popt


theta_s = dadi.Inference.optimal_sfs_scaling(func_ex(popt, syn_sfs.sample_sizes, pts_l),syn_sfs)

print("This is theta synonymous")
print(theta_s)


#####################Selection inference with fitdadi##########################################

#compute Nanc for DFE rescaling
mu = 1.5e-8 
Ls = 1454346 * 30 / 3.31  #sum of total exon length in slim
Nanc = theta_s/(4*mu*Ls)

#compute nonsynonymous theta
mut_ratio = 2.31 # mu_ns/mu_s, right
theta_ns = theta_s * mut_ratio

#set int breaks
int_breaks=[0.1, 1e-2, 1e-3, 1e-4, 1e-5] 

#convert to gamma scale
int_breaks = [x*2*Nanc for x in int_breaks]
s1, s2, s3, s4, s5 = int_breaks

#set eps, tiny relative to gamma -- minimize trapezoid integration error
eps = 1e-5

#custom vector of gammas
gamma_list = []
gamma_list = numpy.append(gamma_list, -numpy.logspace(numpy.log10(s1),numpy.log10(s2), 50))
gamma_list = numpy.append(gamma_list, -numpy.logspace(numpy.log10(s2-eps),numpy.log10(s3), 50))
gamma_list = numpy.append(gamma_list, -numpy.logspace(numpy.log10(s3-eps),numpy.log10(s4), 50))
gamma_list = numpy.append(gamma_list, -numpy.logspace(numpy.log10(s4-eps),numpy.log10(s5), 50))
gamma_list = numpy.append(gamma_list, -numpy.logspace(numpy.log10(s5-eps),numpy.log10(eps), 50))

#print gamma list
print("This is the gamma list")
print(gamma_list)

#parameters for building selection spectra
ns = syn_sfs.sample_sizes

#inferring SFS
spectra = Selection.spectra(popt, ns, two_epoch_sel, pts_l=pts_l, 
							gamma_list=gamma_list, echo=True, mp=False, n=Nanc) #, cpus=30)

#could see WARNING:Numeric:Extrapolation at this point

spectra.spectra = numpy.array([[0]*sfs_size, *spectra.spectra])
spectra.gammas = numpy.append(-1*2*Nanc, spectra.gammas)

print(spectra)

#pickle.dump(spectra, open('spectra_h_{coef}.sp'.format(coef=h_coefficient),'wb'))

#load nonsynonymous
nonsyn_sfs = numpy.genfromtxt(nonsyn_SFS,delimiter="\n")
nonsyn_sfs = dadi.Spectrum(nonsyn_sfs)

let_points = numpy.arange(0.0, 0.52, 0.02)

for l in let_points :
	
    neugammalet_vec = numpy.frompyfunc(neugammalet, 5, 1) #----5, 1

    # parameters = pneu, alpha, beta, let 
    # fixed pneu = 0 for gammalet model
    lower_bound=[0.0, l, 1e-3, 1e-2]         
    upper_bound=[0.0, l, 1.0, 50000.] 
    params = (0.0, l, 0.2, 1000.) 
    
    gamma_max_likelihoods = []
    gamma_guesses = dict()
    
    p0 = dadi.Misc.perturb_params(params, upper_bound=upper_bound)
    popt = Selection.optimize_log(p0, nonsyn_sfs, spectra.integrate, neugammalet_vec,
    								 theta_ns, lower_bound=lower_bound, 
    								 upper_bound=upper_bound, verbose=len(p0),
    								 fixed_params = [0.0, l, None, None], maxiter=25)
    
    gamma_max_likelihoods.append(popt[0])
    gamma_guesses[popt[0]] = popt
    
    modelsfs = spectra.integrate(popt[1], neugammalet_vec, theta_ns) 
    numpy.savetxt('log_inference/computed_SFS_{l}_{let_prop}_{h}_{ss}_{replicate}.txt'.format(l = l, let_prop = let_prop, ss=sfs_size - 2, replicate=replicate, h =h), modelsfs, delimiter=',') 
    
    results_gamma = []
    print(gamma_guesses.keys())
    for i in range(len(gamma_guesses.keys())):
        best_popt_gamma = gamma_guesses[gamma_max_likelihoods[i]]
        results_gamma.append([best_popt_gamma[0], best_popt_gamma[1][0], best_popt_gamma[1][1], best_popt_gamma[1][2], best_popt_gamma[1][3]])
    
    df_gamma = pd.DataFrame(results_gamma, columns =['Likelihood', 'pneu', 'plet','alpha','beta'], dtype = float) #, 'alpha', 'beta'
    df_gamma['Nanc'] = str(Nanc) 
    df_gamma['nu'] = str(demo[0])
    df_gamma['T'] = str(demo[1])
    df_gamma['scale'] = 'na'
    df_gamma.to_csv('log_inference/DFE_inference_{l}_{let_prop}_{h}_{ss}_{replicate}.csv'.format(l = l, let_prop = let_prop, ss = sfs_size - 2, h=h,replicate=replicate), index=False)
    


