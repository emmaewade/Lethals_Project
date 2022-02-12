
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 09:51:33 2021

@author: ew126
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 18:37:45 2021

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

###############Demography stuff####################################

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

#Demographic model with selection, from Bernard's
def two_epoch_sel(params, ns, pts):
    nu, T, gamma = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx, gamma=gamma)
    phi = dadi.Integration.one_pop(phi, xx, T, nu, gamma=gamma)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

#pneu = neutal paramter
#alpha + beta = gamma parameter
#let = lethal parameter

# < 1e-4 = pneu
# 1e-4 < gamma < (1e-2????)
#so, 1e-1 = let

def neugammalet(mgamma, pneu,plet, alpha, beta):
    mgamma = -mgamma
    
    #neutral
    if ( 0 <= mgamma ) and ( mgamma < 1e-4 ) :
        return pneu/(1e-4) + (1-pneu-plet)*Selection.gamma_dist(-mgamma, alpha, beta)
    
    #gamma
    else:
        return Selection.gamma_dist(-mgamma, alpha, beta) * (1 - pneu-plet)
    
def neuexpolet(mexpo, pneu, plet, scale):
    mexpo = -mexpo
    
    #neutral
    if ( 0 <= mexpo ) and ( mexpo < 1e-4 ) :
        return pneu/(1e-4) + (1-pneu-plet)*Selection.exponential_dist(-mexpo, scale)
    
    #expo
    else:
        return Selection.exponential_dist(-mexpo,scale) * (1 - pneu-plet) 

#load synonymous SFS
#infilename = syn #######MAKE################
syn_sfs = numpy.genfromtxt(syn_SFS,delimiter="\n")
syn_sfs = dadi.Spectrum(syn_sfs)

#initial parameter guesses, convert to coalescent units
Nanc = 10000. ###added
N1 = 1 #means the population size doesn't change???
T1 = 0.01 #right???

#setup dadi stuff
pts_l = [1200, 1600, 2000] #bernard has [600, 800, 1000], kept Isabel's because Nanc was closer
func = two_epoch
func_ex = dadi.Numerics.make_extrap_log_func(func)
params = [N1, T1] #N2, T2, NC, TC] #changed
lower_bound = [1e-2, 0.0001]#, 1e-2, 1e-3, 1, 1e-3]
upper_bound = [10,50.]#, 10, 0.5, 200, 0.5] ###########bounds for T1???##################
fixed_params = [None, None]#, None, None, None, None] #made none

# fit demographic model
# need to make sure parameters are a good fit
# i usually run 25-30 runs and make sure they've converged

# randomize starting point
p0 = dadi.Misc.perturb_params(params, upper_bound=upper_bound)
#took out p0[2] b/c t isn't fixed, not bottleneck

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


#####################Selection stuff##########################################

#compute Nanc for DFE rescaling
mu = 1.5e-8 #if mutation rate, correct
Ls = 1454346 * 30 / 3.31 #8058343 # length ##found by adding total exon length in slim, difference everytime though
Nanc = theta_s/(4*mu*Ls)

#compute nonsynonymous theta
mut_ratio = 2.31 # mu_ns/mu_s, right
theta_ns = theta_s * mut_ratio


#set int breaks, Bernard does this part differently
int_breaks=[0.1, 1e-2, 1e-3, 1e-4, 1e-5] 

#convert convert to gamma scale
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
#print(gamma_list)


#parameters for building selection spectra
ns = syn_sfs.sample_sizes

#**copy Selection.py before running
#inferring SFS
spectra = Selection.spectra(popt, ns, two_epoch_sel, pts_l=pts_l, 
							gamma_list=gamma_list, echo=True, mp=False, n=Nanc) #, cpus=30)

#could see WARNING:Numeric:Extrapolation at this point

#manual add lethal bound (|s| = 1) <------ commented out on Isabel's right???
spectra.spectra = numpy.array([[0]*sfs_size, *spectra.spectra]) #*******************
spectra.gammas = numpy.append(-1*2*Nanc, spectra.gammas)

#print spectra, pickle for later?
#print(spectra)

#pickle.dump(spectra, open('spectra_h_{coef}.sp'.format(coef=h_coefficient),'wb'))

#load nonsynonymous
#infilename = nonsyn
nonsyn_sfs = numpy.genfromtxt(nonsyn_SFS,delimiter="\n")
nonsyn_sfs = dadi.Spectrum(nonsyn_sfs)
#nonsyn_sfs = nonsyn_sfs.fold() #don't fold??

 #fit a gamma DFE, Bernard's use a neutral+gamma DFE
dfe=Selection.gamma_dist 
###parameteres = pneu, alpha, beta, let <- #added parameter for neutrals, lethals
lower_bound=[1e-3, 1]         
upper_bound=[100, 100000] 
params = (0.2, 1000) 
gamma_max_likelihoods = []
gamma_guesses = dict()

p0 = dadi.Misc.perturb_params(params, upper_bound=upper_bound)
popt = Selection.optimize_log(p0, nonsyn_sfs, spectra.integrate, dfe,
								 theta_ns, lower_bound=lower_bound, 
								 upper_bound=upper_bound, verbose=len(p0),
								 maxiter=25)
print("Selection")
print(popt)      

gamma_max_likelihoods.append(popt[0])
gamma_guesses[popt[0]] = popt                                                            

modelsfs = spectra.integrate(popt[1], dfe, theta_ns)  ###just added
#numpy.savetxt(model_file, modelsfs,
      #delimiter = ",")
numpy.savetxt('inference/computed_SFS_gamma_{let_prop}_{h}_{ss}_{replicate}.txt'.format(let_prop=let_prop, ss=sfs_size - 2, replicate=replicate, h=h), modelsfs, delimiter=',') 

results_gamma = []
print(gamma_guesses.keys())
for i in range(len(gamma_guesses.keys())):
    best_popt_gamma = gamma_guesses[gamma_max_likelihoods[i]]
    results_gamma.append([best_popt_gamma[0], best_popt_gamma[1][0], best_popt_gamma[1][1]])

df_gamma = pd.DataFrame(results_gamma, columns =['Likelihood','alpha','beta'], dtype = float) #, 'alpha', 'beta'
df_gamma['Nanc'] = str(Nanc) 
df_gamma['nu'] = str(demo[0])
df_gamma['T'] = str(demo[1])
df_gamma['plet'] = 'na'
df_gamma['pneu'] = 'na'
df_gamma['scale'] = 'na'
# Saving all the other parameters to a file
df_gamma.to_csv('inference/DFE_inference_gamma_{let_prop}_{h}_{ss}_{replicate}.csv'.format(let_prop = let_prop, ss = sfs_size - 2, h=h,replicate=replicate), index=False)

#fit a gamma DFE, Bernard's use a neutral+gamma DFE
###dfe=Selection.gamma_dist -------> changed!!!!!!!!!!!!
neugammalet_vec = numpy.frompyfunc(neugammalet, 5, 1) #----5, 1
###parameteres = pneu, alpha, beta, let <- #added parameter for neutrals, lethals
# pneu,plet, alpha, beta
lower_bound=[1e-3, 1e-3, 1e-3, 1e-2]         
upper_bound=[1, 1, 1., 50000] 
params = (.2, 0.2, 0.2, 10000.) 
gamma_max_likelihoods = []
gamma_guesses = dict()

p0 = dadi.Misc.perturb_params(params, upper_bound=upper_bound)
popt = Selection.optimize_log(p0, nonsyn_sfs, spectra.integrate, neugammalet_vec,
								 theta_ns, lower_bound=lower_bound, 
								 upper_bound=upper_bound, verbose=len(p0),
								 maxiter=25)

gamma_max_likelihoods.append(popt[0])
gamma_guesses[popt[0]] = popt

modelsfs = spectra.integrate(popt[1], neugammalet_vec, theta_ns)  ###just added
#numpy.savetxt(model_file, modelsfs,
      #delimiter = ",")
numpy.savetxt('inference/computed_SFS_neugammalet_{let_prop}_{h}_{ss}_{replicate}.txt'.format(let_prop = let_prop, ss=sfs_size - 2, replicate=replicate, h =h), modelsfs, delimiter=',') 

results_gamma = []
print(gamma_guesses.keys())
for i in range(len(gamma_guesses.keys())):
    best_popt_gamma = gamma_guesses[gamma_max_likelihoods[i]]
    results_gamma.append([best_popt_gamma[0], best_popt_gamma[1][0], best_popt_gamma[1][1], best_popt_gamma[1][2], best_popt_gamma[1][3]])

df_gamma = pd.DataFrame(results_gamma, columns =['Likelihood', 'pneu', 'plet', 'alpha','beta'], dtype = float) #, 'alpha', 'beta'
df_gamma['Nanc'] = str(Nanc) 
df_gamma['nu'] = str(demo[0])
df_gamma['T'] = str(demo[1])
df_gamma['scale'] = 'na'
# Saving all the other parameters to a file
df_gamma.to_csv('inference/DFE_inference_neugammalet_{let_prop}_{h}_{ss}_{replicate}.csv'.format(let_prop = let_prop, ss = sfs_size - 2, h=h,replicate=replicate), index=False)

neugamma_vec = numpy.frompyfunc(neugammalet, 5, 1) #----5, 1
###parameteres = pneu, alpha, beta, let <- #added parameter for neutrals, lethals
# pneu,plet, alpha, beta
lower_bound=[1e-3, 0.0, 1e-3, 1e-2]         
upper_bound=[1.0, 0.0, 1.0, 50000.] 
params = (.2, 0.0, 0.2, 10000.)  
gamma_max_likelihoods = []
gamma_guesses = dict()

p0 = dadi.Misc.perturb_params(params, upper_bound=upper_bound)
popt = Selection.optimize_log(p0, nonsyn_sfs, spectra.integrate, neugamma_vec,
								 theta_ns, lower_bound=lower_bound, 
								 upper_bound=upper_bound, verbose=len(p0), fixed_params = [None, 0.0, None, None],
								 maxiter=25)

gamma_max_likelihoods.append(popt[0])
gamma_guesses[popt[0]] = popt

modelsfs = spectra.integrate(popt[1], neugamma_vec, theta_ns)  ###just added
#numpy.savetxt(model_file, modelsfs,
      #delimiter = ",")
numpy.savetxt('inference/computed_SFS_neugamma_{let_prop}_{h}_{ss}_{replicate}.txt'.format(let_prop = let_prop, ss=sfs_size - 2, replicate=replicate, h=h), modelsfs, delimiter=',') 

results_gamma = []
print(gamma_guesses.keys())
for i in range(len(gamma_guesses.keys())):
    best_popt_gamma = gamma_guesses[gamma_max_likelihoods[i]]
    results_gamma.append([best_popt_gamma[0], best_popt_gamma[1][0], best_popt_gamma[1][1], best_popt_gamma[1][2], best_popt_gamma[1][3]])

df_gamma = pd.DataFrame(results_gamma, columns =['Likelihood', 'pneu', 'plet', 'alpha','beta'], dtype = float) #, 'alpha', 'beta'
df_gamma['Nanc'] = str(Nanc) 
df_gamma['nu'] = str(demo[0])
df_gamma['T'] = str(demo[1])
df_gamma['scale'] = 'na'
# Saving all the other parameters to a file
df_gamma.to_csv('inference/DFE_inference_neugamma_{let_prop}_{h}_{ss}_{replicate}.csv'.format(let_prop = let_prop, ss = sfs_size - 2, h=h,replicate=replicate), index=False)

gammalet_vec = numpy.frompyfunc(neugammalet, 5, 1) #----5, 1
# pneu,plet, alpha, beta
lower_bound=[0, 1e-3, 1e-3, 1e-2]         
upper_bound=[0, 1, 1., 50000] 
params = (0, 0.2, 0.2, 10000.) 
gamma_max_likelihoods = []
gamma_guesses = dict()

p0 = dadi.Misc.perturb_params(params, upper_bound=upper_bound)
popt = Selection.optimize_log(p0, nonsyn_sfs, spectra.integrate, gammalet_vec,
								 theta_ns, lower_bound=lower_bound, 
								 upper_bound=upper_bound, verbose=len(p0), fixed_params = [0.0, None, None, None],
								 maxiter=25)

gamma_max_likelihoods.append(popt[0])
gamma_guesses[popt[0]] = popt

modelsfs = spectra.integrate(popt[1], gammalet_vec, theta_ns)  ###just added
#numpy.savetxt(model_file, modelsfs,
      #delimiter = ",")
numpy.savetxt('inference/computed_SFS_gammalet_{let_prop}_{h}_{ss}_{replicate}.txt'.format(let_prop = let_prop, ss=sfs_size - 2, h=h, replicate=replicate), modelsfs, delimiter=',') 

results_gamma = []
print(gamma_guesses.keys())
for i in range(len(gamma_guesses.keys())):
    best_popt_gamma = gamma_guesses[gamma_max_likelihoods[i]]
    results_gamma.append([best_popt_gamma[0], best_popt_gamma[1][0], best_popt_gamma[1][1], best_popt_gamma[1][2], best_popt_gamma[1][3]])

df_gamma = pd.DataFrame(results_gamma, columns =['Likelihood', 'pneu', 'plet', 'alpha','beta'], dtype = float) #, 'alpha', 'beta'
df_gamma['Nanc'] = str(Nanc) 
df_gamma['nu'] = str(demo[0])
df_gamma['T'] = str(demo[1])
df_gamma['scale'] = 'na'
# Saving all the other parameters to a file
df_gamma.to_csv('inference/DFE_inference_gammalet_{let_prop}_{h}_{ss}_{replicate}.csv'.format(let_prop = let_prop, ss = sfs_size - 2, h=h, replicate=replicate), index=False)

neuexpolet_vec = numpy.frompyfunc(neuexpolet, 4, 1) #----5, 1
###pneu, plet, scale
lower_bound=[1e-3, 1e-3, 1e-3]         
upper_bound=[1,  1, 1] 
params = (.2, .2, .2) 
gamma_max_likelihoods = []
gamma_guesses = dict()

p0 = dadi.Misc.perturb_params(params, upper_bound=upper_bound)
popt = Selection.optimize_log(p0, nonsyn_sfs, spectra.integrate, neuexpolet_vec,
								 theta_ns, lower_bound=lower_bound, 
								 upper_bound=upper_bound, verbose=len(p0),
								 maxiter=25)

gamma_max_likelihoods.append(popt[0])
gamma_guesses[popt[0]] = popt

modelsfs = spectra.integrate(popt[1], neuexpolet_vec, theta_ns)  ###just added
#numpy.savetxt(model_file, modelsfs,
      #delimiter = ",")
numpy.savetxt('inference/computed_SFS_neuexpolet_{let_prop}_{h}_{ss}_{replicate}.txt'.format(let_prop = let_prop, ss=sfs_size - 2, h=h, replicate=replicate), modelsfs, delimiter=',') 

results_gamma = []
print(gamma_guesses.keys())
for i in range(len(gamma_guesses.keys())):
    best_popt_gamma = gamma_guesses[gamma_max_likelihoods[i]]
    results_gamma.append([best_popt_gamma[0], best_popt_gamma[1][0], best_popt_gamma[1][1], best_popt_gamma[1][2]])

df_gamma = pd.DataFrame(results_gamma, columns =['Likelihood', 'pneu','plet', 'scale'], dtype = float) #, 'alpha', 'beta'
df_gamma['Nanc'] = str(Nanc)
df_gamma['nu'] = str(demo[0])
df_gamma['T'] = str(demo[1])
df_gamma['alpha'] = 'na'
df_gamma['beta'] = 'na'
# Saving all the other parameters to a file
df_gamma.to_csv('inference/DFE_inference_neuexpolet_{let_prop}_{h}_{ss}_{replicate}.csv'.format(let_prop = let_prop, ss = sfs_size - 2, h=h, replicate=replicate), index=False)

