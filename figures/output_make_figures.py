#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:54:30 2023
@author: Marine Laporte - Université Claude Bernard Lyon 1
@last-update: 24/10/2024
@contact: marinelaporte.ensg@gmail.com
"""

## Libraries we use 
## ------------------
import sys 
import os
import time

import plotly.express as px
import numpy as np
from math import exp,sqrt, erf
from mpmath import mp
from datetime import date 
import scipy.signal as ss
import scipy.stats as st
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

from matplotlib import cm
from matplotlib.pyplot import cm
from cmcrameri import cm
## ------------------
from subfunctions.imports_outfiles import FileReader 



def compute_weighted_mean_std(dist, likelihood):
    """
    Calcule la moyenne pondérée et l'écart-type pour chaque colonne d'une matrice de vraisemblance.

    Parameters:
    - dist: np.ndarray, les valeurs du paramètre (axe vertical, ex: b, mu, sigma)
    - likelihood: np.ndarray, la matrice de vraisemblance de forme (len(dist), N)

    Returns:
    - mean_vals: list of float, moyenne pondérée pour chaque colonne
    - std_vals: list of float, écart-type pondéré pour chaque colonne
    """
    mean_vals = []
    std_vals = []

    for i in range(likelihood.shape[1]):
        weights = likelihood[:, i]
        total_weight = np.sum(weights)
        if total_weight == 0:
            mean_vals.append(np.nan)
            std_vals.append(np.nan)
            continue

        mean = np.sum(weights * dist) / total_weight
        std = sqrt(np.sum(weights * (dist - mean) ** 2) / total_weight)

        mean_vals.append(mean)
        std_vals.append(std)

    return mean_vals, std_vals

def B_Positive_function(times, mags, w):
    """Fonction qui calcule la valeur de la b-value
    a partir de la formule B positive du papier de Van Der Elst (2021)
    utilisant les différences positives de magnitudes suivant la distribution de Laplace
    Inputs :
        - mags : liste ou array de magnitudes (floats)
        - times : liste ou array des temps (floats) associés aux magnitudes 
        - w : taille de la fenetre (nb d'evenements/fenetre)
    Outputs :
        - b_positive : liste ou array de valeurs de b-positives
        - times_b_positive : temps associés aux valeurs de b_positive """
    # INITIALISATION
    n_events =len(mags)
    differences_magnitudes = np.zeros(n_events-1)
    differences_absolues = np.zeros(n_events-1)
    # init outputs
    b_positive= np.zeros(n_events-w)
    times_b_positive= np.zeros(n_events-w)
    
    # Calcul des differences de magnitudes sur tout le catalogue
    for i in range (0, n_events-1):
        differences_magnitudes[i] = mags[i+1] - mags[i]       
        differences_absolues[i] = abs(mags[i+1] - mags[i])
        
    for m in range(len(differences_magnitudes)-(w-1)):
       differences_positives=[]
       flag_pos =0
       flag =0
       Minimum_difference = min(differences_absolues[m:m+w]+0.4)
       for valeur in differences_magnitudes[m:m+w]:
           flag+=1
           if valeur > Minimum_difference:
               flag_pos+=1 # Pour savoir sur combien d evenements on calcule effectivement la b-value
               differences_positives.append(valeur)
       windowed_mean = sum(np.ones(len(differences_positives)) * differences_positives) / len(differences_positives)
       
       
       # d'apres la formule de Laplace (Van der Elst 2021)
       windowed_beta = 1/float(windowed_mean-Minimum_difference) 
       b_positive[m]=windowed_beta/np.log(10)  # On retrouve b_value = beta/2.3

    print("Valeur de b-value obtenue sur %s events sur fenetre de %s events" %(str(flag_pos),str(w)))
    times_b_positive = times[int(w/2):-int(w/2)]
    return times_b_positive, b_positive

def load_all_data(year0, Prop_out, Nc_out, Dataset, Mc_out, b_out, S_out, Acc_out, Espace_out):
    reader = FileDataReader()
    
    prop_vec = reader.read_single_column_file(Prop_out)
    nc_vecs = reader.read_matrix_file(Nc_out)
    times, mags = reader.read_dataset(Dataset)
    times2 = [t + year0 for t in times]
    Likelihood_mu = reader.read_matrix_file(Mc_out)
    Likelihood_b = reader.read_matrix_file(b_out)
    Likelihood_sigma = reader.read_matrix_file(S_out)
    n_acc, n_tot, pct = reader.read_accuracy_file(Acc_out)
    Mu_prior, b_prior, Sigma_prior, it_hyp = reader.read_hypercube_file(Espace_out)

    return {
        "prop_vec": prop_vec,
        "nc_vecs": nc_vecs,
        "times": times,
        "times2": times2,
        "mags": mags,
        "Likelihood_mu": Likelihood_mu,
        "Likelihood_b": Likelihood_b,
        "Likelihood_sigma": Likelihood_sigma,
        "n_acc": n_acc,
        "n_tot": n_tot,
        "pct": pct,
        "Mu_prior": Mu_prior,
        "b_prior": b_prior,
        "Sigma_prior": Sigma_prior,
        "it_hyp": it_hyp
    }

if __name__ == "__main__" :
    
    
    # Get today's date and time
    today = time.strftime('%d_%m_%Y_%Hh')
    
    # Define constants
    name = 'Kumamoto_2025'
    year0 = 0
    synthetics = 'False'
    Num_run = 1
    

# =============================================================================
#   Default outfiles names in OUT/
# =============================================================================
    Prop_out ='../outputs/Probability_discontinuities.txt'
    Nc_out ='../outputs/Number_of_discontinuities.txt'
    Dataset ='../outputs/input.txt'
    Mu_out ='../outputs/Marginal_Mu.txt'
    b_out ='../outputs/Marginal_b.txt'
    Sigma_out ='../outputs/Marginal_Sigma.txt'
    Acc_out ='../outputs/Acceptances_Rates.txt'
    Espace_out ='../outputs/prior.txt'
# =============================================================================
#   If run on synthetics - give synthetic parameters to retrieve
# =============================================================================


    # T_vec = [0, 25, 55, 100 ]
    # N_vec = [200,600,900]
    # b_vec = np.array([0.8,1.2,0.9])*2.3
    # mu_vec = np.array([0.5, 0.5, 2])
    # sigma_vec = np.array([0.3, 0.2, 0.25])
    # Ntot = sum(N_vec)
# =============================================================================
#     Reading outfiles  
# =============================================================================
    reader = FileReader()
    prop_vec = reader.read_single_column_file(Prop_out)
    Nc_vecs = reader.read_matrix_file(Nc_out)
    times, mags = reader.read_dataset(Dataset)
    #times2 = [t + year0 for t in times]
    Likelihood_mu = reader.read_matrix_file(Mu_out)
    Likelihood_b = reader.read_matrix_file(b_out)
    Likelihood_sigma = reader.read_matrix_file(Sigma_out)
    N_acc, N_tot, pct = reader.read_acceptance_file(Acc_out)
    Mu_prior, b_prior, Sigma_prior, it_params = reader.read_hypercube_file(Espace_out)
    
    Nc_vec = Nc_vecs[:,0]
    N_iter = len(Nc_vec)
    N_tirs = int(it_params[1])
    N_mags = len(mags)
    tmax = max(times)
    if N_mags ==1 :
        tmax= max(times)*2

  
    # Define common time range
    time_range = np.linspace(year0, year0 + tmax, len(prop_vec))
    bins = time_range
    bins_20 = np.linspace(year0, year0 + tmax, 20)

      

# =============================================================================
#     # Create output path
# =============================================================================    
    
    out_path = f'Run_{today}'
    os.makedirs(out_path, exist_ok=True)
    
    # Get today's date again for formatted output
    today_date = date.today()
    Date = today_date.strftime("%d%m")
    
    # Define titles
    Title_1 = f'Transdimensional RJ-MCMC - Synthetic dataset {N_mags} evts'
    Subtitle = f'{today_date.strftime("%d/%m")} : Run {Num_run} - {N_iter} models'
    
    # Output filename construction
    out_name = f'{name}_{N_iter}it_{N_tirs}rd_'
    

# =============================================================================
# Retrieve the location of the discontinuities
# =============================================================================
    Nc_counts, Nc_edges = np.histogram(Nc_vec, bins=np.arange(0, 20, 1))
    Nc_max = Nc_edges[np.argmax(Nc_counts)]
    nbins_cat = max(times) / len(prop_vec)
    
    disc_calc = []
    peaks, properties = ss.find_peaks(prop_vec, height=1, prominence=1)
    prominences = properties['prominences']
    
    for _ in range(int(Nc_max) - 1):
        max_idx = np.argmax(prominences)
        disc_calc.append(peaks[max_idx] * nbins_cat)
        prominences[max_idx] = 0  # Suppress the current max to find the next one
        
        
# =============================================================================
# Compute weighted mean and standard deviation of marginal distributions
# =============================================================================        
    # Pour mu
    Mu_dist = np.linspace(Mu_prior[0], Mu_prior[1], int(Mu_prior[2]))
    Mean_mu, Std_mu = compute_weighted_mean_std(Mu_dist, Likelihood_mu)
    
    # Pour b
    B_dist = np.linspace(b_prior[0]/np.log(10), b_prior[1]/np.log(10), int(b_prior[2]))
    Mean_b, Std_b = compute_weighted_mean_std(B_dist, Likelihood_b)
    
    # Pour sigma
    Sigma_dist = np.linspace(Sigma_prior[0], Sigma_prior[1], int(Sigma_prior[2]))
    Mean_S, Std_S = compute_weighted_mean_std(Sigma_dist, Likelihood_sigma)
     
# =============================================================================
#  Decompose the input catalogue at each discontinuities
# =============================================================================  

    fig0, ax0 = plt.subplots(figsize=(10, 8))
    ax0b = ax0.twinx()
    # Scatter plot: magnitudes
    ax0b.scatter(times, mags,
                 c=times,
                 s=[exp(m) * 2 for m in mags],
                 edgecolor='black',
                 linewidth=0.25,
                 cmap=cm.vik,
                 alpha=0.8)
    ax0b.set_ylabel('Magnitudes', fontsize=17, color='dimgrey')
    ax0b.set_ylim(0, 8)
    ax0b.tick_params(axis='y', colors='dimgrey', labelsize=17)
    
    # Normalize and plot probability vector
    Prop_norm = [y / max(prop_vec) for y in prop_vec]
    ax0.plot(bins, Prop_norm, lw=3, color='black', label='Inversion McMC', zorder=2)
    # Plot synthetic discontinuities if applicable
    if synthetics == 'True':
        for i, t in enumerate(T_vec):
            ax0.vlines(t + year0, 0, 1,
                       linestyle=':', lw=2, color='green',
                       label='Fixed discontinuities' if i == 0 else None)
    
    # Formatting axes
    ax0.set_xlim(min(times), max(times))
    ax0.set_ylim(0, 1)
    ax0.set_xlabel('Time', fontsize=17)
    ax0.set_ylabel('Probability of discontinuities', fontsize=17, color='black')
    ax0.tick_params(axis='x', colors='black', labelsize=17, rotation=20)
    ax0.tick_params(axis='y', colors='black', labelsize=17)
    ax0.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    
    # Legend and save
    ax0.legend(loc='upper right', fontsize=17)
    fig0.savefig(os.path.join(out_path, out_name + 'Fig0.pdf'))
# =============================================================================
#     FIGURE 1 : Informative plot on acceptance rates
# =============================================================================
    fig1, axs = plt.subplots(2, 2, figsize=(10, 8))
    ax1, ax2, ax3, ax4 = axs[0, 0], axs[1, 0], axs[0, 1], axs[1, 1]
    plt.subplots_adjust(wspace=0.4, hspace=0.2)
    
    # --- Panel 1: Probability of discontinuities over time ---
    ax1.plot(bins, Prop_norm, color='black', zorder=2)
    ax1.set_xlabel('Dataset times')
    ax1.set_ylabel('Nb of acceptance for that time')
    ax1.grid(lw=0.3)
    
    if synthetics == 'True':
        ax1b = ax1.twinx()
        for t in T_vec:
            ax1b.vlines(t + year0, 0, 1, lw=1.5, linestyle='--', color='green')
    
    # --- Panel 2: Evolution of Nc over iterations ---
    bins_Nc = np.linspace(0, len(Nc_vec), len(Nc_vec))
    color_iter = iter(cm.grayC(np.linspace(0, 1, Nc_vecs.shape[1])))
    for ni in range(Nc_vecs.shape[1]):
        ax2.plot(bins_Nc, Nc_vecs[:, ni], color=next(color_iter), lw=0.6, alpha=0.5)
    
    mean_Nc = np.mean(Nc_vecs, axis=1)
    # Uncomment below if you want to include mean/std shading
    # std_Nc = np.std(Nc_vecs, axis=1)
    # ax2.plot(bins2, mean_Nc, color='black', lw=0.8, alpha=0.9)
    # ax2.fill_between(bins2, mean_Nc - std_Nc, mean_Nc + std_Nc, color='dimgrey', alpha=0.6)
    
    ax2.set_xlabel('Iteration number')
    ax2.set_ylabel('Nb of accepted discontinuities Nc')
    ax2.grid(lw=0.6)
    
    # --- Panel 3: Bar plot of total vs accepted models ---
    labels = ['Births', 'Deaths', 'Perts']
    x_pos = np.arange(len(N_tot))
    
    ax3.bar(labels, N_tot, color='darkslateblue', edgecolor='black', alpha=0.5, label=f'Nb total : {sum(N_tot)} models')
    ax3.bar(labels, N_acc, color='greenyellow', edgecolor='black', lw=2, alpha=0.6, label=f'Nb accepted : {sum(N_acc)} models')
    
    for i in range(3):
        ax3.text(i - 0.25, N_tot[i] / 3, f'Na={N_acc[i]}')
        ax3.text(i - 0.25, N_tot[i] / 5, f'Ta={round(pct[i] * 100, 1)}%')
    
    ax3.set_ylabel('Nb of tot/acc models')
    ax3.legend(loc='upper center')
    
    # --- Panel 4: Histogram of Nc values ---
    bins_Nc = np.arange(0, max(Nc_vec), 1)
    ax4.hist(
        Nc_vecs.flatten(),
        bins=bins_Nc,
        color='firebrick',
        alpha=0.6,
        edgecolor='darkslateblue',
        lw=2,
        orientation='horizontal',
        align='left',
        label='Nc'
    )
    ax4.set_ylim(0, max(Nc_vec) + 2)
    ax4.set_ylabel('Nb of accepted discontinuities Nc')
    ax4.set_xlabel('Number of models')
    ax4.legend(loc='upper right')
    
    # --- Final adjustments ---
    fig1.suptitle(f"{Title_1}\n{Subtitle}", fontweight="bold")
    fig1.savefig(os.path.join(out_path, out_name + 'Fig1.pdf'))
    plt.show()
    
# =============================================================================
#     FIGURE 2 - Likelihood Mc 
# =============================================================================
    fig2, ax = plt.subplots(figsize=(8, 8))
    
    # --- Main heatmap: Likelihood of μ ---
    im = ax.imshow(
        Likelihood_mu,
        aspect='auto',
        cmap='Reds',
        extent=(year0, tmax + year0, Mu_prior[1], Mu_prior[0])
    )
    ax.set_ylim(Mu_prior[0], Mu_prior[1])
    ax.set_ylabel(r'Values of $\mu$', color='firebrick', fontsize=18)
    ax.set_xlabel('Time', fontsize=18)
    ax.tick_params(axis='y', colors='firebrick', labelsize=18)
    ax.tick_params(axis='x', colors='black', labelsize=18)
    ax.spines['left'].set_color('red')
    ax.grid(lw=0.3)
    
    # --- Overlay: Mean and std of μ ---
    ax.plot(bins, Mean_mu, color='red', linestyle='-.', lw=1, alpha=0.7, label='Mean μ')
    ax.plot(bins, [m + s for m, s in zip(Mean_mu, Std_mu)], color='black', linestyle='--', lw=1.5, alpha=0.7)
    ax.plot(bins, [m - s for m, s in zip(Mean_mu, Std_mu)], color='black', linestyle='--', lw=1.5, alpha=0.7)
    
    # --- Optional: Synthetic discontinuities ---
    if synthetics == 'True':
        for i in range(len(T_vec) - 1):
            ax.hlines(mu_vec[i], T_vec[i], T_vec[i + 1], color='red', linestyle='--', lw=3)
            ax.vlines(T_vec[i] + year0, Mu_prior[0], Mu_prior[1], color='lightgreen', linestyle=':', lw=3, label='Fixed discontinuities')
    
    # --- Secondary y-axis: Magnitude scale ---
    ax2 = ax.twinx()
    ax2.set_ylim(0, 7)
    ax2.set_ylabel('Magnitudes MLv', color='dimgrey', fontsize=18)
    # ax2.tick_params(axis='y', colors='dimgrey', labelsize=14)  # Uncomment if styling needed
    
    # --- Title and save ---
    Title2 = fr'Likelihood $\mu$ - Synthetic dataset {N_mags} evts'
    fig2.suptitle(f"{Title2}\n{Subtitle}", fontweight="bold")
    fig2.savefig(os.path.join(out_path, out_name + '_Fig2_likelihoodMu.pdf'))
    plt.show()
    
# # =============================================================================
# #     FIGURE 3 - Likelihood b
# # =============================================================================   
    # Compute true b-values
    b_true = Likelihood_b / np.log(10)
    
    # Create figure and axis
    fig3, axb = plt.subplots(figsize=(8, 11))
    
    # Heatmap of b-values
    imb = axb.imshow(
        b_true,
        aspect='auto',
        cmap='Reds',
        extent=(year0, tmax + year0, b_prior[1]/np.log(10), b_prior[0]/np.log(10))
    )
    
    # Plot mean and ±1σ std deviation
    axb.plot(bins, Mean_b, c='red', linestyle='-.', lw=1, alpha=0.9, label='Mean b-value')
    axb.plot(bins, [m + s for m, s in zip(Mean_b, Std_b)], c='black', linestyle='--', lw=1, alpha=0.7, label='±1σ')
    axb.plot(bins, [m - s for m, s in zip(Mean_b, Std_b)], c='black', linestyle='--', lw=1, alpha=0.7)
    
    # Optional synthetic data annotations
    if synthetics == 'True':
        for i in range(len(T_vec) - 1):
            label = 'Imposed values' if i == 0 else None
            axb.hlines(b_vec[i] / np.log(10), T_vec[i], T_vec[i + 1], color='red', linestyle='--', lw=3, label=label)
            axb.vlines(T_vec[i] + year0, b_prior[0]/np.log(10), b_prior[1]/np.log(10), color='lightgreen', linestyle=':', lw=3)
    
    # Secondary axis for magnitudes
    axb2 = axb.twinx()
    axb2.set_ylim(0, 7)
    axb2.set_ylabel('Magnitudes MLv', color='dimgrey', fontsize=18)
    
    # Axis styling
    axb.set_ylim(b_prior[0]/np.log(10), b_prior[1]/np.log(10))
    axb.set_ylabel('b values', color='firebrick', fontsize=14)
    axb.set_xlabel('Time', fontsize=18)
    axb.grid(lw=0.3)
    axb.spines['left'].set_color('firebrick')
    axb.tick_params(axis='y', colors='firebrick', labelsize=18)
    axb.tick_params(axis='x', colors='black', labelsize=18)
    
    # Add legend and colorbar
    axb.legend(loc='upper right', fontsize=12)
    fig3.colorbar(imb, ax=axb, orientation='horizontal')
    
    # Add title and save
    Title3 = f'Likelihood b - Synthetic dataset {N_mags} evts'
    fig3.suptitle(f"{Title3}\n{Subtitle}", fontweight="bold")
    fig3.savefig(os.path.join(out_path, out_name + '_Fig3_likelihoodb.pdf'))
    plt.show()
    
# # =============================================================================
# #     FIGURE 4 - Likelihood S
# # =============================================================================   

    fig4, ax = plt.subplots(figsize=(8, 8))
    
    # Heatmap of likelihood for Sigma
    im = ax.imshow(
        Likelihood_sigma,
        aspect='auto',
        cmap='Reds',
        extent=(year0, tmax + year0, Sigma_prior[1], Sigma_prior[0])
    )
    
    # Plot mean and ±1σ
    ax.plot(bins, Mean_S, c='red', linestyle='-.', lw=1, label='Mean Sigma')
    ax.plot(bins, [m + s for m, s in zip(Mean_S, Std_S)], c='black', linestyle='--', lw=1, alpha=0.7, label='+1σ')
    ax.plot(bins, [m - s for m, s in zip(Mean_S, Std_S)], c='black', linestyle='--', lw=1, alpha=0.7, label='-1σ')
    
    # Secondary y-axis for magnitudes
    ax2 = ax.twinx()
    ax2.set_ylim(0, 7)
    ax2.set_ylabel('Magnitudes MLv', color='dimgrey', fontsize=18)
    
    # Set plot limits and grid
    ax.set_ylim(Sigma_prior[0], Sigma_prior[1])
    ax.set_xlabel('Time', fontsize=18)
    ax.set_ylabel('Values of Sigma', color='firebrick', fontsize=18)
    ax.grid(lw=0.3)
    
    # Aesthetic adjustments
    ax.spines['left'].set_color('firebrick')
    ax.tick_params(axis='y', colors='firebrick', labelsize=18)
    ax.tick_params(axis='x', colors='black', labelsize=18)
    
    # Plot synthetic references if applicable
    if synthetics == 'True':
        for i in range(len(T_vec) - 1):
            ax.hlines(sigma_vec[i], T_vec[i], T_vec[i + 1], color='red', linestyle='--', lw=3)
            ax.vlines(T_vec[i] + year0, Sigma_prior[0], Sigma_prior[1], linestyle=':', lw=3, color='lightgreen', label='Fixed discontinuities')
    
    # Title and save
    Title4 = f'Likelihood Sigma - Synthetic dataset {N_mags} evts'
    fig4.suptitle(f"{Title4}\n{Subtitle}", fontweight="bold")
    fig4.savefig(os.path.join(out_path, out_name + '_Fig4_likelihood_Sigma.pdf'))
    plt.show()