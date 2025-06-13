# b_Bayesian_binned - v2 - 2025 -

## Probabilistic Estimate of b-value Temporal Variations  
**Last update:** 11/06/2025 
**Author:** Marine Laporte  
**Contact:** [marine.laporte@univ-lyon1.fr](mailto:marine.laporte@univ-lyon1.fr)  
**Version:** Version 2 - added magnitude binning and quicker computation time 
---

## 📚 Citation  
If you use this code, please cite:  

> Laporte, M., Durand, S., Bodin, T., Gardonio, B., & Marsan, D. (2025). *b-bayesian: The full Probabilistic estimate of b-value Temporal variations for non-truncated catalogs.* Journal of Geophysical Research: Solid Earth, 130, e2024JB029973. https://doi.org/10.1029/2024JB029973 
---

## 🧭 Description

**b-Bayesian** performs a **transdimensional Bayesian inversion** of the temporal evolution of:

- **b-value** from the Gutenberg-Richter law (beta= b-value/ln(10))
- **μ (mu)** and **σ (sigma)**, which describe detectability using the error function (Ogata and Katsura, 1993)

<<<<<<< HEAD

=======
>>>>>>> 7757129f85fd96c67424ce5ea41d625d7b95af40
>Difference of this new version *binned*:
*b-Bayesian_binned* is a similar version of b-Bayesian that considers the binning of an earthquake catalog (usually 0.1), this parameterization reduces considerably the computational time for earthquake catalogs of several thousands of events (the likelihood is not the product over the number of events but over the number of bins). 
---

### 🔁 Outputs

The code returns:

- `Probability_discontinuities.txt`: Probability of a temporal discontinuity for each time bin
- `Number_discontinuities.txt`: Number of discontinuities at each MCMC iteration
- `Marginal_b.txt`: Full probability distribution of b-value for each time bin
- `Marginal_Mu.txt`: Full probability distribution of mu for each time bin
- `Marginal_Sigma.txt`: Full probability distribution of sigma for each time bin

---

## ⚙️ Installation

### Dependencies

- **To run b-Bayesian**: a Fortran90 compiler (choose it in the Makefile)
- **To generate figures**: Python 3 with the following libraries:
  - `numpy`
  - `matplotlib`
  - `crameri`
  - `math`

---

## 🔄 Running the Inversion

### 1. Set Input File and Parameters

Edit `/inversion/b_Bayesian.f90`:

- Set your dataset path:  
  `Magnitude_file = "your/data/path"`
- Adjust MCMC parameters (iterations, burn-in, etc.)
- Set prior distributions for **β**, **μ**, and **σ** (start with broad priors)

### 2. Compile and Run

In the `/inversion/` folder:

```bash
make
./run
```

---

## 📂 Output Files

b-Bayesian generates 8 text files in the `/outputs/` directory:

| File             		    | Description |
|-----------------------------------|-------------|
| `input.txt`      		    | Copy of the input earthquake catalog |
| `prior.txt`	    		    | Parameters used in the run (priors, iterations) |
| `Probability_discontinuities.txt` | Probability of discontinuity per time bin (size: T_bins) |
| `Number_discontinuities.txt`      | Number of discontinuities per iteration (size: itmax) |
| `Marginal_b.txt`		    | Probability density of b-value (T_bins × B_bins) |
| `Marginal_mu.txt`		    | Probability density of mu (T_bins × Mu_bins) |
| `Marginal_sigma.txt`		    | Probability density of sigma (T_bins × Sigma_bins) |
| `Acceptances_Rates.txt` 	    | 3×3 matrix: birth/death/move acceptance statistics:  
  Columns = [Birth, Death, Move]  
  Rows = [Accepted, Proposed, Acceptance rate] |

---

## 📈 Plotting Figures

The script `output_figures.py` (in `/figures`) reads the inversion outputs and generates visualizations. 

### Steps:

1. Rename output files accordingly in the main section of `output_figures_bMuSigma.py`
2. Run the script:

```bash
python3 output_figures_bMuSigma.py
```

### Output:

Creates a new folder in `/figures/` named `/Run_day_month_year_hour` containing:

- **Fig0**: Magnitudes over time with overlaid probability of discontinuities  
- **Fig1**: 4-panel summary:
  - (Top-left): Final probability of discontinuities over time bins
  - (Top-right): Accepted vs. proposed models (birth/death/move)
  - (Bottom-left): Number of discontinuities per iteration
  - (Bottom-right): Histogram of the number of discontinuities (peak = most probable number)
- **Fig2**: Probability density of **b-value** over time
- **Fig3**: Probability density of **mu** over time
- **Fig4**: Probability density of **sigma** over time

---

## 📝 License

Academic use only. Please contact the author for permission to use, modify, or distribute this code.

