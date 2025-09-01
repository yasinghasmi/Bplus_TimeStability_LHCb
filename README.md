# Bplus_TimeStability_LHCb

Study of the time stability of the yield ratio between the  
$B^+ \to \overline{D}^0 \pi^+$ and $B^+ \to J/\psi K^+$ decay modes in the **2024 LHCb dataset**.

This work investigates the feasibility of using these decays for a precise determination of the branching fraction ratio

$$
\frac{\mathcal{B}(B^+ \to \overline{D}^0 \pi^+)}{\mathcal{B}(B^+ \to J/\psi K^+)}
$$

focusing on two complementary aspects:

### Temporal stability
- Yield ratios are studied across different **blocks** and **fills** of the 2024 run.  
- The ratios are fitted to a constant value *(0th-degree polynomial)*.  
- The reduced chi-squared statistic is used as the stability metric.

### Branching fraction extraction
- The absolute value of the ratio is evaluated using **efficiencies from Monte Carlo simulations**.

### Additional notes
- Mass fits are performed using the sum of two Crystal Ball functions.  
- Separation of $B^+ \to \overline{D}^0 \pi^+$ from $B^+ \to J/\psi K^+$ is achieved with independent Crystal Ball fits.  
- Stability is checked over four data-taking blocks and their corresponding fills.  

### Data structure

The repository is organized with an empty `data/` folder, which mirrors the layout exisitng on EOS at  
`/eos/lhcb/user/y/yghasemi/B2JpsiKs/2024/`.  
Raw files placed in the correct directories can be used to regenerate all derived datasets.

```
data/
 ├─ processed/              # 1. Files separated based on fill and block, first step of cleaning
 ├─ processed_clean_bp_p/   # 2. same files after duplicate-momentum cleaning
 ├─ fitted_data/            # 3. mass-fit outputs and RooFit workspaces
 ├─ monte_carlo/            # MC.1. Original MC data
 ├─ monte_carlo_processed/  # MC.1. Merged MC data, to form blocks 
 └─ real_5to8_raw/          # 0. Raw files based on the given naming for blocks 5 to 8
```
Each file contains a corresponding ```md``` block with the explanation.

Codes descriptions for External Users:
- Blocks_Fitting.ipynb: Performs block-level mass fitting analysis on B+ meson decay data. Implements simultaneous fits for B+ → D⁰π+ and combined B+ → J/ψK+ + J/ψπ+ channels using double Crystal Ball PDFs with exponential backgrounds. Processes multiple ROOT files corresponding to different data-taking periods (blocks), extracting yields, fit parameters, and quality matrices.

- Blocks_Yield_Stability.ipynb: Analyzes temporal stability of B+ yields across different data-taking blocks. Examines systematic variations in signal yields, background fractions, and fit quality metrics as a function of time/block number.

- Branching_Fraction_Calculation.ipynb: Computes branching fraction ratios BF(B+ → D⁰π+)/BF(B+ → J/ψK+) using fitted yields from the mass fitting analysis. Implements proper error propagation, efficiency corrections, and systematic uncertainty evaluation. Compares results with world averages (PDG) and other experiments (Belle, BaBar).

- Fills_Fitting.ipynb: Performs fill-level mass fitting analysis - a more granular analysis than block-level fitting. Each fill represents a shorter data-taking period within a block, providing higher time resolution for stability studies. Uses similar fitting procedures as Blocks_Fitting.ipynb but processes individual fills.

- MC_File_Merger.ipynb: Handles Monte Carlo data preprocessing by merging multiple MC ROOT files into consolidated datasets. Ensures proper event selection, applies necessary cuts, and organizes MC samples by decay mode (B2OC, B2CC, J/ψπ). Critical for preparing clean MC samples used in parameter extraction and systematic studies.

- MC_Fitting.ipynb: Performs Monte Carlo mass fitting to extract signal shape parameters (Crystal Ball parameters, resolutions, tail parameters) from pure MC samples without background contamination. These parameters are then used to constrain real data fits, improving stability and reducing systematic uncertainties.









**License**: MIT
