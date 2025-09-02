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
 ├─ real_5to8_raw/          # 0. Raw files based on the given naming for blocks 5 to 8
 ├─ processed/              # 1. Files separated based on fill and block, first step 
 |                            of cleaning.
 ├─ processed_clean_bp_p/   # 2. same files after duplicate-momentum cleaning, after running 
 |                            the fitting code it will contain the fit result as well.
 ├─ monte_carlo/            # MC.1. Original MC data
 ├─ monte_carlo_processed/  # MC.2. Merged MC data, to form blocks 
 └─ block5_analysis/        # 3. This file is used for block 5 investigation analysis
    └─ outputs                that is epxlained more deeply in the md of file. The data
       ├─ fit_plots           beofre fitting should be places here (next to outputs), and result
       └─ histograms          of analysis will be saved in the output file, in two subdirectories.
```

Each file comes with its own accompanying md block explaining the details, but here’s a clearer overview of their purpose:

- **split_fills_clean.ipynb**: The first step of data preparation. It takes the raw data, splits it into block- and fill-level files, and removes bad runs.

- **clean_bp_p_duplicates.ipynb**: Uses the split block/fill files to clean up duplicate candidates. It clusters entries with the same event number based on the Bp_P branch and keeps only one representative per cluster.

- **blocks_fitting.ipynb**: Performs block-level mass fits for B⁺ decays. It applies simultaneous fits to B⁺ → D⁰π⁺ and the combined B⁺ → J/ψK⁺ + J/ψπ⁺ channels using double Crystal Ball functions with exponential backgrounds. It processes multiple ROOT files from different data-taking blocks and extracts yields, fit parameters, and fit quality metrics. This code UPDATEs the existing files in the same directory.

- **fills_fitting.ipynb**: Similar to blocks_fitting.ipynb, but at the fill level. Since a fill is a shorter data-taking period within a block, this analysis provides finer time resolution for stability studies.

- **blocks_yield_stability.ipynb**: Studies how B⁺ yields evolve over different data-taking blocks. It looks at trends in signal yields, background fractions, and fit quality over time.

- **MC_file_merger.ipynb**: Prepares Monte Carlo samples by merging multiple ROOT files into unified datasets. It applies event selection and cuts, and organizes the data by decay channel (B2OC, B2CC, J/ψπ). This creates clean MC samples for later parameter extraction and systematic studies.

- **MC_fitting.ipynb**: Fits the Monte Carlo datasets to extract the signal shape parameters (Crystal Ball means, widths, and tail parameters). Since MC samples are free of background, these parameters can later be used to constrain fits in real data, improving robustness and reducing systematics.

- **branching_ratio_block_plotting.ipynb**: Plots the branching ratio and total error (statistical + external) for all 4 blocks.

- **efficiency_ratio_block.ipynb**: Plots efficiency ratios for each block, using values manually provided from branching_ratio_calculation.ipynb.

- **nPV_per_fill.ipynb**: Checks the number of primary vertices in each fill. It produces plots of the average and standard deviation of nPV across blocks.

- **branching_ratio_calculation.ipynb**: Calculates the branching ratio and its uncertainties (both statistical and external). It also computes yield and efficiency ratios based on the fitted data.

- **branching_ratio_vs_world_plot.ipynb**: Compares the branching fraction results against world averages (PDG) and measurements from other experiments such as Belle and BaBar.

- **block5_analysis.ipynb**: Loads Block 5 and per fill ROOT data, performs B⁺ mass fits (2CB+exp) to write sWeights, computes χ² diagnostics, plots per fill normalized PID histograms, and generates a by fill signal/background table and graph.

**License**: MIT
