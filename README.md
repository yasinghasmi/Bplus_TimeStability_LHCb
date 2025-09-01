# Bplus_TimeStability_LHCb

Study of the time stability of the yield ratio between the
B+ → D̅0 π+ and B+ → J/ψ K+ decay modes in the 2024 LHCb dataset.

This work investigates the feasibility of using these decays for a precise determination of the branching fraction ratio

    B(B+ → D̅0 π+) / B(B+ → J/ψ K+)

focusing on two complementary aspects:

Temporal stability
- Yield ratios are studied across different blocks and fills of the 2024 run.
- The ratios are fitted to a constant value (0th-degree polynomial).
- The reduced chi-squared statistic (χ²/ndf) is used as the stability metric.

Branching fraction extraction
- The absolute value of the ratio is evaluated using efficiencies from Monte Carlo simulations.

Additional notes
- Mass fits are performed using the sum of two Crystal Ball functions.
- Separation of B+ → D̅0 π+ from B+ → J/ψ K+ is achieved with independent Crystal Ball fits.
- Stability is checked over four data-taking blocks and their corresponding fills.

Data structure

The repository is organized with an empty `data/` folder, which mirrors the layout used on EOS at
/eos/lhcb/user/y/yghasemi/B2JpsiKs/2024/.
Raw files placed in the correct directories can be used to regenerate all derived datasets.

data/
 ├─ processed/              # original input ROOT files (raw analysis data)
 ├─ processed_clean_bp_p/   # same files after duplicate-momentum cleaning
 ├─ fitted_data/            # mass-fit outputs and RooFit workspaces
 ├─ monte_carlo/            # simulated datasets for efficiency evaluation
 └─ real_5to8_raw/          # raw sub-samples of 2024 real data (fills 5–8)

Each subfolder is automatically created by the scripts when needed, and the contents are regenerated from upstream data.

License: MIT
