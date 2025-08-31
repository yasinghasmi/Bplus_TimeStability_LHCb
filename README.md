# Bplus_TimeStability_LHCb

Study of the time stability of the yield ratio between the  
**$B^+ \to \overline{D}^0 \pi^+$** and **$B^+ \to J/\psi K^+$** decay modes in the 2024 LHCb dataset.

---

## Goal

This project evaluates the feasibility of using these decays for a precise measurement of the branching fraction ratio:

\[
\frac{\mathcal{B}(B^+ \to \overline{D}^0 \pi^+)}{\mathcal{B}(B^+ \to J/\psi K^+)}
\]

The analysis focuses on two aspects:

1. **Temporal stability**  
   – Studying the yield ratio across different blocks and fills of the 2024 run.  
   – Fitting the ratio to a constant value (0th-degree polynomial).  
   – Using the reduced chi-squared, $\chi^2/\text{ndf}$, as the stability metric.

2. **Branching fraction extraction**  
   – Evaluating the absolute value of the ratio using efficiencies determined from Monte Carlo simulations.

---

## Method

- **Mass fits** are performed using the sum of two Crystal Ball functions.  
- For separation of $B^+ \to D^0 \pi^+$ from $B^+ \to J/\psi K^+$, independent Crystal Ball fits are applied.  
- Stability is checked over **four data-taking blocks** in 2024 and their corresponding fills.  

---

## Data Directory

The `data/` directory in this repository is empty by default.  
It mirrors the structure where datasets should be placed locally.  

The same structure exists at:  
`/eos/lhcb/user/y/yghasemi/B2JpsiKs/2024/`  

Data can be fetched per step or per file from EOS.  
Using the raw data in the correct subdirectories, the other derived datasets can be produced.

---

## Data Structure (expected)

Each folder under `data/` represents a different processing stage.  
Explanations of their purpose and expected contents will be added here.  
Currently the hierarchy is empty but serves as a template.

---

## File Documentation

Each code file in this project starts with a **Markdown block** at the top that explains:
- the details of the analysis step it performs, and
- the main functions implemented.

This ensures reproducibility and clarity of workflow.

---

## License

This project is licensed under the MIT License.
