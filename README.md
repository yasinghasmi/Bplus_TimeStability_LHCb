# Bplus_TimeStability_LHCb

Study of the time stability of the yield ratio between the  
$B^+ \to \overline{D}^0 \pi^+$ and $B^+ \to J/\psi K^+$ decay modes in the **2024 LHCb dataset**.

---

## Goal

This project evaluates the feasibility of using these decays for a precise measurement of the branching fraction ratio:

\[
\frac{\mathcal{B}(B^+ \to \overline{D}^0 \pi^+)}{\mathcal{B}(B^+ \to J/\psi K^+)}
\]

The analysis focuses on two main aspects:

### 1. Temporal stability
- Studying the yield ratio across different **blocks** and **fills** of the 2024 run.  
- Fitting the ratio to a constant value (0th-degree polynomial).  
- Using the reduced chi-squared statistic,  
  \[
  \chi^2/\text{ndf},
  \]  
  as the stability metric.

### 2. Branching fraction extraction
- Evaluating the absolute value of the ratio using **efficiencies determined from Monte Carlo simulations**.

---

## Notes

- **Mass fits** are performed using the sum of two Crystal Ball functions.  
- For separating $B^+ \to \overline{D}^0 \pi^+$ from $B^+ \to J/\psi K^+$, independent Crystal Ball fits are applied.  
- Stability checks are performed across **four data-taking blocks** and their corresponding fills.  

---

## License

This project is licensed under the **MIT License**.
