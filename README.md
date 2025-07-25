# Bplus_TimeStability_LHCb

Study of the time stability of the yield ratio between the  
**B⁺ → D̄⁰π⁺** and **B⁺ → J/ψK⁺** decay modes in the 2024 LHCb dataset.

This project evaluates the feasibility of using these decays for a precise measurement of the branching fraction ratio:

\[
\frac{\mathcal{B}(B^+ \to \overline{D}^0\pi^+)}{\mathcal{B}(B^+ \to J/\psi K^+)}
\]

by analyzing the temporal stability of their yield ratio across multiple data-taking blocks.

Mass fits are performed using the sum of two Crystal Ball functions, and time stability is evaluated by binning the dataset into quantile-based time intervals and fitting the resulting yield ratios to a constant value using a 0th-degree polynomial:

\[
\chi^2/\text{ndf} \quad \text{is used as a stability metric for each block.}
\]

The aim is to establish a stable baseline that enables future high-precision normalization measurements within the LHCb experiment.

---

### Key Features
- B⁺ mass fits using Crystal Ball functions
- RooFit + extended maximum likelihood
- Quantile-based time binning
- Yield ratio stability analysis across LHCb 2024 blocks
