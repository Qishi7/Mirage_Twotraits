# Mirage_Twotraits

This repository hosts the MIRAGE Two Trait Extension: a two-trait extension of MIRAGE for joint gene-level rare-variant analysis. The motivating first use case is autism as trait 1 and schizophrenia as trait 2.

Website: https://qishi7.github.io/Mirage_Twotraits/

Repository: https://github.com/Qishi7/Mirage_Twotraits

## Model Direction

The first implementation uses a gene-level four-state model:

- `00`: neither trait
- `10`: trait 1 only
- `01`: trait 2 only
- `11`: both traits

State `11` represents gene-level pleiotropy. It does not require the same variants to affect both traits. Conditional on gene state, the two trait-specific likelihoods are independent products of trait-specific MIRAGE likelihoods.

## Repository Contents

- `R/`: first R implementation of the two-trait MIRAGE EM model and toy simulation helper.
- `tests/testthat/`: unit tests and synthetic simulation smoke checks.
- `analysis/`: workflowr source pages.
- `docs/`: rendered GitHub Pages website.
- `Mirage_Twotraits.Rproj`: RStudio project file.

## Data Privacy

No real ASD, schizophrenia, WES, dbGaP, SPARK, protected, private, or large data files should be committed. Tests and examples use toy or simulated data only. Local reference documents such as papers, technical specifications, TODO files, VCF exploration reports, and implementation plans are excluded from version control.

## Running Tests

From the repository root, run:

```powershell
Rscript --vanilla tests/testthat.R
```

On the development Windows machine used for the first build, the explicit local command was:

```powershell
& 'D:\software\R-4.4.2\bin\Rscript.exe' --vanilla tests/testthat.R
```
