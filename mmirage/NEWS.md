# mmirage 0.2.0

- Publish the installable package alongside the repository scripts.
- Add automatic EM termination metrics, stopping reasons, initial parameters,
  objective and likelihood change traces, fixed-block metadata, boundary flags,
  simplex/box stationarity residuals, and diagnostic messages.
- Add `mirage_convergence()` with optional numerical observed information,
  eigenvalues and conditioning. Boundary nonregularity and numerical failures
  are reported explicitly; no automatic standard errors are produced.
- Preserve the existing EM updates, stopping definition, likelihood, priors,
  BFs, parameters and posterior results. Add an exact regression fixture from
  the previous core implementation.
- Publish the existing separate optional native-simplex ML refinement and its
  stored synthetic tests. It does not replace EM or claim a global joint optimum.
- Retain remote main's expanded benchmark metric utilities in the package.
- Add current manuscript aggregate summaries with convergence, calibration and
  input-harmonization limitations, keeping historical outputs unchanged.
- Add repository and package tests to GitHub Actions.
