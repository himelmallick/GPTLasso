# GPTLasso 0.99.0

- Supports both the exact legacy PCL list contract and `MultiAssayExperiment`
  across fitting, cross-validation, and prediction, including mixed-format
  train/test workflows.
- Returns PCL-style containers from simulation helpers and packaged examples.
- Supports multiple- and single-study data with one or multiple views through
  stable `gptLasso` and `cv.gptLasso` result classes.
- Adds Gaussian, binomial, and Cox workflows, including study-aware Cox folds
  and prediction semantics.
- Adds an evaluated BiocStyle vignette and a testthat test suite.
