# twinsvm 0.0.4

* Test suite no longer compares support-vector counts at a fixed tolerance.
  The count thresholds dual variables that can sit exactly at the cost bound, so
  it varies with the BLAS the checking machine links against; this made
  `test-baselines.R` fail on the CRAN MKL and OpenBLAS check flavors while
  passing everywhere else. Counts are now compared up to a small fraction of the
  sample size, and the exact-value regression test is `skip_on_cran()`.
  Per-point prediction agreement against `e1071` is unchanged and remains the
  strict correctness check.
* Reference-fit tolerances relaxed from 1e-6 to 1e-5 for the same reason.
* No model-fitting or prediction behavior changed.

# twinsvm 0.0.3

* New `resample_accuracy()` evaluates `tsvm()` and `svms()` over repeated
  stratified random holdout splits and reports the mean test accuracy together
  with its between-split standard deviation and standard error, so method
  comparisons can be judged against sampling noise. Splits are drawn once and
  shared across models for a paired comparison. Adds `print` and `plot` (boxplot)
  methods for the returned `resample_accuracy` object.
* The animation test for `morph_boundary()` now skips cleanly when the suggested
  `gganimate` package is absent, so `R CMD check` passes without Suggests
  installed. No model-fitting or prediction behavior changed.

# twinsvm 0.0.2

* Fixed cross-platform test portability on macOS by avoiding bitwise equality
  checks for generated floating-point data.
* No model-fitting or prediction behavior changed.
