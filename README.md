# AstroPDM.jl

Native Julia phase dispersion minimization (PDM) for unevenly sampled astronomical
light curves. No C compiler, shared library, build step, or runtime dependencies
are required. Julia 1.8 or later is supported.

From this checkout, run `julia --project=.` and load `using AstroPDM`.
To use the checkout in another Julia environment:

```julia
using Pkg
Pkg.develop(path="/path/to/AstroPDM.jl")
using AstroPDM
```

```julia
using AstroPDM

t = collect(0.0:0.13:40.0)
y = sin.(2pi .* t ./ 2.7)
frequencies, theta = pdm(t, y; f_min=0.2, f_max=0.6, delf=0.001,
                         nbin=10, covers=3)
best_period = 1 / frequencies[argmin(theta)]

# An explicit frequency grid, preserving the supplied order:
frequencies, theta = pdm(t, y, [0.35, 0.37, 0.39])

# Optional measurement uncertainties:
sigma = fill(0.05, length(t))
frequencies, theta = pdm(t, y; f_min=0.2, f_max=0.6, delf=0.001, sigma)

# Former wrapper's positional call signature remains available:
frequencies, theta = pdm2(t, y, Float64[], 0.2, 0.6, 0.001, 10)
```

PDM folds observations at each trial frequency and compares the scatter within
phase bins with the scatter in the full light curve. A small theta indicates a
candidate period; aliases and harmonics can also produce minima. Frequency is
in inverse units of the input times. The search grid is `f_min:delf:f_max`;
`f_max` is included only if it lies on that grid.

For unweighted observations, the implemented statistic is

```math
\Theta = \frac{\sum_{c,b}\sum_{i\in(c,b)}(y_i-\bar y_{c,b})^2}
                   {\sum_{c,b}(n_{c,b}-1)}
          \bigg/ \frac{\sum_i(y_i-\bar y)^2}{N-1},
```

where the bin sums include only occupied bins. Each of `covers` partitions the
unit phase interval into `nbin` equal bins, shifting the origin by
`c / (nbin * covers)`. Bins start at phase zero in the first cover. Phases use
the earliest time as origin. Empty bins are ignored and singleton bins have
zero residual degrees of freedom. A trial with no residual degrees of freedom
returns `NaN`; exclude such trials before finding a minimum. Zero frequency is
allowed and gives theta equal to one.

Positive `sigma` selects inverse-variance weighted bin and global means and
weighted residual sums, retaining count-based degrees of freedom. This is a
weighted extension of the unweighted statistic, not a significance probability.
Multiplying all uncertainties by a common factor leaves theta unchanged.
Omitted, empty, or all-zero uncertainties select the unweighted statistic.
Mixed zero and positive uncertainties are rejected.

Inputs may be real-valued vectors, ranges, or views and are converted to
Float64. At least three finite observations, a nonzero time span, and
nonconstant values are required. Invalid dimensions and parameters throw
exceptions. Inputs remain unchanged and returned arrays belong to the caller;
independent calls share no mutable state. Runtime is proportional to the number
of frequencies times covers times (observations + bins); workspace is linear
in observations, bins, and frequencies.

## Migration from PDM

The package/module is now `AstroPDM` and retains the original package UUID.
Replace `using PDM` with `using AstroPDM`. The `pdm2` positional signature and
pair of result arrays are retained. `test_pdm` was a debugging helper and has
been replaced by the automated suite.

Results intentionally differ from the old C PDM2 variant: this implementation
uses the classical pooled-bin statistic and explicit shifted covers. It does
not apply automatic gap segmentation, uncertainty subtraction/clipping,
spline fits, or subharmonic averaging. The C wrapper had fixed array limits,
unsafe shared result buffers, input mutation, and a mismatch between Julia's
zero-offset pointers and C's one-based data accesses. Those behaviors are not
reproduced. The old 100-observation minimum and 100-bin maximum are removed.

The initial native release is version 0.1.0. Legacy C/Python sources and
historical notebooks are excluded from this repository. The package has no
runtime dependencies; local Manifest files are ignored.

## Validation and references

Run `julia --project=. -e 'using Pkg; Pkg.test()'`.
Tests compare against an independent two-pass definition, check weighted and
multiple-cover statistics, recover a synthetic period, exercise invalid and
sparse inputs, and scan the bundled Tu Cas observations.

- [Stellingwerf (1978), Period determination using phase dispersion minimization](https://scixplorer.org/abs/1978ApJ...224..953S/abstract)
- [Stingray's PDM documentation and statistic](https://docs.stingray.science/en/v2.2.3/notebooks/Pulsar/Phase%20Dispersion%20Minimization.html)
