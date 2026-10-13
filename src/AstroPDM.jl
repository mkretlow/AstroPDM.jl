module AstroPDM

export pdm, pdm2

function prepare(t, y, sigma)
    length(t) == length(y) || throw(DimensionMismatch("t and y must have equal lengths"))
    length(t) >= 3 || throw(ArgumentError("at least three observations are required"))
    tt, yy = Float64.(collect(t)), Float64.(collect(y))
    all(isfinite, tt) && all(isfinite, yy) || throw(ArgumentError("observations must be finite"))
    tt .-= minimum(tt)
    all(isfinite, tt) || throw(ArgumentError("time span overflows Float64"))
    maximum(tt) > 0 || throw(ArgumentError("times must span a nonzero interval"))
    weights = ones(length(yy))
    if sigma !== nothing && !isempty(sigma)
        length(sigma) == length(yy) || throw(DimensionMismatch("sigma must match t and y"))
        ss = Float64.(collect(sigma))
        all(isfinite, ss) && all(>=(0), ss) || throw(ArgumentError("sigma must be finite and nonnegative"))
        if !all(iszero, ss)
            all(>(0), ss) || throw(ArgumentError("sigma must be all zero or strictly positive"))
            weights .= (minimum(ss) ./ ss).^2
            all(>(0), weights) || throw(ArgumentError("sigma dynamic range is too large"))
        end
    end
    # Center and scale before accumulating variances to avoid cancellation.
    yy .-= yy[1]
    all(isfinite, yy) || throw(ArgumentError("value span overflows Float64"))
    scale = maximum(abs, yy)
    scale > 0 || throw(ArgumentError("PDM is undefined for constant values"))
    yy ./= scale
    yy .-= sum(weights .* yy) / sum(weights)
    variance = sum(weights .* yy.^2) / (length(yy) - 1)
    variance > 0 && isfinite(variance) || throw(ArgumentError("PDM requires positive finite variance"))
    return tt, yy, weights, variance
end

"""
    pdm(t, y, frequencies; nbin=10, covers=1, sigma=nothing)

Compute phase dispersion minimization at explicit trial frequencies, returning
`(frequencies, theta)` as independent `Vector{Float64}` arrays. Smaller theta
indicates a tighter folded light curve. Frequencies are in inverse time units.

Each cover partitions phase into `nbin` equal bins; successive covers shift the
origin by `1 / (nbin * covers)`. Residual sums and degrees of freedom are pooled
across covers. Empty bins are ignored; singleton bins contribute zero degrees
of freedom. A trial with no residual degrees of freedom returns `NaN`.

With no uncertainties, theta is pooled within-bin sample variance divided by
total sample variance. Positive `sigma` uses inverse-variance weighted means
and residuals, with count-based degrees of freedom in both variances. Empty or
all-zero sigma selects unweighted PDM. Mixed zero/positive sigma is invalid.
Inputs must be finite, contain at least three points and have nonconstant values
and times. Inputs are never mutated. Zero frequency is allowed (theta = 1).
"""
function pdm(t::AbstractVector{<:Real}, y::AbstractVector{<:Real},
             frequencies::AbstractVector{<:Real}; nbin::Integer=10,
             covers::Integer=1, sigma::Union{Nothing,AbstractVector{<:Real}}=nothing)
    nbin >= 2 || throw(ArgumentError("nbin must be at least two"))
    covers >= 1 || throw(ArgumentError("covers must be positive"))
    ff = Float64.(collect(frequencies))
    all(isfinite, ff) && all(>=(0), ff) || throw(ArgumentError("frequencies must be finite and nonnegative"))
    tt, yy, weights, variance = prepare(t, y, sigma)
    isempty(ff) || isfinite(maximum(tt) * maximum(ff)) || throw(ArgumentError("phase calculation overflows Float64"))
    counts = zeros(Int, nbin)
    means, sums = zeros(nbin), zeros(nbin)
    bins = Vector{Int}(undef, length(tt))
    phases = similar(tt)
    theta = similar(ff)
    for (k, f) in pairs(ff)
        phases .= mod.(tt .* f, 1.0)
        residual, dof = 0.0, 0
        for cover in 0:(covers - 1)
            fill!(counts, 0); fill!(means, 0); fill!(sums, 0)
            shift = (cover / covers) / nbin
            for i in eachindex(tt)
                b = min(floor(Int, mod(phases[i] + shift, 1.0) * nbin) + 1, nbin)
                bins[i] = b
                counts[b] += 1
                sums[b] += weights[i]
                means[b] += weights[i] * yy[i]
            end
            for b in eachindex(counts)
                if counts[b] > 0
                    means[b] /= sums[b]
                    dof += counts[b] - 1
                end
            end
            for i in eachindex(yy)
                residual += weights[i] * (yy[i] - means[bins[i]])^2
            end
        end
        theta[k] = dof > 0 ? residual / dof / variance : NaN
    end
    return ff, theta
end

"""
    pdm(t, y; f_min=0.0, f_max=1.0, delf=0.1, nbin=10, covers=1, sigma=nothing)

Scan `f_min:delf:f_max`. The upper endpoint is included only when it lies on
the grid. See the explicit-frequency method for statistic and input semantics.
"""
function pdm(t::AbstractVector{<:Real}, y::AbstractVector{<:Real};
             f_min::Real=0.0, f_max::Real=1.0, delf::Real=0.1, kwargs...)
    lo, hi, step = Float64(f_min), Float64(f_max), Float64(delf)
    all(isfinite, (lo, hi, step)) || throw(ArgumentError("frequency bounds and step must be finite"))
    0 <= lo <= hi || throw(ArgumentError("require 0 <= f_min <= f_max"))
    step > 0 || throw(ArgumentError("delf must be positive"))
    return pdm(t, y, collect(lo:step:hi); kwargs...)
end

"""
    pdm2(t, y, s=Float64[], f_min=0.0, f_max=1.0, delf=0.1, nbin=10; covers=1)

Compatibility call signature for the former C wrapper. Computes native PDM;
does not reproduce legacy segmentation, sigma clipping or subharmonic smoothing.
"""
function pdm2(t::AbstractVector{<:Real}, y::AbstractVector{<:Real},
              s::AbstractVector{<:Real}=Float64[], f_min::Real=0.0,
              f_max::Real=1.0, delf::Real=0.1, nbin::Integer=10; covers::Integer=1)
    return pdm(t, y; f_min, f_max, delf, nbin, covers, sigma=s)
end

end
