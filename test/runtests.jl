using AstroPDM, Test, Random

# Independent, two-pass definition for one cover.
function reference(t, y, f; nbin=10, covers=1, sigma=ones(length(y)))
    w = 1 ./ sigma.^2
    mean = sum(w .* y) / sum(w)
    total = sum(w .* (y .- mean).^2) / (length(y) - 1)
    rss, dof = 0.0, 0
    phase = mod.((t .- minimum(t)) .* f, 1)
    for c in 0:covers-1
        labels = floor.(Int, mod.(phase .+ c / (covers * nbin), 1) .* nbin)
        for b in 0:nbin-1
            ii = findall(==(b), labels)
            isempty(ii) && continue
            mu = sum(w[ii] .* y[ii]) / sum(w[ii])
            rss += sum(w[ii] .* (y[ii] .- mu).^2)
            dof += length(ii) - 1
        end
    end
    dof == 0 ? NaN : rss / dof / total
end

@testset "Native statistic" begin
    rng = MersenneTwister(42)
    t = 40 .* rand(rng, 250)
    y = sin.(2pi .* 0.37 .* t) .+ 0.1 .* randn(rng, length(t))
    sigma = 0.1 .+ rand(rng, length(t))
    for covers in (1, 3), s in (ones(length(t)), sigma)
        f, theta = pdm(t, y, [0.0, 0.17, 0.37, 0.9]; covers, sigma=s)
        @test theta ≈ [reference(t, y, fi; covers, sigma=s) for fi in f]
    end
    f, theta = pdm(t, y; f_min=0.2, f_max=0.6, delf=0.001, covers=3)
    @test abs(f[argmin(theta)] - 0.37) <= 0.002
    @test minimum(theta) < 0.1
    @test pdm(t, y, [0.0])[2][1] ≈ 1
    original_t, original_y = copy(t), copy(y)
    old_f, old_theta = pdm2(t, y, Float64[], 0.2, 0.6, 0.001)
    @test (old_f, old_theta) == pdm(t, y; f_min=0.2, f_max=0.6, delf=0.001)
    @test t == original_t && y == original_y
    old_copy = copy(old_theta)
    pdm(t, y, [0.3])
    @test old_theta == old_copy
    perm = randperm(rng, length(t))
    @test pdm(t[perm], y[perm], f)[2] ≈ pdm(t, y, f)[2]
    @test pdm(t .+ 1000, 5 .* y .+ 100, f)[2] ≈ pdm(t, y, f)[2]
    @test pdm(view(t, :), Float32.(y), [0.37])[2][1] ≈ pdm(t, y, [0.37])[2][1] rtol=1e-6
    @test pdm(t, y, Float64[]) == (Float64[], Float64[])
    @test pdm(t, y, [0.3]; sigma=zeros(length(t))) == pdm(t, y, [0.3])
    @test length(pdm(t, y; f_min=0.1, f_max=0.3, delf=0.1)[1]) == 3
    @test pdm(t, y; f_min=0.3, f_max=0.3)[1] == [0.3]
end

@testset "Sparse bins and boundaries" begin
    @test isnan(pdm([0, 1, 2], [0, 1, 0], [0.1]; nbin=10)[2][1])
    t, y = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75], [1, 2, 3, 4, 1, 2, 3, 4]
    @test pdm(t, y, [1]; nbin=4)[2] == [0.0]
    @test pdm([0, 0, 1, 1], [0, 1, 2, 3], [0.5]; nbin=2)[2][1] ≈ 0.3
end

@testset "Invalid input" begin
    t, y = [0., 1., 2., 3.], [0., 1., 0., 2.]
    @test_throws DimensionMismatch pdm(t, y[1:3], [1])
    @test_throws DimensionMismatch pdm(t, y, [1]; sigma=[1.])
    @test_throws ArgumentError pdm(t[1:2], y[1:2], [1])
    @test_throws ArgumentError pdm(t, ones(4), [1])
    @test_throws ArgumentError pdm(ones(4), y, [1])
    for bad in (NaN, Inf)
        @test_throws ArgumentError pdm([0, 1, 2, bad], y, [1])
        @test_throws ArgumentError pdm(t, [0, 1, 2, bad], [1])
        @test_throws ArgumentError pdm(t, y, [bad])
    end
    @test_throws ArgumentError pdm(t, y, [-1])
    @test_throws ArgumentError pdm(t, y, [1]; sigma=[0, 1, 1, 1])
    @test_throws ArgumentError pdm(t, y, [1]; sigma=[-1, 1, 1, 1])
    @test_throws ArgumentError pdm(t, y, [1]; nbin=1)
    @test_throws ArgumentError pdm(t, y, [1]; covers=0)
    @test_throws ArgumentError pdm(t, y; delf=0)
    @test_throws ArgumentError pdm(t, y; f_min=2, f_max=1)
    @test_throws ArgumentError pdm(t, y; f_min=-1)
    @test_throws ArgumentError pdm(t, y; f_max=Inf)
end

@testset "Bundled Tu Cas observations" begin
    rows = readlines(joinpath(@__DIR__, "TuCas.dat"))[2:end]
    data = [parse.(Float64, split(row)) for row in rows if !isempty(strip(row))]
    t, y = getindex.(data, 1), getindex.(data, 2)
    f, theta = pdm(t, y; f_min=0.4, f_max=0.6, delf=0.0001, nbin=5, covers=3)
    @test all(isfinite, theta)
    @test 0.45 < f[argmin(theta)] < 0.5
end
