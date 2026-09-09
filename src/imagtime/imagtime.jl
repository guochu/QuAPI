include("util.jl")
include("imagcorr.jl")
include("fermionic.jl")
include("bosonic.jl")
include("../deprecated/imagtime/fermionic.jl")
include("../deprecated/imagtime/bosonic.jl")


"""
    Δiw_to_Δτ(Δiw::AbstractVector{<:Number}; β::Real, N::Int)

Convert a hybridization function given in the Matsubara-frequency domain to
the discrete imaginary-time hybridization matrix Δ[τᵢ, τⱼ].

`Δiw` holds the values Δ(iωₙ) at the first `length(Δiw)` Matsubara
frequencies (an even number, ordered as returned by `ifrequencies`).  The
discrete Fourier transform

    Δ(τ) = Σₙ -(2/β) Δ(iωₙ) (1 - cos(ωₙ δτ)) e^{-iωₙ k δτ} / ωₙ²

is evaluated on a grid of `N` steps δτ = β/N (detailed derivation in the
author's notes, `docs/quapi.pdf`, not distributed with the repository),
returning an `ImagCorrelationFunction` wrapping a real
`CorrelationMatrix`.  The imaginary part is checked to vanish to numerical
precision, and an `error` is raised otherwise.

Arguments:
- `Δiw`:  hybridization values at the Matsubara frequencies (even length);
- `β`:  the inverse temperature;
- `N`:  number of imaginary-time steps, δτ = β/N.
"""
function Δiw_to_Δτ(Δiw::AbstractVector{<:Number}; β::Real, N::Int)
    iseven(length(Δiw)) || throw(ArgumentError("even number of frequencies expected"))
    δτ = β/N
    nmax = div(length(Δiw), 2) - 1

    ηⱼₖ = zeros(ComplexF64, N)
    ηₖⱼ = zeros(ComplexF64, N)

    for (i, ωₙ) in enumerate(ifrequencies(β, nmax)) 
        a = -(2/(β*ωₙ^2)) * Δiw[i] * (1 - cos(ωₙ * δτ))
        ηⱼₖ[1] += a / 2
        ηₖⱼ[1] += a / 2
        for k in 1:N-1
            ηⱼₖ[k+1] += a * exp(-im*ωₙ*k*δτ)
            ηₖⱼ[k+1] += a * exp(im*ωₙ*k*δτ)
        end
    end
    for i in 1:N
        (abs(imag(ηⱼₖ[i])) < 1.0e-8) || error("imag part of ηⱼₖ[$i]=$(imag(ηⱼₖ[i])) is nonzero")
        (abs(imag(ηₖⱼ[i])) < 1.0e-8) || error("imag part of ηₖⱼ[$i]=$(imag(ηₖⱼ[i])) is nonzero")
    end

    ImagCorrelationFunction(CorrelationMatrix{Float64}(real(ηⱼₖ), real(ηₖⱼ)))
end