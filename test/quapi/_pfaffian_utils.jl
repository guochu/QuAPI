# 共享 Pfaffian/Grassmann 求积工具（由 quapi/bcsreal.jl 与 test/dense 对比测试共用）

# =====================================================================
# Pfaffian of an antisymmetric matrix (Parlett–Reid, O(n³))
# =====================================================================
function pfaffian(A::AbstractMatrix)
    n = size(A, 1)
    (iseven(n)) || throw(ArgumentError("Pfaffian requires even dimension"))
    n == 0 && return one(eltype(A))
    B = Matrix{ComplexF64}(A)
    p = 1.0 + 0.0im
    for k in 1:2:n-1
        if iszero(B[k, k+1])
            found = false
            for i in k+2:n
                if !iszero(B[k, i])
                    # swap rows/columns i <-> k+1 (each swap flips Pf by -1)
                    B[:, [k+1, i]] .= B[:, [i, k+1]]
                    B[[k+1, i], :] .= B[[i, k+1], :]
                    p = -p
                    found = true
                    break
                end
            end
            found || return zero(ComplexF64)
        end
        p *= B[k, k+1]
        for i in k+2:n, j in i:n
            B[i, j] = (B[i, j] * B[k, k+1] - B[i, k] * B[j, k+1] + B[i, k+1] * B[j, k]) / B[k, k+1]
            B[j, i] = -B[i, j]
        end
    end
    return p
end

# Gaussian Grassmann integral with monomial insertions:
#   ∫ exp(½ξᵀBξ) ξ_{i₁}···ξ_{i_{2k}} = Pf(B) · Pf(C),  C_{rs} = (B⁻¹)_{iᵣ,iₛ}
function gauss_integral(B::AbstractMatrix, monomial::Vector{Int})
    m = length(monomial)
    iseven(m) || return 0.0 + 0.0im
    m == 0 && return pfaffian(B)
    Binv = inv(B)
    C = Binv[monomial, monomial]
    return pfaffian(B) * pfaffian(C)
end

# =====================================================================
