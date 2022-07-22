function optiminv_T(λ, M, HH, H, α, β, σ, η, μ, ϕ, τ_w, τ_e, t, a, init_s, init_e)

    s_value=μ*ϕ/(μ*ϕ+σ/β-η)
    e_value=( (1+τ_e)/((1-t)*(1-τ_w)*λ*((2*H/M)^β)*(a^(α*β/σ))*(s_value^(ϕ*β/σ))*(M/2/HH)*η*β/σ) )^(1/(η*β/σ-1))
    obj_value=log(1-s_value)+μ*log((1-t)*(1-τ_w)*λ*((2*H/M)^β)*(a^(α*β/σ))*(s_value^(ϕ*β/σ))*(e_value^(η*β/σ))*(M/2/HH)-(1+τ_e)*e_value)

    return obj_value,s_value,e_value

end
