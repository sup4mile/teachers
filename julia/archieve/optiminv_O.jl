function optiminv_O(A, M, H, α, σ, η, μ, ϕ, τ_w, τ_e, t, a, init_s, init_e)

    s_value=μ*ϕ/(μ*ϕ+1-η)
    e_value=( (1+τ_e)/((1-t)*(1-τ_w)*A*((2*H/M)^σ)*(a^α)*(s_value^ϕ)*η) )^(1/(η-1))
    obj_value=log(1-s_value)+μ*log((1-t)*(1-τ_w)*A*((2*H/M)^σ)*(a^α)*(s_value^ϕ)*(e_value^η)-(1+τ_e)*e_value)

    return obj_value,s_value,e_value
end
