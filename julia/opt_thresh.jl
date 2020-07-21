function opt_thresh(a_grid,vT,aO,vO)
    # 'vO' is a scalar, 'vT' and 'a_grid' are one-dimensional arrays; 'aO' is an initial guess (used in 'find_zero').

    # Fit spline to 'a_grid' and 'vT':
    spl = Spline1D(a_grid,vT)
    # Solve for the threshold 'aT' such that 'spl(aT)' = 'vO':
    function f(x)
        if x <= a_grid[1]
            # Non-linear extrapolation to guarantee existence of an inverse function (to characterize a threshold of a^O as a function of a^T):
            # Functional form: V = α + βa^γ
            α = -1e6
            γ = derivative(spl,a_grid[1])*a_grid[1]/(vT[1]-α)
            β = derivative(spl,a_grid[1])/(γ*a_grid[1]^(γ-1))
            vO - (α+β*x^γ)
        elseif x >= a_grid[end]
            mu = (x-a_grid[end-1])/(a_grid[end]-a_grid[end-1])
            vO - (mu*vT[end] + (1-mu)*vT[end-1])
        else
            vO - spl(x)
        end
    end
    r = find_zero(f,(0,a_grid[end])) # bisection method
    return r
end
