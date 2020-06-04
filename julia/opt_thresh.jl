function opt_thresh(a_grid,vT,aO,vO)
    # 'vO' is a scalar, 'vT' and 'a_grid' are one-dimensional arrays; 'aO' is an initial guess (used in 'find_zero').

    # Fit spline to 'a_grid' and 'vT':
    spl = Spline1D(a_grid,vT)
    # Solve for the threshold 'aT' such that 'spl(aT)' = 'vO':
    function f(x)
        if x <= a_grid[1]
            # Linear extrapolation such that ability threshold is strictly positive:
            mu_lo = (-a_grid[1])/(a_grid[2]-a_grid[1])
            vT_lo = min(mu_lo*vT[2] + (1-mu_lo)*vT[1],vO)
            mu = x/a_grid[1]
            vO - (mu*vT[1] + (1-mu)*vT_lo)
        elseif x >= a_grid[end]
            mu = (x-a_grid[end-1])/(a_grid[end]-a_grid[end-1])
            vO - (mu*vT[end] + (1-mu)*vT[end-1])

        else
            vO - spl(x)
        end
    end

    r = find_zero(f,aO)
    return r
end
