#############################
#
# Authors: David Fleming (github: dflemin3) and Jake Lustig-Yaeger (github: jlustigy)
# Date: April 3rd, 2016
#
# Functions to solve Kepler's Equation
#
############################

#@stest keplerEquation(0.5, 0.5)

function keplerEquation(M, e)
    # Given mean anomaly M in radians and eccentricity e, return the eccentric anomoly
    # by solving Kepler's Equation via the Newton-Raphson Method
    
    const TOL = 1.0e-12 # Solution tolerance
    const MAX_ITER = 30
    
    # If e > 1, E doesn't make sense
    if e >= 1.0
        return -1.0
    end
    
    # If e == 0, E == M
    if e < TOL
        return M
    end
    
    # Initial guesses for E
    E0 = M .+ 0.85.*e .* sign(M) 
    E = E0
    
    # dE = |E - E0| for convergence check
    dE = TOL + 1.0
    
    # Iteration count to make sure loop doesn't go on forever
    iter = 0
    
    # Solve!
    while dE > TOL
        E = E0 .- (E0 .- e .* sin(E0) .- M)./(1.0 .- e .* cos(E0))
        dE = abs(E .- E0)
        E0 = E
    
        # Has solution been running for too long?
        if iter >= MAX_ITER
            println("Error: keplerEqn reached MAX_ITER: $MAX_ITER")
            break
        end
        
        iter += 1
    end
        
    return E
                    
end

# end function