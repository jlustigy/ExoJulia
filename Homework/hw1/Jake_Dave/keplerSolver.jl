#############################
#
# Authors: David Fleming (github: dflemin3) and Jake Lustig-Yaeger (github: jlustigy)
# Date: April 3rd, 2016
#
# Functions to solve Kepler's Equation
#
############################

# stest keplerEquation(0.5, 0.5)
#@stest [[keplerEquation!(M,ecc) for M in linspace(0,2pi,100)] for ecc in linspace(0,0.999,100)]

function keplerEquation!(M::Float64, e::Float64;TOL=1.0e-12,MAX_ITER=30)
    # Given mean anomaly M in radians and eccentricity e, return the eccentric anomoly
    # by solving Kepler's Equation via the Newton-Raphson Method
   
    # Make M range from [0,2pi)
    Mnew = mod(M,2.0pi)

    # If e > 1, E doesn't make sense
    if e >= 1.0
        return -1.0
    end
    
    # If e == 0, E == M
    if e < TOL
        return Mnew
    end
    
    # Initial guesses for E
    E0 = Mnew + 0.85 * e * sign(sin(Mnew)) 
    E = E0
    
    # dE = |E - E0| for convergence check
    dE = 1.0
    
    # Iteration count to make sure loop doesn't go on forever
    iter = 0
    
    # Solve!
    while dE > TOL
        E0 = E
        E = E0 - ((E0 - e * sin(E0) - Mnew)/(1.0 - e * cos(E0)))
        dE = abs(E - E0)
    
        # Has solution been running for too long?
        if iter >= MAX_ITER
            println("Error: keplerEqn reached MAX_ITER: $MAX_ITER")
            return -1
        end
        
        iter += 1
    end
        
    return E
                    
end

# end function
