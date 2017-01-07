module nicole_hw1_solver 

# Solving kepler's equations
# E   = eccentric anomaly
# ecc = eccentricity 
# M   = mean anomaly

# Testing stuff
#M   = pi/24
#ecc = 0.1

# Functions
init(M) = M
next(E,M,ecc) = M + (ecc * sin(init(M)))
g(E,M,ecc)    = E - M - (ecc * sin(E))  # where g_E = 0 
dg(E,ecc)   = 1 - (ecc * cos(E))

function solver_loop(M,ecc)
    E_init = init(M)
    E_next = next(E_init,M,ecc)
    i = 0
    while i <= 5
        E_prev = E_next
        E_next = E_prev - ( g(E_next,M,ecc) / (dg(E_next,ecc)) )
        #println(E_next)
        i += 1
    end
    E_final = E_next

end
    
end