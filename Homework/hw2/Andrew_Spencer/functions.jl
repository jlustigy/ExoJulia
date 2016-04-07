# Kepler's equation
# ecc = eccentricity
# E = eccentric anomaly
# M = mean motion
function g(ecc,E,M)
    retval = E-ecc*sin(E) - M
    retval
end

# Derivative of Kepler's equation
function dg_dE(ecc,E)
    retval = 1-ecc*cos(E)
    retval
end

# Newton-raphson function
function newt_kepler(ecc,E,M,eps)
    h = -g(ecc,E,M)/dg_dE(ecc,E)
    E += h
    if(abs(h) > eps)
        newt_kepler(ecc,E,M,eps)
    else
        E
    end
end