# Kepler's equation
# ecc = eccentricity
# E = eccentric anomaly
# M = mean motion  "$HW" == "./$1"
#@stest [[time_newt(ecc,M) for M in linspace(0.0,2*pi,100)] for ecc in linspace(0.0,0.999,100)]
function g(ecc,E,M)
    E-ecc*sin(E) - M
end

# Derivative of Kepler's equation
function dg_dE(ecc,E)
    1-ecc*cos(E)
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

function time_newt(ecc,M)
    eps = 1e-15
    E0 = M + 0.85*ecc*sign(sin(M))
    newt_kepler(ecc,E0,M,eps)
end
