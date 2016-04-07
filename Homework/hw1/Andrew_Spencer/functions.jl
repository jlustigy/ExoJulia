# Kepler's equation
# ecc = eccentricity
# E = eccentric anomaly
# M = mean motion  "$HW" == "./$1"
#@stest time_newt()
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

function time_newt()
    eps = 1e-14

    #Timing test
    #Construct parameters
    ecc = linspace(0,0.999,100)
    M = linspace(0,2*pi,100)

    #Calculate eccentric anomaly E
    for (i,valecc) in enumerate(ecc)
        for (j,valM) in enumerate(M)
            E = valM + 0.85*valecc*sign(sin(valM))
            newt_kepler(valecc,E,valM,eps)
        end
    end
end
