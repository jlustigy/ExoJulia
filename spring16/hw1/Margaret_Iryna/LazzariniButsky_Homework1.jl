#@stest [[Kepler(M,ecc) for M in linspace(0,20pi,1000)] for ecc in linspace(0,0.999,100)]

function Kepler(M, e)
    E_old = M + 0.85 * e
    if(sin(M) < 0)
      E_old = M - 0.85 * e
    end
    epsilon = 1.0
    while epsilon > 1.0E-12
        g = E_old - e*sin(E_old) - M
        g_prime = 1.0 - e*cos(E_old)
        E_new = E_old - (g/g_prime)
        epsilon = abs(E_old - E_new)/E_old
        E_old = E_new
    end
    return E_old
end
