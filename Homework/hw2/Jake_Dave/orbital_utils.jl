######################
#
#
#
#
######################


function mean_anomaly(P::Float64, t::Float64, t0::Float64)
    (2.0pi / P) * (t-t0)
end

function f_from_E(E::Float64, ecc::Float64)
    2.0 * atan2(sqrt(1.0 - ecc) * cos(E/2.0), sqrt(1.0 + ecc) * sin(E/2.0))
end

function r_from_f(a::Float64, ecc::Float64, f::Float64)
    a * (1.0 - ecc*ecc) / (1.0 - ecc*cos(f))
end

function f_from_M(ecc::Float64, M::Float64)
    f_from_E(ExoJulia.Orbit.kepler_solve!(M, ecc), ecc)
end

function f_from_t(P::Float64, ecc::Float64, t::Float64, tp::Float64)
    f_from_M(ecc, mean_anomaly(P, t, tp))
end 

function f_func(semi::Float64,dE::Float64,r0::Float64)
    (semi/r0)*(cos(dE) + 1.0)
end

function g_func(t::Float64,t0::Float64,P::Float64,dE::Float64)
    (t-t0) + (P/(2.0*pi))*(sin(dE) - dE) 
end

function f_dot_func(semi::Float64,n::Float64,r0::Float64,r::Float64,dE::Float64)
    ((-a*a*n)/(r*r0))*sin(dE)
end

function g_dot_func(semi::Float64,r::Float64,dE::Float64)
    (a/r)*(cos(dE)-1.0) + 1.0
end