######################
#
#
#
#
######################


function mean_anomaly!(P::Float64, t::Float64, t0::Float64)
    (2.0pi / P) * (t-t0)
end

function f_from_E!(E::Float64, ecc::Float64)
    2.0 * atan2(sqrt(1.0 - ecc) * cos(E/2.0), sqrt(1.0 + e) * sin(E/2.0))
end

function r_from_f!(a::Float64, ecc::Float64, f::Float64)
    a * (1.0 - ecc*ecc) / (1.0 - ecc*cos(f))
end

function f_from_M!(ecc::Float64, M::Float64)
    f_from_E!(ExoJulia.Orbit.kepler_solve!(M, ecc), ecc)
end
