# Implementing the radial velocity/doppler formula for
# a planet with an elliptical orbit
# Read in data table for the mystery planet with columns:
#   time
#   RV (m/s)
#   error (m/s)
include("../../hw1/nicole/nicole_hw1_solver.jl")
using nicole_hw1_solver
using DataFrames
using PyPlot
using LsqFit

# Functions
function get_period(time,RV)
    t_N = 5000
    per_array = linspace(10.0,500.0,t_N)
    chi = zeros(t_N)
    for i in range(1,t_N)
        per_rems  = mod(time,per_array[i])
        per_where = sortperm(per_rems)
        rv_rems   = RV[per_where]
        #println(i)
        sum = 0.0
        for j in range(2,length(time)-2)
            sum += (rv_rems[j]-rv_rems[j-1])^2.
            #println(j)
        end
        chi[i] = sum
    end
    return per_array[indmin(chi)]
end

function get_M(per,time,t_peri)
    M   = 2 * pi/per *(time-t_peri)
    return M
end

function get_f(ecc,E)
    f   = 2.0 * atan(((1.0 + ecc)/(1.0 - ecc))^0.5 * tan(E/2.0))
    return f
    end
    
function get_vrad(time,init_param)
    vrad = zeros(length(time))
    ecc = init_param[2]
    if (ecc >= 1) || (ecc < 0)
    vrad = 0
    end
    for k in range(1,length(time)-1)
        #println(k)
        M = get_M(init_param[1],time[k],init_param[3])
        E = nicole_hw1_solver.solver_loop(M,ecc)
        f = get_f(ecc,E) 
        vrad[k] = init_param[4]*cos(f) + init_param[5]*sin(f) + init_param[6]
        #println(vrad[k])
    end
    return vrad
end