include("nicole_hw1_solver.jl")
using nicole_hw1_solver

# Testing that importing nicole_hw1_solver works
#M   = pi/6
#ecc = 0.6
#E   = 2.

#println(nicole_hw1_solver.next(E,M,ecc))

# Test my final E value by putting it into our original equation
M   = 5.58505
ecc = 1.0
#E   = 2.

E_final = nicole_hw1_solver.solver_loop(M,ecc)
println(E_final)
E_test = M + (ecc * sin(E_final))
println(E_test)

if (E_final) == (E_test)
    println("Initial test case: These values are equal; the solver worked!")
    println(" ")
end

# Test to see if this works for 0 < e < 1 and 0 < M < 2pi
#ecc = collect(0:0.1:1) #This is a neat tool, but not using here
ecc = linspace(0.,1.,50)
println("Test eccentricity, e, range: ",collect(ecc))

M = linspace(0,(2*pi),50)
println("Test mean anomaly, M, range: ",collect(M))

count = 0
for i in eachindex(ecc)
    for j in eachindex(M)
        if j == length(M)
            break
        end
        E_final = nicole_hw1_solver.solver_loop(M[j],ecc[i])
        E_test = M[j] + (ecc[i] * sin(E_final))
        if (round(E_final,12)) == (round(E_test,12))
            #println("These values are equal; the solver worked!")
            #count += 1
            #println(count)
        else
            #println("DIDN'T WORK.")
            #println("M: ",M[j])
            #println("ecc: ",ecc[i])
            #println("Final E: ",E_final)
            #println("Test E: ",E_test)
            count += 1
            #println(count)
        end
    end
end
println(count/2500)

count = 0
for i in eachindex(ecc)
    for j in eachindex(M)
        if j == length(M)
            break
        end
        E_final = nicole_hw1_solver.solver_loop(M[j],ecc[i])
        E_test = M[j] + (ecc[i] * sin(E_final))
        if (E_final) == (E_test)
            #println("These values are equal; the solver worked!")
            #count += 1
            #println(count)
        else
            #println("DIDN'T WORK.")
            #println("M: ",M[j])
            #println("ecc: ",ecc[i])
            #println("Final E: ",E_final)
            #println("Test E: ",E_test)
            count += 1
            #println(count)
        end
    end
end
println(count/2500)

println("When truncated to 12 numbers after the decimal, all (except < 0.5%) M and e values result in a solution to Kepler's equations. Also, excluding M = 0.0 and e = 1.0. Otherwise, ~ 20% of the M and e pairs break (including e = 1.0 and M = 6.28319).")