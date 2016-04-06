#############################
#
# Authors: David Fleming (github: dflemin3) and Jake Lustig-Yaeger (github: jlustigy)
# Date: April 3rd, 2016
#
# Functions test Kepler's Equation solution
#
############################

using PyPlot

include("keplerSolver.jl")

function testKeplerEquation1()
    # This function tests the ExoJulia keplerEquation solver.
    # Kepler's Equations is solved for e ranging from [0, 0.9]
    # over all values of mean anomaly M.  The resulting eccentric anomaly E
    # is then plotted vs M 
    
    # Inits
    len = 9
    len_M = 100
    e = linspace(0,0.9,len);
    M = linspace(0,2pi,len_M);
    E = zeros(M);
    ecc = 0.0
    
    # Loop over eccentricities
    for i=1:len
        # Loop over mean anomalies
        for j=1:len_M
            E[j] = keplerEquation(M[j], e[i])
        end

        # Plot result for single eccentricity
        ecc = round(e[i],1)
        plot(M, E, lw=2, label="e = $ecc")
        
        # Format plot
        ylabel("E [rad]")
        xlabel("M [rad]")
        xlim(0,2pi)
        ylim(0,2pi)
        legend(loc="upper left")
        grid()
    end
end

# end function

function testKeplerEquation2()
    # This function tests the ExoJulia keplerEquation solver's accuracy.
    # Kepler's Equations is solved for e ranging from [0, 0.9]
    # over all values of mean anomaly M.  The residuals M - (E - esinE) are plotted vs
    # M to probe the solver's accuracy.  Typically, we expect residuals ~ 1.0e-15
    
    # Inits
    len = 9
    len_M = 100
    e = linspace(0,0.9,len);
    M = linspace(0,2pi,len_M);
    E = zeros(M);
    residuals = zeros(E);
    ecc = 0.0
    
    # Loop over eccentricities
    for i=1:len
        # Loop over mean anomalies
        for j=1:len_M
            E[j] = keplerEquation(M[j], e[i])
        end

        # Compute, plot residuals vs M
        residuals = M .- (E .- e[i] .* sin(E))
        ecc = round(e[i],1)
        plot(M, residuals, lw=2, label="e = $ecc")
        
        # Format plot
        ylabel("M - (E - esinE)")
        xlabel("M [rad]")
        
        legend(bbox_to_anchor=[1.3, 1.02])
        grid()
    end
end