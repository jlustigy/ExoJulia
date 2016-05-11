
using PyPlot
using LsqFit


data = readdlm("mystery_planet2.txt") ;

time = data[:,1]
flux = data[:,2]
err = data[:,3] ;

flux

mean(flux)

scatter(time,flux/mean(flux))
#xlim(0,10)
#ylim(0.994,1.001)

function Period(time, rv)
    num = 3000
    periods = linspace(0.0, 15.0, num)
    chisq = ones(num)
    for i in range(1, num)
        mods = mod(time, periods[i])
        ind = sortperm(mods)
        newrv = rv[ind]
        sum = 0.0
        for k in range(1,length(time)-1)
            sum += (newrv[k+1]-newrv[k])^2.
        end
        chisq[i] = sum
    end
    #scatter(periods,chisq)
    #yscale("log")
    #xlim(0,100)
    return periods[indmin(chisq)] 
end

function CircleOverlap(d, r1, r2)
    r = min(r1, r2)
    R = max(r1, r2)
    if d >= (r + R)
        return 0.0
    end
    if d <= (R-r)
        return pi*r^2
    end

    root_term = 0.5*sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R))
    r_term = r^2*acos((d^2 + r^2 - R^2)/(2.0*d*r))
    R_term = R^2*acos((d^2 + R^2 - r^2)/(2.0*d*R))
    return r_term + R_term - root_term
end


#p[b, P, rho, k, t_off]
#p[1] = b impact parameter
#p[2] = P period
#p[3] = rho stellar density
#p[4] = k r_planet/r_star
#p[5] = t_0 midpoint of first transit

function transit(t, p)
    G = 6.67259E-8 #cm^3 g^-1 cm^-2
    b = p[1]
    P = p[2]
    rho = p[3]
    k = p[4]
    t_0 = p[5]
    f0 = p[6]
    if rho <= 0.0 
        return Inf
        elseif rho > 5
        return Inf
    end
    v = (8.0 * pi^2.0 * G * rho / (3.0 * P*24.0*3600.0))^(1.0/3.0) *3600.*24.
    #T = (3.0 / pi^2.0)^(1.0/3.0) * sqrt(1.0-b^2.0) * (P/(2.0*pi))^(1.0/3.0) * (G*rho)^(-1.0/3.0)
    
    #Define array for the distance between the center of the planet and center of star at each time in t
    flux = ones(length(t))
    
    #Check for testing parameters ouside of allowed values which are:
    # 0 =< b < 1
    # 0.01 < k < 0.07
    # 0 < rho < 5 just to give some limit on the bounds of stellar density should be ~ 1
    if abs(b) > 1.0 
        return Inf
    end
    
    for i in range(1, length(t))
        del_t = mod((t[i]-t_0+P/2),P)-P/2
        d = sqrt((del_t*v)^2+b^2)
        flux[i] = 1.0 - CircleOverlap(d, k, 1.0)/pi
    end
    return flux*f0
end

index = zeros(length(time))
for i in range(1,length(time))
    if time[i] < 2.0
        index[i] = i
    else
        index[i] = 0
    end
end

transit_bottom_index = maximum(index)
j = round(Int,transit_bottom_index)
t_off_guess = 0.0
for i in range(1, j)
    if flux[i] == minimum(flux[1:j])
        t_off_guess = time[i]
    end
end
per_guess = Period(time,flux)
println("Guess for timing offset = ", t_off_guess)
println("Guess for period = ", per_guess)

#test = [0.0, per_guess, 1.0, 0.05, t_off_guess, f0]
#F = transit(time, test);

errors = 1.0 ./ (err ./ mean(flux)).^2.0

#p[b, P, rho, k, t_off, f0]
test = [0.5, per_guess, 0.2, 0.041, 0.85, 1.0]
fit = curve_fit(transit, time, flux/mean(flux), errors, test)
tfit = transit(time, fit.param);

#Fit parameters 
println("Impact parameter b = ", fit.param[1])
println("Period P = ", fit.param[2])
println("Stellar density rho = ", fit.param[3])
println("Relative planet radius k = ", fit.param[4])
println("Time until first transit = ", fit.param[5])

scatter(time, flux/mean(flux), color = "b", label = "data", s=0.5)
plot(time, tfit, color = "r", label = "fit")
ylim(0.995, 1.002)
legend()

