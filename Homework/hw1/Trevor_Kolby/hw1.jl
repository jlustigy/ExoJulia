#.jl version of hw1.ipynb
#@stest newtons(g,g_prime,E0(1.23456),0.5,1.23456)


# From time to M
M(t::Float64,t_p::Float64,P::Float64) = (2*pi/P)*(t-t_p)

# this should evaluate to zero when E is correct
g(E::Float64,e,m) = E - e*sin(E) - m

# g'(E)
g_prime(E::Float64,e) = 1.0 - e*cos(E)

#+1 for positive x, -1 for negative, with proper handling of zero division
function sign(x::Float64)
    if x!= 0.0
        return x/abs(x)
    else
        return 0.0
    end
end

#Best guess at initial E0
E0(m::Float64) = m + 0.85*sign(sin(m))

#Newtons takes g, g', a guess at E0 to make the g(E0)=0 and the parameters to the function. delta can be specified depending on your machine's precision. Output is the true E
function newtons(f::Function,f_prime::Function,E_0::Float64,e::Float64,m::Float64;delta=1e-14)
    
    E_new = E_0
    dE = 1.
    
    while abs(dE)>delta
        E_old = E_new
        E_new = E_old - (f(E_old,e,m)/f_prime(E_old,e))
        dE = E_new-E_old
    end
    
    return E_new
    
end

#How to convert from E to f
EtoF(E::Float64,e::Float64) = 2.*atan((((1.+e)/(1.-e))^(1./2.))*tan(E/2.))

#How to convert from E to r
EtoR(E::Float64,e::Float64,a::Float64) = a*(1-e*cos(E))

function orbit(e::Float64,sma::Float64,period::Float64;numpoints=1000)
    
    #returns a tuple of two arrays, the radii and f values of the whole orbit (0<=M<2pi)
    #This will let us check to see that the orbits look right
    
    #e is eccentricity
    
    #sma is the smi major axis
    
    #period is what it sounds like\
    
    #numpoints is an optional keyword parameter, which tells it how many points to compute
    #along the orbit. While all the points computed are correct, it looks better when you plot
    #with lots of points.
    
    rs = []
    fs = []
    
    times = linspace(0,period+(1.0/numpoints),numpoints)
    
    for time in times
        em = M(time,0.,period)
        E = newtons(g,g_prime,E0(em),e,em)
        push!(fs,EtoF(E,e))
        push!(rs,EtoR(E,e,sma))
    end
    
    return (rs,fs)
    
end

