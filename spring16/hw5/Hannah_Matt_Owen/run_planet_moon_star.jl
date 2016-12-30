include("planet_moon_star.jl")

yplanet = 0.
y1moon = 0
y2moon = 0
eps1 = 1e-5
eps2 = 0.5
map = planet_moon_star(yplanet,y1moon,y2moon,eps1,eps2)
print(size(map))

#print(x)
radii = collect(0.25:0.7:4)
max = zeros(length(radii))

for (k,radius) in enumerate(radii)
    x = zeros(Float64, (1000,1000))
    for (i=1:1000)
        for(j=1:1000)
            if ((i-500)^2+(j-500)^2 < (radius*100)^2)
                x[i,j] = 1
            end
        end
    end
    println(k)
    max[k] = maximum(conv2(map,x))/(pi*radius^2)
end


using PyPlot
figure(figsize=(12,8))
plot(radii,max)
#show(radii,max,extent=[-4,4,-1,1])
show()