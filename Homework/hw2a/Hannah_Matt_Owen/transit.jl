using PyPlot

"""
    circle_overlap(d, r_s, r_p)
Computer the area of overlap of two circles as a function of the distance
between them.

#Arguments
* `d`: the straight-line distance between the centers of the two bodies
* `r_s`: the radius of the star
* `r_p`: the radius of the planet, assumed that r_p < r_s
"""
function circle_overlap(d,r_s,r_p)
  if d <= r_s - r_p
    #the planet is completely within the star
    area = pi*r_p^2
  elseif d >= r_s + r_p
    #there is no overlap
    area = 0
  else
    phi_s = 2*(acos((r_s^2+d^2-r_p^2)/(2*r_s*d)))
    phi_p = 2*(acos((r_p^2+d^2-r_s^2)/(2*r_p*d)))
    area = 0.5*(phi_p*r_p^2-r_p^2*sin(phi_p)+phi_s*r_s^2-r_s^2*sin(phi_s))
  end

  return area
end


function read_data_file()
  data = readdlm("mystery_planet2.txt")
  time = data[:,1]
  flux = data[:,2]
  err = data[:,3]
  return (time, flux, err)
end

time, flux, err = read_data_file()
plot(time,flux)
show()
