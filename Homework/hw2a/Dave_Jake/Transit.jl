module Transit

# Import necessary modules
using LsqFit
using Optim

# Define global variables for the data
global time = nothing
global flux = nothing
global err = nothing

############# Create Setter functions #############

function set_time(val)
   global time
   time = val
   println("Set Transit.time")
end

function set_flux(val)
   global flux
   flux = val
   println("Set Transit.flux")
end

function set_err(val)
   global err
   err = val
   println("Set Transit.err")
end

############# Include Fitting functions #############


end
