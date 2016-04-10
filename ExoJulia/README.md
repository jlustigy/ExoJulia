# `ExoJulia`

To use `ExoJulia`: 
```julia
# Add this directory to your path
# path = "/PATH/TO/ExoJulia/ExoJulia"
path = "." # ...if in this directory
push!(LOAD_PATH, path)

# import ExoJulia
using ExoJulia
```
### Kepler's Equation
```julia
mean_anomaly = 0.5
eccentricity = 0.1
eccentric_anomaly = ExoJulia.Orbit.kepler_solve!(mean_anomaly, eccentricity)
```
