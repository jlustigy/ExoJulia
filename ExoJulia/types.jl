######## RV ########

type RadialVelocityData
    time::Array{Float64, 1}
    rv::Array{Float64, 1}
    error::Array{Float64, 1}
end

function RadialVelocityData(time::Array{Float64, 1}, rv::Array{Float64, 1}, error::Array{Float64, 1})
    RadialVelocityData(time, rv, error)
end

######## Transit ########

type TransitData
    time::Array{Float64, 1}
    flux::Array{Float64, 1}
    error::Array{Float64, 1}
end

function TransitData(time::Array{Float64, 1}, flux::Array{Float64, 1}, error::Array{Float64, 1})
    TransitData(time, flux, error)
end
