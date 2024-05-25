using AstroLib

"""
Converts a θ and ϕ from TAMBO coordinates to the coordinate system used by AstroLib.jl
This has the ẑ vector pointing perpendicular to Earth's surface and the 
azimuth 0 lying along the North direction with positive angles resulting
from left-handed rotations

# Examples
```jldoctest
julia> to_astrolib(π / 4, π / 2)
(0.7853981633974483, 3.492899615108751)
```
"""
function to_astrolib(θ::Real, ϕ::Real)
    ϕ = π / 2 - ϕ
    ϕ = ϕ < 0 ? 2π + ϕ : ϕ
    return θ, ϕ
end

"""
Converts a TAMBO direction to the coordinate system used by AstroLib.jl
This has the ẑ vector pointing perpendicular to Earth's surface and the 
azimuth 0 lying along the North direction with positive angles resulting
from left-handed rotations

# Examples
```jldoctest
julia> using Tambo

julia> d = Tambo.Direction(rand() * π, rand() * 2 * π)

julia> to_astrolib(d);
```
"""
function to_astrolib(d::Tambo.Direction)
    to_astrolib(d.θ, d.ϕ)
end

"""
Convert zenith to altitude

# Examples
```jldoctest
julia> zenith_to_latitude(π / 4)
0.7853981633974483
"""
function zenith_to_latitude(θ::Real)
    return π / 2 - θ
end 


function latitude_to_zenith(lat::Real)
    return π / 2 - lat
end

function tambo_to_radec(θ::Real, ϕ::Real, jd::Real, lat::Real, long::Real)
    θ, ϕ = to_astrolib(θ, ϕ)
    latitude = zenith_to_latitude(θ)
    ra, δ, _ = hor2eq.(rad2deg(latitude), rad2deg(ϕ), jd, rad2deg(lat), rad2deg(long))
    return deg2rad(ra), deg2rad(δ)
end
