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
    ra = π / 2 - ϕ
    ra = ϕ < 0 ? 2π + ϕ : ϕ
    return θ, ra
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


function zenazi_2_localcodecha(obs, θ, ϕ)
    ha, δ = deg2rad.(altaz2hadec(rad2deg(π / 2 - θ), rad2deg(ϕ), rad2deg(obs.latitude)))
    coδ = π / 2 - δ
    ha += obs.longitude
    # Make sure it is in range 0 <= ha < 360
    ha = mod(ha + 2π, 2π)
    # Move hour angle to -180 < ha <= 180
    ha = π < ha ? ha - 2π : ha
    return [coδ, ha]
end

function stereographic(θ, ϕ; rmin=1e-5, rmax=1e5)
    r = 2 * tan((π - θ) / 2)
    r = minimum([r, rmax])
    r = maximum([r, rmin])
    x, y = r * cos(ϕ), r * sin(ϕ)
    return [x, y]
end

function inv_stereographic(x, y; rmin=1e-5, rmax=1e5, branch=0.0)
    ϕ = rebranch(atan(y, x), branch)
    r = sqrt(x^2 + y^2)
    θ = π - 2 * atan(r / 2)
    return [θ, ϕ]
end

function sph_to_cart(θ, ϕ)
    return [cos(ϕ) * sin(θ), sin(ϕ) * sin(θ), cos(θ)]
end

function cart_to_sph(x, y, z)
    z /= sqrt(x^2 + y^2 + z^2)
    return [acos(z), atan(y, x)]
end

function rebranch(v::Real, branch::Real)
    while v < branch
        v += 2π
    end
    while branch + 2π <= v
        v -= 2π
    end
    return v
end

function rebranch(pts::Vector{Vector{Float64}}, branch::Real)
    for pt in pts
        pt[2] = rebranch(pt[2], branch)
    end
    return pts
end
