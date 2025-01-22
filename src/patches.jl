abstract type Patch end

struct OpenPatch <: Patch
    bound::Vector{Vector{Float64}}
end


function OpenPatch(patch::OpenPatch, branch::Real)
    bound = rebranch(patch.bound, branch)
    return OpenPatch(bound)
end

function OpenPatch(
    obs::AstroLib.Observatory,
    θs::Tuple{T, T},
    ϕs::Vector{T};
    ngeodesic=100
):: OpenPatch where T <: Real 
    
    pts = Vector{Float64}[]

    idx = 1
    θ = θs[1]
    for ϕ in ϕs
        push!(pts, zenazi_2_localcodecha(obs, θ, ϕ))
        idx += 1
    end

    θ2, ϕ2 = coδ, ha = zenazi_2_localcodecha(obs, θs[2], ϕs[end])
    geo_pts = geodesic_on_sphere(pts[idx-1][1], pts[idx-1][2], θ2, ϕ2)[2:end-1]

    for geo_pt in geo_pts
        push!(pts, geo_pt)
        idx +=1
    end
    
    θ = θs[2]
    for ϕ in reverse(ϕs)
        push!(pts, zenazi_2_localcodecha(obs, θ, ϕ))
        idx += 1
    end

    θ2, ϕ2 = coδ, ha = zenazi_2_localcodecha(obs, θs[1], ϕs[1])
    geo_pts = geodesic_on_sphere(pts[idx-1][1], pts[idx-1][2], θ2, ϕ2)[2:end-1]

    for geo_pt in geo_pts
        push!(pts, geo_pt)
        idx +=1
    end
    return OpenPatch(pts)
end

struct ClosedPatch <: Patch
    bound1::Vector{Vector{Float64}}
    bound2::Vector{Vector{Float64}}
end

function ClosedPatch(patch::ClosedPatch, branch::Real)
    bound1 = rebranch(patch.bound1, branch)
    bound2 = rebranch(patch.bound2, branch)
    return ClosedPatch(bound1, bound2)
end

function ClosedPatch(
    obs::AstroLib.Observatory,
    θs::Tuple{T, T}
) where T <: Real
    ϕs = LinRange(0, 2π-0.001, 200)
    bound1 = [zenazi_2_localcodecha(obs, θs[1], ϕ) for ϕ in ϕs]
    bound2 = [zenazi_2_localcodecha(obs, θs[2], ϕ) for ϕ in reverse(ϕs)]
    return ClosedPatch(bound1, bound2)
end

normalize_vector(v) = v / norm(v)

function geodesic_on_sphere(θ1, ϕ1, θ2, ϕ2; num_points=100)

    p1 = sph_to_cart(θ1, ϕ1)
    p2 = sph_to_cart(θ2, ϕ2)

    # Compute the axis of rotation and angle
    cross_prod = cross(p1, p2)
    dot_prod = dot(p1, p2)
    angle = acos(dot_prod)

    # Create a rotation matrix
    function rotation_matrix(axis, θ)
        x, y, z = axis
        c, s = cos(θ), sin(θ)
        R = [
            c + x^2*(1-c)      x*y*(1-c) - z*s  x*z*(1-c) + y*s
            y*x*(1-c) + z*s    c + y^2*(1-c)    y*z*(1-c) - x*s
            z*x*(1-c) - y*s    z*y*(1-c) + x*s  c + z^2*(1-c)
        ]
        return R
    end

    # Generate points along the geodesic
    geodesic_points = []
    for t in range(0, 1, length=num_points)
        θ = t * angle
        R = rotation_matrix(normalize_vector(cross_prod), θ)
        push!(geodesic_points, cart_to_sph((R * p1)...))
    end
    return geodesic_points
end
