using Base.Iterators: cycle

struct Crossing
    pt0::Vector{Float64}
    pt1::Vector{Float64}
    x0::Float64
end

function Crossing(pt0, pt1)
    m = (pt0[2] - pt1[2]) / (pt0[1] - pt1[1])
    x0 = -pt0[2] / m  + pt0[1]
    return Crossing(pt0, pt1, x0)
end

function find_next_crossing(crossing, crossings)
    current_idx = findfirst([c==crossing for c in crossings])
    R = RotMatrix2(-π / 2)
    tangent = crossing.pt0 - crossing.pt1
    s = Int(sign((R * tangent)[1]))
    if s + current_idx==0
        return "origin"
    elseif s + current_idx > length(crossings)
        return "infinity"
    else
        return crossings[s + current_idx]
    end
end

function crossed_branch(pt0, pt1)
    d = (sign(pt0[2]) - sign(pt1[2]))==0
    if d
        return false
    end
    m = (pt1[2] - pt0[2]) / (pt1[1] - pt0[1])
    #@show pt0, pt1
    #@show m
    x0 = (m * pt0[1] - pt0[2]) / m
    #@show x0
    return x0 >= 0
end

function split_along_branch(patch::ClosedPatch, branch::Real; small=1e-5, big=1e5, n=50)
    holes = []
    solids = []
    crossings = []
    R = RotMatrix2(-branch)
    for (idx, bound) in enumerate([patch.bound1, patch.bound2])
        stereo_bound = map(x->R * stereographic(x[1], x[2]), bound)
        stereo_bound = filter(x->x[2]!=0.0, stereo_bound)

        m = map(
            xy->crossed_branch(xy[1], xy[2]),
            zip(stereo_bound, circshift(stereo_bound, -1))
        )

        if sum(m)==0
            current_patch = map(x->R' * x, stereo_bound)
            push!(holes, map(x->inv_stereographic(x[1], x[2]; branch=branch), current_patch))
            continue
        end

        for jdx in findall(m)
            x0 = stereo_bound[jdx]
            x1 = jdx < length(m) ? stereo_bound[jdx+1] : stereo_bound[1]
            push!(crossings, (Crossing(x0, x1), idx))
        end
        push!(solids, stereo_bound)
    end

    if length(solids)==0
        return [], holes
    end

    sort!(crossings, by=x->x[1].x0)
    out, current_patch = [], []

    itr, state = cycle(solids[1]), 1
    while sum(map(x->length(x), solids)) > 0
        pt, state = iterate(itr, state)
        if R' * pt in current_patch
            push!(out, map(x->inv_stereographic(x[1], x[2]; branch=branch), current_patch))
            for solid in solids
                filter!(x->(R' * x) ∉ current_patch, solid)
                if length(solid)>1
                    itr, state = cycle(solid), 1
                end
            end
            current_patch = []
            continue
        end

        push!(current_patch, R' * pt)

        if pt in map(x->x[1].pt0, crossings)
            current_crossing, _ = crossings[findfirst(map(x->all(x[1].pt0.==pt), crossings))]
            next_crossing = find_next_crossing(current_crossing, map(x->x[1], crossings))
            sgn = sign(current_crossing.pt0[2])
            if next_crossing=="infinity"
                crossing = current_crossing
                pts1 = [R' * [x, sgn * small] for x in LinRange(crossing.x0, big, n)]
                pts2 = [R' * [big*cos(ψ), big*sin(ψ)] for ψ in LinRange(small, 2π - small, n)]
                if sgn==1
                    reverse!(pts2)
                end
                pts3 = [R' * [x, -sgn * small] for x in LinRange(big, crossing.x0, n)]
                current_patch = vcat(current_patch, pts1, pts2, pts3)
            elseif next_crossing=="origin"
                crossing = current_crossing
                pts1 = [R' * [x, sgn * small] for x in LinRange(crossing.x0, -small, n)]
                pts2 = [R' * [x, -sgn * small] for x in LinRange(-small, crossing.x0, n)]
                x = map(x->inv_stereographic(x[1], x[2]), pts1)
                y = map(x->inv_stereographic(x[1], x[2]), pts2)
                current_patch = vcat(current_patch, pts1, pts2)
            else
                pt0 = [current_crossing.x0, sgn * small]
                pt1 = [next_crossing.x0, sgn * small]
                itp = linear_interpolation([0, 1], [pt0, pt1])
                pts = [R' * pt for pt in itp.(LinRange(0, 1, n))]
                current_patch = vcat(current_patch, pts)
                next_crossing, idx = crossings[findfirst(map(x->x[1]==next_crossing, crossings))]
                solid = solids[idx]
                #@show solid
                if length(solid)==0
                    continue
                end
                idx = findfirst([all(isapprox.(x, next_crossing.pt1)) for x in solid])
                #@show pt
                #@show next_crossing
                #@show idx
                circshift!(solid, -idx+1)
                itr, state = cycle(solid), 1
            end
        end
    end
    return out, holes
end

function split_along_branch(patch::OpenPatch, branch::Real; small=1e-5, big=1e5, n=50)

    R = RotMatrix2(-branch)
    solid = map(x->R * stereographic(x[1], x[2]), patch.bound)
    solid = filter(x->x[2]!=0.0, solid)
    m = map(
        xy->crossed_branch(xy[1], xy[2]),
        zip(solid, circshift(solid, -1))
    )

    crossings = []
    for idx in findall(m)
        x0 = solid[idx]
        x1 = idx!=length(m) ? solid[idx+1] : solid[1]
        push!(crossings, Crossing(x0, x1))
    end
    if length(crossings)==0
        return [patch.bound]
    end

    sort!(crossings, by=x->x.x0)

    out, current_patch = [], []
        
    itr, state = cycle(solid), 1
    pt, _  = iterate(itr, state)
    pt, _  = iterate(itr, 2)
    while length(solid) > 0
        pt, state = iterate(itr, state)
        if R' * pt in current_patch
            push!(out, map(x->inv_stereographic(x[1], x[2]; branch=branch), current_patch))
            filter!(x->(R' * x) ∉ current_patch, solid)
            if length(solid) > 1
                itr, state = cycle(solid), 1
            end
            current_patch = []
            continue
        end
            
        push!(current_patch, R' * pt)
            
        if pt in map(x->x.pt0, crossings)
            current_crossing = crossings[findfirst(map(x->all(x.pt0.==pt), crossings))]
            next_crossing = find_next_crossing(current_crossing, map(x->x, crossings))
            sgn = sign(current_crossing.pt0[2])
            if next_crossing=="infinity"
                crossing = current_crossing
                pts1 = [R' * [x, sgn * small] for x in LinRange(crossing.x0, big, n)]
                pts2 = [R' * [big*cos(ψ), big*sin(ψ)] for ψ in LinRange(small, 2π - small, n)]
                if sgn==1
                    reverse!(pts2)
                end         
                pts3 = [R' * [x, -sgn * small] for x in LinRange(big, crossing.x0, n)]
                current_patch = vcat(current_patch, pts1, pts2, pts3)
            elseif next_crossing=="origin"
                crossing = current_crossing
                pts1 = [R' * [x, sgn * small] for x in LinRange(crossing.x0, -small, n)]
                pts2 = [R' * [x, -sgn * small] for x in LinRange(-small, crossing.x0, n)]
                current_patch = vcat(current_patch, pts1, pts2)
            else
                pt0 = [current_crossing.x0, sgn * small]
                pt1 = [next_crossing.x0, sgn * small]
                itp = linear_interpolation([0, 1], [pt0, pt1])
                pts = [R' * pt for pt in itp.(LinRange(0, 1, n))]
                current_patch = vcat(current_patch, pts)
                idx = findfirst([all(x.==next_crossing.pt1) for x in solid])
                circshift!(solid, -idx+1)
                itr, state = cycle(solid), 1
            end
        end
    end
    return out                        
end
