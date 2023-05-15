using StaticArrays
using LinearAlgebra

function sample_random(N)

    α = π * (2 * rand(N) .- 1) # uniform [-pi,pi]
    β = acos.(2 * rand(N) .- 1) # cos-uniform [-1,1]

    positions = @. SVector(cos(α) * sin(β), sin(α) * sin(β), cos(β))
    return positions

end


# spherical interpolation
# https://en.wikipedia.org/wiki/Slerp
function slerp(p0, p1, d)
    R = norm(p0)
    dt = d / R # angle corresponding to desired spacing
    p0n = p0 / R
    p1n = p1 / norm(p1)
    Omega = acos(dot(p0n, p1n))
    t = dt / Omega
    p = sin((1 - t) * Omega) / sin(Omega) * p0 + sin(t * Omega) / sin(Omega) * p1
    return p
end


function sample_landings(N, threshold_d, dimer_d; maxiter=1e3, k=0)

    if sqrt(4 * pi * 1^2 / N) < (pi * threshold_d^2)
        @warn "The requested number of points will not fit"
    end


    #initial sample
    s = sample_random(N + k)
    sold = s

    indices = trues(N + k) # all assumed good
    dimers = .!indices

    # first pass, checking distances
    for i in 1:(N+k) # points to test
        for j in (i+1):(N+k)
            dist = norm(s[i] - s[j])
            if (dist < threshold_d)  # this i point is bad
                indices[i] = false
                break # bad point, no need to test further
            end
        end
    end

    todo = (sum(indices) < N)

    # if more than N, we're done without dimers, return first N positions 
    if !todo
        @info "no iteration needed, zero dimers"
        return (s=s[1:N], dimers=dimers[1:N])
    end

    # otherwise, move pairs too close to set dimer distance and try again
    iter = 0
    while todo

        number_bad = sum(.!indices)

        @info "iteration $iter, $number_bad bad"
        for i in 1:(N+k) # points to test
            for j in (i+1):(N+k)
                dist = norm(s[i] - s[j])
                if (dist < threshold_d)  # this pair of points is too close, make it a dimer
                    dimers[i] = true
                    dimers[j] = true
                    # now assume this pair is fine until proven otherwise below
                    indices[i] = true
                    indices[i] = true
                    # shift j along the great circle to fixed separation d
                    newp = slerp(s[i], s[j], dimer_d)
                    s[j] = newp
                end
            end
        end

        # redo pass, checking distances
        for i in 1:(N+k) # points to test
            for j in (i+1):(N+k)
                dist = norm(s[i] - s[j])
                if (dist < threshold_d)  # this i point is bad
                    indices[i] = false
                    break # bad point, no need to test further
                end
            end
        end

        number_bad = sum(.!indices)
        @info "now $number_bad bad"
        # if more than N, we're done
        if sum(indices) >= N
            @info "iterations successful"
            pick = findall(indices)
            return (s=s[pick[1:N]], dimers=dimers[pick[1:N]])
        end

        if iter >= maxiter
            @warn "max number of iterations reached"
        end

        todo = (sum(indices) < N) && (iter < maxiter)
        iter = iter + 1
    end

    # we should not get here ideally
    return (s=s[1:N], dimers=dimers[1:N])

end

R = 14
threshold_d = 0.3
dimer_d = 0.5
sample_landings(3000, threshold_d / R, dimer_d / R, k=30)
