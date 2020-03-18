# ---------------------------------
# PIECEWISE LINEAR UPPER ENVELOPE |
# ____________________ __________ |
# Code written by Gustavo Pereira |
#   Columbia University, 2020     |
# ---------------------------------


mutable struct PiecewiseLinear{T}
    xcoords::Vector{T}
    ycoords::Vector{T}
end

# These methods contribute very little, except for readability later on.
@inline Base.eltype(seg::PiecewiseLinear{T}) where T=T
@inline Base.length(seg::PiecewiseLinear)=length(seg.xcoords)
@inline Base.first(seg::PiecewiseLinear)=first(seg.xcoords)
@inline Base.last(seg::PiecewiseLinear)=last(seg.xcoords)
@inline Base.zero(seg::PiecewiseLinear{T}) where T=zero(T)
@inline @inbounds get_x(seg::PiecewiseLinear, i)=seg.xcoords[i]
@inline @inbounds get_y(seg::PiecewiseLinear, i)=seg.ycoords[i]
@inbounds function set_knot!(seg, i, x)
    seg.xcoords[i]=x
    nothing
end
@inline @inbounds get_y(seg::PiecewiseLinear, i)=seg.ycoords[i]
@inbounds function set_value!(seg, y, i)
    seg.values[i]=y
    nothing
end
@inline get_point(seg, i)=Point2(get_x(seg, i), get_y(seg, i))

# Create a Point2 type for the sole purpose of overloading base arithmetics for
# 2-Tuples without actually doing so :-)
struct Point2{T}
    x::T
    y::T
end
@inline get_x(pt::Point2)=pt.x
@inline get_y(pt::Point2)=pt.y
@inline function set_point!(seg, pt, i)
    @inbounds seg.xcoords[i]=pt.x
    @inbounds seg.ycoords[i]=pt.y
end


Base.Tuple(u::Point2)=(u.x, u.y)
Base.:+(u::Point2, v::Point2)=Point2(u.x+v.x, u.y+v.y)
Base.:*(u::Point2, k)=Point2(k*u.x, k*u.y)
Base.:*(k, u::Point2)=Base.:*(u, k)
Base.:-(u::Point2, v::Point2)=Point2(u.x - v.x, u.y - v.y)
rot90cc(u::Point2)=Point2(-u.y, u.x) # rotate 90 deg counter clockwise
pdot(u::Point2, v::Point2)=u.x*v.x+u.y*v.y
inner_rotate(v::Point2, w::Point2)=pdot(v, rot90cc(w))



# * Main function

function compute_envelope(segments)
    envsize_max = length(segments[1]) + length(segments[2]) +
        max(length(segments[1]), length(segments[2]))
    envelope = PiecewiseLinear(map(k -> fill(NaN, envsize_max), (1,2))...)


    @inbounds seg_state = determine_first(segments)
    v_idx_visit = [1, 1]

    i_env = 1


    while is_lower_x(segments, seg_state, v_idx_visit)
        if v_idx_visit[seg_state] < length(segments[seg_state])

            this_segment = segments[seg_state]
            this_ind = v_idx_visit[seg_state]
            pt = get_point(this_segment, this_ind)
            set_point!(envelope, pt, i_env); i_env += 1

            # Increase counter for current segment & insert points in envelope

            inc_counter!(v_idx_visit, seg_state)
        else
            break
        end
    end

    while get_x(segments[1], v_idx_visit[1]) == get_x(segments[2], v_idx_visit[2])
        if get_y(segments[1], 1) > get_y(segments[2], 1)
            set_point!(envelope, get_point(segments[1], v_idx_visit[1]), i_env); i_env += 1
        else
            set_point!(envelope, get_point(segments[2], v_idx_visit[2]), i_env); i_env += 1
        end
        inc_counter!(v_idx_visit)
    end


    @inbounds seg_state = determine_first(segments, v_idx_visit)


    @inbounds while (v_idx_visit[1] <= length(segments[1])) && (v_idx_visit[2] <= length(segments[2]))
        i_curr = v_idx_visit[seg_state]
        i_other = v_idx_visit[other(seg_state)]

        pt_curr = get_point(segments[seg_state], i_curr)
        pt_other = get_point(segments[other(seg_state)], i_other)

        pt_other_prev = get_point(segments[other(seg_state)], i_other - 1) # TODO: This could be saved from a previous interaction
        sco = inner_rotate(pt_curr - pt_other_prev, pt_other - pt_other_prev)

        if sco > 0
            set_point!(envelope, pt_curr, i_env); i_env += 1
        end

        if i_curr + 1 > length(segments[seg_state])
            seg_state = other(seg_state)
            break
        end

        pt_curr_next = get_point(segments[seg_state], i_curr+1)

        i_env += segment_intersect!(envelope, i_env, sco, pt_other_prev, pt_other, pt_curr, pt_curr_next)
        inc_counter!(v_idx_visit, seg_state)

        if get_x(pt_curr_next) >= get_x(pt_other)
            seg_state = other(seg_state)
        end

    end


    # Now finalize by adding all leftmost nodes
    while v_idx_visit[seg_state] <= length(segments[seg_state])
        set_point!(envelope, get_point(segments[seg_state], v_idx_visit[seg_state]), i_env); i_env = i_env + 1
        inc_counter!(v_idx_visit, seg_state)
    end

    return envelope
end





# * Auxiliary functions
other(i)=3-i
inc_counter!(v_idx, i)=v_idx[i] += 1
inc_counter!(v_idx)=begin
    v_idx[1] += 1; v_idx[2] += 1;
end

function is_lower_x(segments, seg_state, v_idx_visit)
    seg_other = other(seg_state)

    @inbounds i_curr = v_idx_visit[seg_state]
    @inbounds i_other = v_idx_visit[seg_other]
    get_x(segments[seg_state], i_curr) < get_x(segments[seg_other], i_other)
end



function append_unsafe!(pl1, pl2)
    append!(pl1.xcoords, pl2.xcoords)
    append!(pl1.ycoords, pl2.ycoords)
end


Base.@propagate_inbounds function determine_first(segments, init_ind=(1, 1))
    r = 1
    if get_x(segments[2], init_ind[2]) < get_x(segments[1], init_ind[1])
        r = 2
    end
    return r
end


function segment_intersect!(v_store, pos_store, sco, ptA1, ptA2, ptB1, ptB2)
    # Returns boolean for whether intersection exists
    #
    # Note that ptA1, ptA2, ptB1 and ptB2 are NOT general points in the plane.
    # Because of the problem's structure, the following inequality holds:
    # coord_x(ptA1) <= coord_x(pt_B1)
    # where `coord_x` stands for the first coordinate...

    # Note that `sco` already stores a value whose sign depends on the half
    # space of point ptB1 relative to the line formed by the segment ptA1-ptA2
    # Hence, the first way of ruling out the intersection is to compute the
    # same value for ptB2;


    # Formally, `sco` is the inner product of `R` and the 90 deg
    # counter-clockwise rotation of `P` as defined below.

    P = ptA2 - ptA1
    Q = ptB2 - ptB1
    R = ptB2 - ptA1


    sco_cmp = inner_rotate(R, P)

    if sign(sco_cmp) == sign(sco)
        return false
    end

    # If sign switches, this means there is _possibly_ an intersection.
    #
    # Use the geometry of the problem to find the coefficient:
    bb = sco_cmp/(sco_cmp - sco)
    intersection_candidate = bb * ptB1 + (1-bb) * ptB2

    r = get_x(intersection_candidate) < get_x(ptA2)

    r && set_point!(v_store, intersection_candidate, pos_store)

    return r
end
