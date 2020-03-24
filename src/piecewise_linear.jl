# ---------------------------------
# PIECEWISE LINEAR UPPER ENVELOPE |
# ____________________ __________ |
#    Author: Gustavo Pereira      |
#   Columbia University, 2020     |
# ---------------------------------


abstract type AbstractPiecewiseLinear{T} end
@inline Base.length(fun::AbstractPiecewiseLinear)=length(get_x(fun))
@inline Base.eltype(fun::AbstractPiecewiseLinear{T}) where T=T
@inline Base.zero(fun::AbstractPiecewiseLinear{T}) where T=zero(T)
Base.Tuple(fun::AbstractPiecewiseLinear)=(get_x(fun), get_y(fun))
Base.resize!(fun::AbstractPiecewiseLinear, n)=begin
    resize!(get_x(fun), n);resize!(get_y(fun), n);
end

"""
# Piecewise linear function

    mutable struct PiecewiseLinear{T}
      xcoords::Vector{T}
      ycoords::Vector{T}
    end

The purpose is to represent "piecewise linear" functions by storing the x and y
coordinates of the points that define them. The i-th point has x coordinate
xcoords[i], and y coordinate `ycoords[i]`

Convenience methods `Base.getindex` and `Base.setindex!` are defined for this
type. Accessing index `i` returns x and y coordinates grouped in a
UpperEnvelope.Point2 type. That type should also be used for insertion.


## Assumptions
**The field `xcoords` should be ordered.** Moreover, `xcoords` and `ycoords`
must be of the same length.


"""
mutable struct PiecewiseLinear{T} <: AbstractPiecewiseLinear{T}
    xcoords::Vector{T}
    ycoords::Vector{T}
    function PiecewiseLinear(xcoords, ycoords)
        length(xcoords) != length(ycoords) && error("x and y coordinates must have the same length")
        new{eltype(xcoords)}(xcoords, ycoords)
    end
end

# These methods contribute very little, except for readability later on.
@inline get_x(fun::PiecewiseLinear, i)=fun.xcoords[i]
@inline get_x(fun::PiecewiseLinear)=fun.xcoords
function set_x!(fun::PiecewiseLinear, i, x)
    fun.xcoords[i]=x
    nothing
end
@inline get_y(fun::PiecewiseLinear, i)=fun.ycoords[i]
@inline get_y(fun::PiecewiseLinear)=fun.ycoords
function set_y!(fun::PiecewiseLinear, y, i)
    fun.values[i]=y
    nothing
end
Base.@propagate_inbounds @inline Base.getindex(fun::PiecewiseLinear, i)=Point2(get_x(fun, i), get_y(fun, i))


mutable struct ExtendedPiecewiseLinear{T, N} <: AbstractPiecewiseLinear{T}
    fun_pl::PiecewiseLinear{T}
    subordinate_vectors::NTuple{N, Vector{T}}
end
@inline get_x(fun::ExtendedPiecewiseLinear)=get_x(fun.fun_pl)
@inline get_x(fun::ExtendedPiecewiseLinear, i)=get_x(fun.fun_pl, i)
@inline set_x!(fun::ExtendedPiecewiseLinear, i, x)=set_x!(fun.fun_pl, i, x)

@inline get_y(fun::ExtendedPiecewiseLinear, i)=get_y(fun.fun_pl, i)
@inline set_y!(fun::ExtendedPiecewiseLinear, i, y)=set_y!(fun.fun_pl, i, y)
Base.@propagate_inbounds @inline Base.getindex(fun::ExtendedPiecewiseLinear, i)=getindex(fun.fun_pl, i)


"""

Created with sole purpose of overloading base arithmetics for 2-Tuples without
messing with Tuples.

The benefit is simply notational convenience.

"""
struct Point2{T}
    x::T
    y::T
end
@inline get_x(pt::Point2)=pt.x
@inline get_y(pt::Point2)=pt.y
@inline Base.setindex!(fun::PiecewiseLinear, pt, i)=begin
    fun.xcoords[i]=pt.x
    fun.ycoords[i]=pt.y
end


Base.Tuple(u::Point2)=(u.x, u.y)
Base.:+(u::Point2, v::Point2)=Point2(u.x+v.x, u.y+v.y)
Base.:*(u::Point2, k)=Point2(k*u.x, k*u.y)
Base.:*(k, u::Point2)=Base.:*(u, k)
Base.:-(u::Point2, v::Point2)=Point2(u.x - v.x, u.y - v.y)

" 90deg counter clockwise rotation "
rot90cc(u::Point2)=Point2(-u.y, u.x) # rotate 90 deg counter clockwise

" Dot product of `Point2` types."
pdot(u::Point2, v::Point2)=u.x*v.x+u.y*v.y

" Rotate the second argument and take dot product with the first one. "
inner_rotate(v::Point2, w::Point2)=pdot(v, rot90cc(w))

# * Main function

"""
# Compute upper envelope of piecewise linear functions

    compute_envelope(pl_funcs::NTuple{2, PiecewiseLinear})

Computes the upper envelope of two piecewise linear functions, as represented by
type PiecewiseLinear.

## Algorithm

Modulo dealing with the non-intersecting part of domains, the idea is to fix a
pair of points (i.e., a segment) in one of the functions -- call it "benchmark"
-- and analyse all segments of the other function -- call it "comparison" --
within that segment. It is easy to check whether the subsegments are "dominated"
or not, and to compute intersections.

When all comparison segments are analysed, switch the states between functions:
the benchmark function becomes the comparison one and vice-versa. If this step
is ignored, intersections **will be** overlooked.
"""
function compute_envelope(pl_funcs::NTuple{2, PiecewiseLinear}, remove_nan=true)
    envsize_max = length(pl_funcs[1]) + length(pl_funcs[2]) +
        max(length(pl_funcs[1]), length(pl_funcs[2]))
    envelope = PiecewiseLinear(map(k -> fill(NaN, envsize_max), (1,2))...)

    @inbounds seg_state = determine_first(pl_funcs)
    v_idx_visit = [1, 1]
    i_env = 1

    # Deal the left part of (possibly) non-intersecting domains
    while is_lower_x(pl_funcs, seg_state, v_idx_visit)
        if v_idx_visit[seg_state] < length(pl_funcs[seg_state])
            this_ind = v_idx_visit[seg_state]
            pt = pl_funcs[seg_state][this_ind]
            envelope[i_env]=pt; i_env += 1
            inc_counter!(v_idx_visit, seg_state)
        else
            break
        end
    end

    # Equal values of x in the beginning ('regular' points)
    n_iter=1
    reached_maxlen = false
    while  !reached_maxlen && (get_x(pl_funcs[1], v_idx_visit[1]) == get_x(pl_funcs[2], v_idx_visit[2]))
        if get_y(pl_funcs[1], v_idx_visit[1]) > get_y(pl_funcs[2], v_idx_visit[2])
            envelope[i_env]=pl_funcs[1][v_idx_visit[1]]; i_env += 1
        else
            envelope[i_env]=pl_funcs[2][v_idx_visit[2]]; i_env += 1
        end

        reached_maxlen = v_idx_visit[1] == length(pl_funcs[1]) || v_idx_visit[2] == length(pl_funcs[2])

        if !reached_maxlen
            i_env += segment_intersect!(envelope, i_env,
                                        pl_funcs[1][n_iter], pl_funcs[1][n_iter+1],
                                        pl_funcs[2][n_iter], pl_funcs[2][n_iter+1])
            inc_counter!(v_idx_visit)
        end
        n_iter += 1
    end

    if reached_maxlen
        return envelope
    end

    # Irregular x values
    @inbounds seg_state = determine_first(pl_funcs, v_idx_visit)
    @inbounds while (v_idx_visit[1] <= length(pl_funcs[1])) && (v_idx_visit[2] <= length(pl_funcs[2]))
        @inbounds i_curr = v_idx_visit[seg_state]
        @inbounds i_other = v_idx_visit[other(seg_state)]

        @inbounds pt_curr = pl_funcs[seg_state][i_curr]
        @inbounds pt_other = pl_funcs[other(seg_state)][i_other]
        @inbounds pt_other_prev = pl_funcs[other(seg_state)][i_other-1] # TODO: Could this be saved from prev iteration?

        sco = inner_rotate(pt_curr - pt_other_prev, pt_other - pt_other_prev)

        if sco > 0
            envelope[i_env]=pt_curr; i_env += 1
        end

        if i_curr + 1 > length(pl_funcs[seg_state])
            seg_state = other(seg_state)
            break
        end

        pt_curr_next = pl_funcs[seg_state][i_curr+1]

        i_env += segment_intersect!(envelope, i_env, sco, pt_other_prev, pt_other, pt_curr, pt_curr_next)
        inc_counter!(v_idx_visit, seg_state)

        if get_x(pt_curr_next) >= get_x(pt_other)
            seg_state = other(seg_state)
        end

    end

    # Now finalize by adding all rightmost nodes
    while v_idx_visit[seg_state] <= length(pl_funcs[seg_state])
        envelope[i_env]=pl_funcs[seg_state][v_idx_visit[seg_state]]; i_env += 1
        inc_counter!(v_idx_visit, seg_state)
    end

    # Remove NaN entries if specified
    remove_nan && resize!(envelope, i_env-1)
    return envelope
end


compute_envelope(tup_pts1, tup_pts2, remove_nan=true)=begin
    Tuple(compute_envelope((PiecewiseLinear(tup_pts1...), PiecewiseLinear(tup_pts2...)),
                           remove_nan))
end




# * Auxiliary functions
" Returns 1 if `i==2`, and 2 if `i==1`."
other(i)=3-i

inc_counter!(v_idx, i)=v_idx[i] += 1
inc_counter!(v_idx)=begin
    v_idx[1] += 1; v_idx[2] += 1;
end


"""
Check that the x value is lower for the current than the alternative state at
indexes given by v_idx_visit.
"""
function is_lower_x(pl_funcs, seg_state, v_idx_visit)
    seg_other = other(seg_state)

    @inbounds i_curr = v_idx_visit[seg_state]
    @inbounds i_other = v_idx_visit[seg_other]
    get_x(pl_funcs[seg_state], i_curr) < get_x(pl_funcs[seg_other], i_other)
end


Base.@propagate_inbounds function determine_first(pl_funcs, init_ind=(1, 1))
    r = 1
    if get_x(pl_funcs[2], init_ind[2]) < get_x(pl_funcs[1], init_ind[1])
        r = 2
    end
    return r
end

"""
    segment_intersect!(v_store, pos_store, sco, ptA1, ptA2, ptB1, ptB2)

If an intersection between pl_funcs [A1, A2] and [B1, B2] exist, insert that in
v_store at position pos_store. Also, returns Boolean for whether the
intersection exists.

Note that ptA1, ptA2, ptB1 and ptB2 are NOT general points in the plane.
Because of the problem's structure, the following inequality holds:
get_x(ptA1) <= get_x(pt_B1)

Note also that `sco` already stores a value whose sign depends on the half space
of point ptB1 relative to the line formed by the segment ptA1-ptA2 Hence, the
first way of ruling out the intersection is to compute the same value for ptB2;

Formally, `sco` is the inner product of the 90 deg counter-clockwise
rotation of `P` (as defined below) with `B1 - A1`.

"""
function segment_intersect!(v_store, pos_store, sco, ptA1, ptA2, ptB1, ptB2)

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

    if r
        v_store[pos_store]=intersection_candidate
    end

    return r
end

segment_intersect!(v_store, pos_store, ptA1, ptA2, ptB1, ptB2)=begin
    segment_intersect!(v_store, pos_store,
                       inner_rotate(ptB1-ptA1, ptA2-ptA1),
                       ptA1, ptA2, ptB1, ptB2)
end
