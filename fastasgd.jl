"""
Copyright (C) 2018  Mark Ponder

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
#-----------------------------------------------------------------------------#
# Spatially and Locally Adaptive Sparse Grid Approximation
# This program is derived from two sources
# [1] "Local and Dimension Adaptive Sparse Grid Interpolation and Quadrature" by
# J.D. Jakeman and S.G. Robert
# and
# [2] "Compact Data Structure and Scalable Algorithms for the Sparse Grid Technique"
# by Pfluger, et. al.
# The original code that this is based off of is "fastsg" and can be found here
# https://github.com/afm/fastsg
# Following [2], I construct a bijection so that the hierarchal coefficients
# are stored in a one dimensional array. Rather than use a regular sparse grid
# and then projecting lower dimensional sparse grids onto the boundary, I
# setup my bijection to use Curtis-Clenshaw grid points
# Following [1], the algorithm is both spatially and locally adaptive. At each
# step, I take the level with the highest error and check its children for
# admissibility.
#-----------------------------------------------------------------------------#

module fastasgd

export ASG, adaptivegrid, iinterpolate

type ASG
    d::Int64
    ds::Int64
    l::Int64
    sg1d::Array{Float64, 2}
    numOfGridPoints::Int64
    lb::Array{Float64, 1}
    ub::Array{Float64, 1}
    wrk::Array{Int64, 1}

    function ASG( d::Int64, ds::Int64, l::Int64, lb, ub )

        d < 1 && return error()
        l < 0 && return error()
        ds < 1 && return error()

        numOfGridPoints = gsize(d, l)
        sg1d = zeros(Float64, ds, numOfGridPoints)
        wrk = zeros(Int, d)
        return new(d, ds, l, sg1d, numOfGridPoints, lb, ub, wrk)
    end
end

#------------------------------------------------------------------------------#
# Helper Functions
#------------------------------------------------------------------------------#
# Find the next level
# levels are ordered as follows
# [3,0] -> [2,1] -> [1,2] -> [0,3] -> [4,0], etc.
function lnext(l̄::Array{Int,1})
    l̄[end] == sum(l̄) && begin
        l̄[1] = sum(l̄) + 1
        for i = 2:length(l̄)
            l̄[i] = 0
        end
        return nothing
    end
    t::Int = 1

    while (l̄[t] == 0)
        t += 1
    end

    l̄[1] = l̄[t] - 1
    t > 1 && (l̄[t] = 0)
    l̄[t+1] = l̄[t+1] + 1
    return nothing
end

# Same thing, mutating different array
function lnext(l̄::Array{Int,1}, r̄::Array{Int,1})
    l̄[end] == sum(l̄) && begin
        r̄[1] = sum(l̄) + 1
        for i = 2:length(l̄)
            r̄[i] = 0
        end
        return nothing
    end
    for i = 1:length(l̄)
        r̄[i] = l̄[i]
    end
    t::Int = 1

    while (l̄[t] == 0)
        t += 1
    end

    r̄[1] = l̄[t] - 1
    t > 1 && (r̄[t] = 0)
    r̄[t+1] = l̄[t+1] + 1
    return nothing
end
#------------------------------------------------------------------------------#
# Bijection
#------------------------------------------------------------------------------#
# Curtis-Clenshaw type bijection. Takes level and index of a point and returns
# the points position in the 1D coefficients array
#------------------------------------------------------------------------------#
# Bijection
#------------------------------------------------------------------------------#
# Curtis-Clenshaw type bijection. Takes level and index of a point and returns
# the points position in the 1D coefficients array
function cc_gp2idx( levels::Array{Int64,1}, indices::Array{Int64, 1}, wrk::Array{Int64, 1}, d::Int )
    index2 = 0
    prod = 1
    index1 = indices[1]
    for i = 1:d-1
        index1 = index1 * (levels[i+1] > 1 ? 2^(levels[i+1]-1) : 2^levels[i+1]) + indices[i+1]
        # index will always be their papers index + 2^sum(level)
    end

    index3 = 0
    for i = 0:Base.sum(levels)-1
        wrk[1] = i
        for j = 2:d
            wrk[j] = 0
        end

        while (wrk[end] != i)

            prod = 1
            for j = 1:d
                wrk[j] < 2 && begin prod *= 2^wrk[j] end
                wrk[j] > 1 && begin prod *= 2^(wrk[j]-1) end
            end
            index3 += prod
            lnext( wrk )
        end

        index3 += (i > 1 ? 2^(i-1) : 2^i)
    end

    wrk[1] = Base.sum(levels)
    for i = 2:d
        wrk[i] = 0
    end

    while (levels != wrk)

        prod = 1
        for j = 1:d
            wrk[j] < 2 && begin prod *= 2^wrk[j] end
            wrk[j] > 1 && begin prod *= 2^(wrk[j]-1) end
        end
        index2 += prod
        lnext( wrk )
    end

    return index1 + index2 + index3
end



function cc_idx2gp( index::Int64, levels::Array{Int64,1}, indices::Array{Int64,1}, wrk::Array{Int64,1}, d::Int64 )

    prod = 1
    index3 = 0
    lidx = 0

    sum = 0

    while ( true )
        wrk[1] = sum
        for i = 2:d
            wrk[i] = 0
        end
        lidx = 0

        while (wrk[end] != sum)
            prod = 1
            for j = 1:d
                wrk[j] < 2 && begin prod *= 2^wrk[j] end
                wrk[j] > 1 && begin prod *= 2^(wrk[j]-1) end
            end
            lidx += prod
            lnext( wrk )

        end

        lidx += (sum > 1 ? 2^(sum-1) : 2^sum)

        index < index3 + lidx && break
        index3 += lidx
        sum += 1

    end

    index -= index3

    wrk[1] = sum
    for i = 2:d
        wrk[i] = 0
    end

    index2 = 0

    while ( true )
        prod = 1
        for j = 1:d
            wrk[j] < 2 && begin prod *= 2^wrk[j] end
            wrk[j] > 1 && begin prod *= 2^(wrk[j]-1) end
        end

        index < index2 + prod && break

        index2 += prod
        lnext( wrk )
    end

    for i = 1:d
        levels[i] = wrk[i]
    end

    index -= index2

    for i = d-1:-1:1
        didx = index % (levels[i+1] > 1 ? 2^(levels[i+1]-1) : 2^levels[i+1])
        indices[i+1] = didx
        index = div(index, (levels[i+1] > 1 ? 2^(levels[i+1]-1) : 2^levels[i+1]))

    end

    indices[1] = index

    return 0

end

#------------------------------------------------------------------------------#
# Moving between grid points and level/index
#------------------------------------------------------------------------------#
function coord2li( coord::Array{Float64,1}, levels::Array{Int64, 1}, indices::Array{Int64, 1}, d::Int64 )
    for i = 0:d-1
        if  (coord[i+1] != 0.0)&(coord[i+1] != 1.0)&(coord[i+1] != 0.5)
            levels[i+1] = 0
            cc = coord[i+1]
            while (cc != floor(cc))
                cc *= 2
                levels[i+1] += 1
            end
            indices[i+1] = Int( floor( (cc - 1.0)/2.0 ) )
        elseif coord[i+1] == 0.0
            levels[i+1] = 1
            indices[i+1] = 0
        elseif coord[i+1] == 0.5
            levels[i+1] = 0
            indices[i+1] = 0
        elseif coord[i+1] == 1.0
            levels[i+1] = 1
            indices[i+1] = 1
        end
    end
    return nothing
end

function coord2li( coord::Float64, levels::Array{Int64, 1}, indices::Array{Int64, 1}, cd::Int64 )
    i = cd
    if  (coord != 0.0)&(coord != 1.0)&(coord != 0.5)
        levels[i+1] = 0
        cc = coord
        while (cc != floor(cc))
            cc *= 2
            levels[i+1] += 1
        end
        indices[i+1] = Int( floor( (cc - 1.0)/2.0 ) )
    elseif coord == 0.0
        levels[i+1] = 1
        indices[i+1] = 0
    elseif coord == 0.5
        levels[i+1] = 0
        indices[i+1] = 0
    elseif coord == 1.0
        levels[i+1] = 1
        indices[i+1] = 1
    end
    return nothing
end

function li2coord(levels::Array{Int64, 1}, indices::Array{Int64, 1}, coords::Array{Float64, 1}, d::Int64)
    for i = 0:d-1
        if levels[i+1] > 1
            coords[i+1] = 1.0 / 2.0^(levels[i+1]-1) * (indices[i+1] + 0.5)
        elseif levels[i+1] == 1
            indices[i+1] == 0 && (coords[i+1] = 0.0)
            indices[i+1] != 0 && (coords[i+1] = 1.0)
        else
            coords[i+1] = 0.5
        end
    end
    return nothing
end

function cc_idx2gp( index::Int64, coords::Array{Float64, 1}, d::Int64 )
    levels = zeros(Int, d)
    indices = zeros(Int, d)
    wrk = zeros(Int, d)
    cc_idx2gp(index, levels, indices, wrk, d);
    li2coord(levels, indices, coords, d);
end

#------------------------------------------------------------------------------#
# Functions to find children points
#------------------------------------------------------------------------------#
function gp2LeftChild( levels::Array{Int64, 1}, indices::Array{Int64, 1}, gp::Array{Float64, 1},
        pgp::Array{Float64, 1}, ord::Int64, cd::Int64, d::Int64 )
    div = (ord + 1 > 1 ? 2^(ord) : 2^(ord+1) )
    li2coord(levels, indices, gp, d)
    for i = 1:d
        pgp[i] = gp[i]
    end
    gp[cd+1] - 1.0/div <  0.0 && return -1
    pgp[cd+1] = gp[cd+1] - 1.0/div
    return 0
end

function gp2RightChild( levels::Array{Int64, 1}, indices::Array{Int64, 1}, gp::Array{Float64, 1},
        pgp::Array{Float64, 1}, ord::Int64, cd::Int64, d::Int64 )
    div = (ord + 1 > 1 ? 2^(ord) : 2^(ord+1) )
    li2coord(levels, indices, gp, d)
    for i = 1:d
        pgp[i] = gp[i]
    end
    gp[cd+1] + 1.0/div >  1.0 && return -1
    pgp[cd+1] = gp[cd+1] + 1.0/div
    return 0
end

# Count Levels

function gsize(d::Int64, n::Int64)
    index3 = 0
    for i = 0:n
        wrk = zeros(Int, d)
        wrk[1] = i
        while (wrk[end] != i)
            prod = 1
            for j = 1:d
                wrk[j] < 2 && begin prod *= 2^wrk[j] end
                wrk[j] > 1 && begin prod *= 2^(wrk[j]-1) end
            end
            index3 += prod
            lnext( wrk )

        end
        index3 += (i > 1 ? 2^(i-1) : 2^i)
    end
    return index3
end

#------------------------------------------------------------------------------#
# Basis Functions
#------------------------------------------------------------------------------#
# linear Basis Function
function basis(left::Float64, right::Float64, coord::Float64, level::Int64)
    if level > 1
        mid = (left + right)/2.0
        dif = abs(coord - mid)
        m_i = 2.0^level + 1.0
        return 1.0 - (m_i-1.0)*dif
    elseif level == 1
        dif = abs(coord - 0.5)
        return 2.0*dif
    else
        return 1.0
    end
end

# Derivative of Linear Basis Function
function dbasis(left::Float64, right::Float64, coord::Float64, level::Int64)
    if level > 1
    mid = (left + right)/2.0
    dif = sign(coord - mid)
    dif == 0.0 && (dif = 1.0)
    m_i = 2.0^level + 1.0
    return -(m_i-1.0)*dif
    elseif level == 1
        dif = sign(coord - 0.5)
        dif == 0.0 && (dif = 1.0)
        return 2.0*dif
    else
        return 0.0
    end
end

# Integral Basis Function
function ibasis(left::Float64, right::Float64, coord::Float64, level::Int64)
    if level > 1
        mid = (left + right)/2.0
        dif = abs(coord - mid)
        m_i = 2.0^level
        if coord <= mid
            lshift = left - m_i*(mid - left*0.5)*left
            val = coord - m_i*(mid - coord*0.5)*coord - lshift
        else
            lshift = -left + m_i*(mid - left*0.5)*left + mid - m_i*(mid - mid*0.5)*mid
            val =   lshift  + coord - m_i*(coord*0.5 - mid)*coord -
                    mid +  m_i*(mid*0.5 - mid)*mid
        end
        return val
    elseif level == 1
        if coord <= 0.5
            val = (0.5 * coord - coord^2/2.0) * 2.0
        else
            val = (coord^2/2.0 - 0.5 * coord) * 2.0 + 0.25
        end
        return val
    else
        return coord
    end
end


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Interpolating Functions
# interpolate -> Interpolation
# dinterpolate -> Derivative of Interpolation
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
function interpolate( coord::Array{Float64, 1},  a::ASG, depth::Int64, state::Int64 )
    d = a.d
    for i = 1:d
        coord[i] = (coord[i] - a.lb[i])/(a.ub[i] - a.lb[i])
    end
    sg1d = view(a.sg1d, state,:)
    wrk = a.wrk
    res = 0.0
    index2 = 0
    for ord = 0:depth
        wrk[1] = ord
        for i = 2:d
            wrk[i] = 0
        end

        while ( true )
            prod = 1.0
            lsum = 1
            index1 = 0
            for j = 1:d
                div = (wrk[j] > 1 ? 2^(wrk[j]-1) : 2^wrk[j])

                if (0.0 <= coord[j]) & (coord[j] <= 1.0)
                  index1 = Int(index1 * div + floor( coord[j]*div ) - (coord[j] == 1.0))
                    left = floor( coord[j]*div  ) / div - (coord[j] == 1.0) / div
                elseif coord[j] < 0.0
                    index1 = Int(index1 * div )
                    left = 0.0
                else
                    index1 = Int(index1 * div + floor( div ) - 1.0)
                    left = floor( div  ) / div - 1.0 / div
                end

                right = left + 1.0/div
                prod *= basis(left, right, coord[j], wrk[j])
                lsum *= (wrk[j] > 1 ? 2^(wrk[j]-1) : 2^wrk[j])
            end

            res += prod * sg1d[index1 + index2 + 1]
            index2 += lsum
            wrk[end] == ord && break
            lnext( wrk )

        end
    end
    res
end


function dinterpolate( coord::Array{Float64, 1},  a::ASG, depth::Int64, state::Int64, dim::Int64 )
    d = a.d
    sg1d = view(a.sg1d, state, :)
    wrk = a.wrk
    for i = 1:d
        coord[i] = (coord[i] - a.lb[i])/(a.ub[i] - a.lb[i])
    end
    res = 0.0
    index2 = 0
    for ord = 0:depth
        wrk[1] = ord
        for i = 2:d
            wrk[i] = 0
        end
        while ( true )
            prod = 1.0
            lsum = 1
            index1 = 0
            for j = 1:d
                div = (wrk[j] > 1 ? 2^(wrk[j]-1) : 2^wrk[j])
                if (0.0 <= coord[j]) & (coord[j] <= 1.0)
                  index1 = Int(index1 * div + floor( coord[j]*div ) - (coord[j] == 1.0))
                    left = floor( coord[j]*div  ) / div - (coord[j] == 1.0) / div
                elseif coord[j] < 0.0
                    index1 = Int(index1 * div )
                    left = 0.0
                else
                    index1 = Int(index1 * div + floor( div ) - 1.0)
                    left = floor( div  ) / div - 1.0 / div
                end
                right = left + 1.0/div
                j == dim && (prod *= dbasis(left, right, coord[j], wrk[j])/ (a.ub[j] - a.lb[j]))
                j != dim && (prod *= basis(left, right, coord[j], wrk[j]))
                lsum *= (wrk[j] > 1 ? 2^(wrk[j]-1) : 2^wrk[j])
            end

            res += prod * sg1d[index1 + index2 + 1]
            index2 += lsum
            wrk[end] == ord && break
            lnext( wrk )

        end

    end
    res
end

function iinterpolate( coord::Array{Float64, 1},  a::ASG, depth::Int64, state::Int64)
    d = a.d
    sg1d = view(a.sg1d, state,:)
    wrk = a.wrk
    levels = zeros(Int, d)
    indices = zeros(Int, d)
    coords = zeros(Float64, d)

    for i = 1:d
        coord[i] = (coord[i] - a.lb[i])/(a.ub[i] - a.lb[i])
    end

    res = 0.0

    for i = 1:length(sg1d)
        sg1d[i] == 0 && continue
        cc_idx2gp(i-1, levels, indices, wrk, d )
        li2coord(levels, indices, coords, d)
        prod = 1.0

            for j = 1:d
                div = (levels[j] > 1 ? 2^(levels[j]-1) : 2^levels[j])
                if (0.0 <= coords[j]) & (coords[j] <= 1.0)
                    left = floor( coords[j]*div  ) / div - (coords[j] == 1.0) / div
                elseif coords[j] < 0.0
                    left = 0.0
                else
                    left = floor( div  ) / div - 1.0 / div
                end
                right = left + 1.0/div

                if (coord[j] >= left)
                    prod *= ibasis(left, right, min(right, coord[j]), levels[j]) * (a.ub[j] - a.lb[j])
                else
                    prod *= 0.0
                end
            end
            res += prod * sg1d[i]

    end
    res
end

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Finds all children points for candidate point and determines which ones
# will be active and which will be redundant
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
function createGrid( l::Array{Int,1}, sg1d::Array{Float64, 2}, Inidx::Array{Int64, 1},
    d::Int64, ϵ::Float64, f::Function, a::ASG, state::Int64)
    plevels = zeros(Int, d)
    for i = 1:d
        plevels[i] = l[i]
    end
    lb = a.lb
    ub = a.ub
    pindices = zeros(Int, d)
    levels = zeros(Int, d)
    indices = zeros(Int, d)
    wrk = zeros(Int, d)
    lwrk = zeros(Int, d)
    gp = zeros(Float64, d)
    pgp = zeros(Float64, d)

    for n = 1:d
        plevels[n] == 0 && continue
        plevels[n] -= 1
        t1 = cc_gp2idx(plevels, pindices, wrk, d) + 1


        lnext(plevels, lwrk)

        t2 = cc_gp2idx(lwrk, pindices, wrk, d)

        #lprev(plevels)

        for j = t1:t2

            Inidx[j] != 0 && begin

                cc_idx2gp( j - 1, plevels, pindices, wrk, d)

                for i = 0:d-1

                gp2RightChild(plevels, pindices, gp, pgp, plevels[i+1] + 1 , i, d) != - 1 && begin
                    coord2li(pgp, levels, indices, d)

                    id = cc_gp2idx(levels, indices, wrk, d) + 1

                    Inidx[id] == 1 && continue
                    for c = 1:d
                        pgp[c] = pgp[c]*(ub[c] - lb[c]) + lb[c]
                    end
                    fval = f(pgp)
                    sg1d[state, id] = fval - interpolate( pgp, a, sum(levels)-1, state );
                    (sum(levels) <= 4) | (abs(sg1d[state, id]) >= ϵ) && (Inidx[id] = 1) # /abs(fval)
                end

                gp2LeftChild(plevels, pindices, gp, pgp, plevels[i+1] + 1, i, d) != -1 && begin
                    coord2li(pgp, levels, indices, d)
                    id = cc_gp2idx(levels, indices, wrk, d) + 1
                    Inidx[id] == 1 && continue
                    for c = 1:d
                       pgp[c] = pgp[c]*(ub[c] - lb[c]) + lb[c]
                    end
                    fval = f(pgp)
                    sg1d[state, id] = fval - interpolate( pgp, a, sum(levels)-1, state );
                    (sum(levels) <= 4) | (abs(sg1d[state, id]) >= ϵ) && (Inidx[id] = 1)
                end

                end

            end

        end

        plevels[n] += 1
    end
end


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Generate Dimension and Local adaptive sparse grids
#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
function adaptivegrid( a::ASG, f::Function, lϵ::Float64, sϵ::Float64, state::Int64 )
    d = a.d
    l = a.l
    lb = a.lb
    ub = a.ub
    ng = a.numOfGridPoints
    adm = 1 # Admissibility check
    # Maximum number of indices based on depth
    totalcnt = 0
    for i = 0:l
        totalcnt += binomial(d - 1 + i, i)
    end

    Aidx = Int64[] # Active Indices location
    err = Float64[] # Errors of Active Indices

    Inidx = zeros(Int64, size(a.sg1d, 2)) # Include index
    Inl = zeros(Bool, size(a.sg1d, 2)) # Include level

    Inl[1] = 1 # Initialize first level
    Inidx[1] = 1 # Initialize first index

    # Working Arrays
    levels = zeros(Int, d)
    indices = zeros(Int, d)

    wrk = zeros(Int, d)
    lwrk = zeros(Int, d)
    # Initialize Function Evaluations
    for i = 1:size(a.sg1d, 2)
        a.sg1d[state, i] = 0.0
    end

    a.sg1d[state, 1] = f(ones(d)*0.5.*(ub - lb) + lb)

    # Intialize Errors and Counts
    push!(Aidx, 0)
    push!(err, abs(a.sg1d[state, 1]) + 1.0)
    currerr = 1.0
    count = 0
    r = 1.0 + abs(a.sg1d[state, 1])

    while (r > sϵ)*(count < totalcnt)

        si = sortperm(err, rev = true) # to determine maximum error
        err = err[si]
        Aidx = Aidx[si]

        aidx = shift!(Aidx)
        currerr = shift!(err)
        r -= currerr

        cc_idx2gp( aidx, levels, indices, wrk, d)

        for i = 1:d
            indices[i] = 0
        end

        # Increment in each dimension
        for k = 1:d
            levels[k] += 1
            sum(levels) > l && begin
                levels[k] -= 1
                continue
            end # Prevents evaluation past the max depth
            adm = 1
            # Check in every prior direction for admissibility
            for n = 1:d
                levels[n] == 0 && continue
                levels[n] -= 1
                idx = cc_gp2idx(levels, indices, wrk, d)
                adm *= (Inl[idx + 1] != 0) # Can make more efficent
                levels[n] += 1
            end

            adm == 1 && begin
                idx = cc_gp2idx(levels, indices, wrk, d)
                Inl[idx + 1] = 1
                push!(Aidx, idx)
                createGrid(levels, a.sg1d, Inidx, d, lϵ, f, a, state)
                lnext(levels, lwrk)

                t2 = cc_gp2idx(lwrk, indices, wrk, d)

                #lprev(levels)
                nerr = sum(abs(a.sg1d[idx+1:t2]) )
                push!(err,  nerr)
                r += nerr
            end
            levels[k] -= 1
        end
        count += 1

    end

    return 0
end
#-------------------------------------------------------------------#
# Wrapper Functions
#-------------------------------------------------------------------#
function(asg::ASG)(x::Array{Float64,1}, state::Int64)
    return interpolate( x , asg, asg.l, state )
end

function(asg::ASG)(x::Vararg{Union{Float64, Int64}})
    pt = zeros(Float64, asg.d)
    for i = 1:asg.d
        pt[i] = x[i]::Float64
    end
    return interpolate( pt , asg, asg.l, x[end]::Int64 )
end

# Wrapper for derivative interpolation
function(asg::ASG)(j::Int64, x::Array{Float64,1}, state::Int64)
    return dinterpolate( x , asg, asg.l, state, j )
end

# Wrapper for arbitrary derivative interpolation
function(asg::ASG)(j::Int64, x::Vararg{Union{Float64, Int64}})
    pt = zeros(Float64, asg.d)
    for i = 1:asg.d
        pt[i] = x[i]::Float64
    end
    return dinterpolate( pt , asg, asg.l, x[end]::Int64, j )
end


end
