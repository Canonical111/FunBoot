using HDF5

# """
#     my_compute_function(x)

# Compute some function of `x`.
# """

struct DeltaTable{T <: Real}
    dtable::Vector{T} #dimension table
    tablelength::Int
    grid::Int #The spacing between different grid is 1/grid
    dmax::T#last element in dtable
end 

"""

OneDimFunctionalTable()


struct OneDimFunctionalTable{T <: Real}
    dphi::T #The higher dim dphi, corresponds to dphi/2 in one dimension. 
    dtable::DeltaTable{T}
    order::Int #α/β_0 to α/β_{order-1} 
    functionalarray::Array{T, 4} #[deltaind, functionalind, alphabeta, plusminus]
    normplus::Vector{T}
    normminus::Vector{T}
    grid::Int #The spacing between different grid is 1/grid
    type::String # "F", "B" or "FB" 
end

"""

struct OneDimFunctionalTable{T <: Real}
    dphi::T #The higher dim dphi, corresponds to dphi/2 in one dimension. 
    dtable::DeltaTable{T}
    order::Int #α/β_0 to α/β_{order-1} 
    functionalarray::Array{T, 4} #[deltaind, functionalind, βB αB βF αF, +-]
    normplus::Vector{T}
    normminus::Vector{T}
    grid::Int #The spacing between different grid is 1/grid
    type::String # "F", "B" or "BF"
    normplusitp::Any
    normminusitp::Any 
end

function griddtable(dtable::Vector{T}) where T <: Real
    length(dtable) < 2 && throw(ArgumentError("dtable must have at least 2 elements"))
    diff = dtable[2] - dtable[1]
    if iszero(diff)
        throw(DomainError(diff, "Cannot divide by zero"))
    end
    return trunc(Int, Float64(1 / diff))
end

function readsinglepoint(hdf5file::String, type::Type=Float64)
    
    tmpdtable::Vector{type}= h5read(hdf5file, "dtable")
    tmpdphi::type= h5read(hdf5file, "dphi")
    tmparray::Array{type, 4}= h5read(hdf5file, "functionaltable")
    newdtable=DeltaTable(tmpdtable, length(tmpdtable), griddtable(tmpdtable), last(tmpdtable) )
    normplus::Vector{type}= h5read(hdf5file, "normplus")
    normminus::Vector{type}= h5read(hdf5file, "normminus")

    normplusitp=interpolate(2*tmpdtable,  normplus./normplus[1], BSplineOrder(8)) 
    normminusitp=interpolate(2*tmpdtable,  normminus./normminus[1], BSplineOrder(8)) 
    newfunctionaltable=OneDimFunctionalTable(tmpdphi, newdtable,  size(tmparray, 2), tmparray, normplus, normminus, griddtable(tmpdtable), "BF", normplusitp, normminusitp)
    return newdtable, newfunctionaltable
end
