using IterTools
import Base: copy

mutable struct LinearProgramCompact{T <: Real}
    coeffs::Vector{T}        # the coefficients in the current solution
    cost::T
    costprevious::T
    costvarmax::T
    err::Int
    vecMatrix::Matrix{T}      # discrete vectors
    AMatrix::Matrix{T}
    
    invA::LU{T, Matrix{T}, Vector{Int64}}                # this is basically the set of vectors in the solution, inverted
    target::Vector{T}            # typically the identity vector
    
    functional::Vector{T}       # the current linear functional, \simeq costs.A^{-1}.
    num_of_functional::Int

    solcostlist::Vector{T}
    veccostlist::SparseVector{T, Int64}

    rclist::Vector{T}
    negindlist::Vector{Int}
    negrclist::Vector{T}
    negnum::Int
    AUXnum::Int

    bbaroutindlist::Vector{Int}
    costvar::Vector{T}

    bbarin::Tuple{Int, Tuple{T, Int}}
    bbarout::Tuple{Int, Tuple{T, Int}}

    gridsol::Vector{Tuple{T, Int}}
    gridVectors::Vector{Tuple{T, Int}}

    label::String               # A description of this linear problem
    
    status::String
    functionallist::Vector{String}
    spinlist::Vector{Int}
    gapchannel::Int
end

function Base.copy(lp::LinearProgramCompact)
    return LinearProgramCompact(
        lp.coeffs,  # reference to the original vector
        lp.cost,  # copy of the scalar value
        lp.costprevious,
        lp.costvarmax,
        lp.err,
        lp.vecMatrix,  # reference to the original matrix
        lp.AMatrix,
        lp.invA,  # reference to the original LU factorization
        lp.target,  # reference to the original vector
        lp.functional,  # reference to the original vector
        lp.num_of_functional,
        lp.solcostlist,  # reference to the original vector
        lp.veccostlist,  # reference to the original SparseVector
        lp.rclist,  # reference to the original vector
        lp.negindlist,  # reference to the original vector
        lp.negrclist,  # reference to the original vector
        lp.negnum,
        lp.AUXnum,
        lp.bbaroutindlist,  # reference to the original vector
        lp.costvar,  # reference to the original vector
        lp.bbarin,  # copy of the Tuple
        lp.bbarout,  # copy of the Tuple
        lp.gridsol,  # reference to the original vector
        lp.gridVectors,  # reference to the original vector
        lp.label,  # copy of the String
        lp.status,  # copy of the String
        lp.functionallist,  # reference to the original vector
        lp.spinlist,  # reference to the original vector
        lp.gapchannel
    )
end

function convert_type(lp::LinearProgramCompact{Float64}, ::Type{T}) where {T <: Real}
    lp_new = LinearProgramCompact{T}(
        T.(lp.coeffs),
        T(lp.cost),
        T(lp.costprevious),
        T(lp.costvarmax),
        lp.err,
        T.(lp.vecMatrix),
        T.(lp.AMatrix),

        lu(T.(lp.invA.factors)), # We assume that invA is always LU. If it's not, you need to handle different cases.
        T.(lp.target),
        T.(lp.functional),
        lp.num_of_functional,
        
        T.(lp.solcostlist),
        SparseVector(lp.veccostlist.n, lp.veccostlist.nzind, T.(lp.veccostlist.nzval)), # convert SparseVector
        
        T.(lp.rclist),
        lp.negindlist,
        T.(lp.negrclist),
        lp.negnum,
        lp.AUXnum,
        
        lp.bbaroutindlist,
        T.(lp.costvar),
        
        (lp.bbarin[1], (T(lp.bbarin[2][1]), lp.bbarin[2][2])),
        (lp.bbarout[1], (T(lp.bbarout[2][1]), lp.bbarout[2][2])),
        
        Tuple{T, Int}[(T(x[1]), x[2]) for x in lp.gridsol],
        Tuple{T, Int}[(T(x[1]), x[2]) for x in lp.gridVectors],

        lp.label,
        lp.status,
        lp.functionallist,
        lp.spinlist,
        lp.gapchannel
    )
    
    return lp_new
end


function addonepoint!(l, delta::T, lp::LinearProgramCompact{T}, table2d::HigherDimFunctionalTable{T}) where T
    lp.vecMatrix=hcat(lp.vecMatrix, functionalop((delta, l), table2d))
    lp.veccostlist=vcat(lp.veccostlist, 0)
    lp.rclist=push!(lp.rclist, zero(T))
    push!(lp.gridVectors, (delta, l))
end


function filteroutop!(lp::LinearProgramCompact{T}) where T
    indfilterout=findall( x->x[2]==-1, lp.gridVectors)
    union!(indfilterout, indAUX)
    lp.vecMatrix=lp.vecMatrix[:, 1:end .∉ [indfilterout]]
    lp.veccostlist=lp.veccostlist[1:end .∉ [indfilterout]]
    lp.rclist=lp.rclist[1:end .∉ [indfilterout]] 
    deleteat!(lp.gridVectors, indfilterout)
end

function filteroutop!(f, lp::LinearProgramCompact{T}, filterAUX=false) where T
    
    indfilterout=findall( f, lp.gridVectors)
    if filterAUX
        indAUX=findall( x->x[2]==-1, lp.gridVectors)
        union!(indfilterout, indAUX)
    end

    lp.vecMatrix=lp.vecMatrix[:, 1:end .∉ [indfilterout]]
    lp.veccostlist=lp.veccostlist[1:end .∉ [indfilterout]]
    lp.rclist=lp.rclist[1:end .∉ [indfilterout]] 
    deleteat!(lp.gridVectors, indfilterout)
end

# function addonepoint!(itp, l, delta::T, lp::LinearProgramCompact{T}) where T

#     lp.vecMatrix=hcat(lp.vecMatrix, )
#     lp.veccostlist.n+=1
#     lp.rclist=push!(lp.rclist, zero(T))
#     push!(lp.gridVectors, (delta, l))
# end


function filterfunctional!(indlist::Vector{Int}, lp::LinearProgramCompact{T}) where T
    @assert lp.status=="Initialized"
    lp.coeffs=lp.coeffs[indlist]
    lp.vecMatrix=lp.vecMatrix[indlist, :]
    lp.AMatrix=lp.AMatrix[indlist, indlist]
    lp.invA=lu(lp.AMatrix)
    lp.target=lp.target[indlist]
    lp.functional=lp.functional[indlist]
    lp.solcostlist=lp.solcostlist[indlist]
    lp.AUXnum=length(indlist)
    lp.gridsol=lp.gridsol[indlist]
    lp.num_of_functional=length(indlist)
    lp.functionallist=lp.functionallist[indlist]
end

function sign_with_zero(x::T) where T
    if x === -zero(T)
        return -one(T)
    elseif x === zero(T)
            return one(T)
        else
        return sign(x)
    end
end

function initialA(target::Vector{T}) where T <: Real
    mat=Matrix{T}(I, length(target), length(target))
    for i in eachindex(target)
        mat[i, i]= sign_with_zero(target[i])
    end
    mat
end

function convolveiterator(lambda::Int, grid::Int, biggrid::Int, dlength::Int)
    iter1=1:div(3*grid, 2)
    iter=chain(iter1, (iter1.stop+1):2:(iter1.stop+grid))
    for i in 3:2*lambda+3
        iter=chain(iter, (iter.it[end].stop+1):(2*i):(iter.it[end].stop+grid))
    end
    return filter!(x->x<=dlength, collect(chain(iter, (iter.it[end].stop+1):biggrid:div(dlength, 3))))
end

function initializeIndices(lengthlist::Vector{Int}, grid::Int, lambda::Int, biggrid::Int)
    # [i==length(lengthlist) ? convolveiterator(gridnum*twistgap, lengthl, biggrip) : collect(1:gridnum*twistgap) for (i, lengthl) in enumerate( lengthlist )]
    [convolveiterator(lambda, grid, biggrid, i)  for i in lengthlist ]
end

function costf(l::Int, delta::T, spintoGap::Int) where T
    if l!=spintoGap
        return zero(delta)
    else  
        cost=T( 1e50 * 10^(-20* (Float64(delta)-l)) )
        # return zero(delta)
        return abs(cost) > T(1e-250) ? cost : zero(delta)
    end
end

# function costf(l::Int, delta::T, spintoGap::Int) where T
#     if l!=spintoGap
#         return zero(delta)
#     else  
#         cost=zero(T)
#         # return zero(delta)
#         return abs(cost) > T(1.) ? zero(cost) : cost
#     end
# end

function solutionRaw(lp::LinearProgramCompact{T}) where T
    println("The solution, with the OPE normalized")
    println("|  L  |     Δ       |      OPE     |")
    println("|:---:|:-----------:|:------------:|") 
    map(x->@printf("| %3d | %11.7f | %7.6e |\n", x[1][2], x[1][1], x[2]),[(lp.gridsol[i], lp.coeffs[i]) for i=1: lp.num_of_functional] |>(x->sort(x,by=x->x[1])));
    println("\n")
end



function initializeLP(functionaltable::HigherDimFunctionalTable{T}, ltable::Vector{Int}, lambda::Int, biggrid::Int=200; spintoGap::Int=0) where T <: Real        
    @assert lambda<=functionaltable.order
    @assert issubset(ltable, functionaltable.ltable)
    
    lindex::Vector{Int}=indexin(ltable, functionaltable.ltable)
    functionalnum=length(functionalpairs(lambda, functionaltable.type))
    indlist=initializeIndices(functionaltable.dtablelengthlist[lindex], functionaltable.grid, lambda, biggrid)
    push!(indlist[1], functionaltable.dtablelengthlist[lindex][1])
    coeffs=abs.(functionaltable.target[1:functionalnum])
    vecMatrix=hcat([transpose(functionaltable[spin][indlist[i],1:functionalnum]) for (i, spin) in enumerate(ltable)]...)
    gridVectors=vcat([functionaltable.gridlist[lindex[i]][indlist[i]] for (i, spin) in enumerate(ltable)]...)
    
    
    vecnum=size(vecMatrix, 2)
    
    
    target=functionaltable.target[1:functionalnum]
    # target+= 1e-16 * randn(length(target))
    AMatrix=initialA(target[1:functionalnum])
    invA=lu(AMatrix)
    functional=zeros(T, functionalnum)
    
    solcostlist=fill(convert(T, 1e128), functionalnum)
    #veccostlist=[l==2 ? T(10^(-3*Float64(delta-l)) ) : zero(T) for (delta, l) in gridVectors]|>sparse
    veccostlist=[costf(l, delta, spintoGap) for (delta, l) in gridVectors]|>sparse
    rclist=Vector{T}(undef, vecnum)
    negindlist=Vector{Int}(undef, 0)
    negrclist=Vector{T}(undef, 0)
    negnum=0

    bbaroutindlist=Vector{Int}(undef, 0)
    costvar=Vector{T}(undef, 0)
    
    gridsol=[(T(-i), -1) for i in 1:functionalnum]
    
    lable="A fucking stupid LP."
    LinearProgramCompact(coeffs,
                         zero(T),
                         zero(T),
                         zero(T),
                         zero(Int),
                         vecMatrix,
                         AMatrix,
                         invA,
                         target,
                         functional,
                         functionalnum,
                         solcostlist,
                         veccostlist,  
                         rclist,
                         negindlist,
                         negrclist,
                         negnum,
                         functionalnum, 
                         bbaroutindlist,
                         costvar,
                         (0, (zero(T), 0)), 
                         (0, (zero(T), -1)), 
                         gridsol,
                         gridVectors,
                         lable,
                         "Initialized",
                         functionalstringpairs(lambda, functionaltable.type),
                         ltable,
                         spintoGap)
    
end
