using SpecialFunctions
import Base: getindex


mutable struct HigherDimFunctionalTable{T <: Real}
    dim::Int
    dphi::T #The higher dim dphi, corresponds to dphi/2 in one dimension. 
    ltable::Vector{Int}
    dtablelist::Vector{Vector{T}}
    dtablelengthlist::Vector{Int}
    gridlist::Vector{Vector{Tuple{T, Int}}}
    order::Int #α/β_0 to α/β_{order-1} 
    functionalDic::Dict{Int, Matrix{T}} #
    target::Vector{T}
    num_of_functional::Int
    grid::Int #The spacing between different grid is 1/grid
    type::String # "F", "B" or "FB"
    itp::Any 
    functionallist::Vector{String}
end


function generateintepolation(functionaltable::HigherDimFunctionalTable{T}, order::Int=10) where T
    return [[interpolate(functionaltable.dtablelist[findfirst(isequal(l), functionaltable.ltable)], functionaltable[l][:,i], BSplineOrder(order)) for i in 1:functionaltable.num_of_functional] for l in functionaltable.ltable]
end



getindex(ft::HigherDimFunctionalTable{T}, i::Int) where T=ft.functionalDic[i]

function lambdapairs(lambda::Int) 
    return sort([[i, j] for i in 1:lambda+1 for j in 1:lambda+1 if i + j <=lambda + 2 ], by=(x->(sum(x),x[1])))
end

function filterBF(BFtype::String)
    if BFtype=="F+F-"
        return (list-> filter(x-> (x[3] in [3, 4])&&(x[4] in [3, 4]),list))
    elseif BFtype=="B+B-"
        return (list-> filter(x-> (x[3] in [1, 2])&&(x[4] in [1, 2]),list))
    elseif BFtype=="F+B-"
        return (list-> filter(x-> (x[3] in [3, 4])&&(x[4] in [1, 2]),list))
    elseif BFtype=="B+F-"
        return (list-> filter(x-> (x[3] in [1, 2])&&(x[4] in [3, 4]),list))
    elseif BFtype=="F+BF-"
        return (list-> filter(x-> ((x[3] in [3, 4])&&(x[4] in [3, 4]))||((x[3] in [3, 4])&&(x[4] in [1, 2])),list))
    elseif BFtype=="BFFB"
        return (list-> filter(x-> ((x[3] in [1, 2])&&(x[4] in [3, 4]))||((x[3] in [3, 4])&&(x[4] in [1, 2])),list))
    else return identity
    end
end

function filterBFwithString(BFtype::String)
    if BFtype=="F+F-"
        return (list-> filter(x-> (x[1][3] in [3, 4])&&(x[1][4] in [3, 4]),list))
    elseif BFtype=="B+B-"
        return (list-> filter(x-> (x[1][3] in [1, 2])&&(x[1][4] in [1, 2]),list))
    elseif BFtype=="F+B-"
        return (list-> filter(x-> (x[1][3] in [3, 4])&&(x[1][4] in [1, 2]),list))
    elseif BFtype=="B+F-"
        return (list-> filter(x-> (x[1][3] in [1, 2])&&(x[1][4] in [3, 4]),list))
    elseif BFtype=="F+BF-"
        return (list-> filter(x-> ((x[1][3] in [3, 4])&&(x[1][4] in [3, 4]))||((x[1][3] in [3, 4])&&(x[1][4] in [1, 2])),list))
    elseif BFtype=="BFFB"
        return (list-> filter(x-> ((x[1][3] in [1, 2])&&(x[1][4] in [3, 4]))||((x[1][3] in [3, 4])&&(x[1][4] in [1, 2])),list))
    else return identity
    end
end


function functionalpairs(lambda::Int, type::String="F+F-") #First two index for α β+ type, the second for α β-
    funcind=lambdapairs(lambda)
    bfind=[[i,j] for i=1:4 for j=1:4]
    tmpind=[vcat(funcind[i], bfind[j]) for i in eachindex(funcind) for j in eachindex(bfind)]
    allind=sort(filter(x-> !((x[3]==1&&x[1]==1)||(x[4]==1&&x[2]==1)), tmpind), by=x->(x[1]+x[2], x)) 
    return filterBF(type)(allind)
end


function functionalstringpairs(lambda::Int, type::String="F+F-") 
    #First two index for subscript, the second for + - type , in the order of βB αB βF αF
    funcind=lambdapairs(lambda)
    bfind=[[i,j] for i=1:4 for j=1:4]
    tmpind=[(vcat(funcind[i], bfind[j]), functionaltype(bfind[j]...)[1]*string_subscript(funcind[i][1]-1)*" "*functionaltype(bfind[j]...)[2]*string_subscript(funcind[i][2]-1)) for i in eachindex(funcind) for j in eachindex(bfind)]
#     filter(x-> !((x[1][3]==1&&x[1][1]==1)||(x[1][4]==1&&x[1][2]==1)), tmpind)
    allind=sort(filter(x-> !((x[1][3]==1&&x[1][1]==1)||(x[1][4]==1&&x[1][2]==1)), tmpind), by=x->(x[1][1]+x[1][2], x[1]))
    return map(x->x[2], filterBFwithString(type)(allind))
end

function spindtable2d(l::Int, dtableinstance::DeltaTable{T}) where T <: Real
    return (2 .* dtableinstance.dtable .+ l)[1:end-l*dtableinstance.grid]
end

function spinlength2d(l::Int, dtableinstance::DeltaTable{T}) where T <: Real
    dtableinstance.tablelength- l*dtableinstance.grid
end

# function spindtable3d(l::Int, dtableinstance::DeltaTable{T}, truncN::Int) where T <: Real
#     return spindtable2d(l, dtableinstance)[1 + div(dtableinstance.grid+1, 2): end-truncN*dtableinstance.grid]
# end

# function spinlength3d(l::Int, dtableinstance::DeltaTable{T}, truncN::Int) where T <: Real
#     dtableinstance.tablelength- l*dtableinstance.grid -truncN*dtableinstance.grid - div(dtableinstance.grid+1, 2)
# end

# function nlnorm(l::Int, dphi::T) where T

#     l==0 ? 1. : (4.)^(-l) * l^(1.5 +2. *dphi)
# end

function nlnorm(l::Int, type::Type)

    l==0 ? type(1.) : type((4. ^ (-l ) ))
end

function nlnorm(l::Int, dp::T) where T

    l==0 ? T(1.) : T((4. ^ (-l ) )* l^(2*dp+1.5))
end

function spinfunctional2d(l::Int, lambda::Int, BFtype::String, functionaltable::OneDimFunctionalTable{T}) where T <: Real
   
    if lambda+1>functionaltable.order 
        error("lambda too high for the table")
    end

    # if functionaltable.type!="F" 
    #     error("Non fermionic functional not defined")
    # end

    if T<: BigFloat
        println(l)
    end

    indlist=functionalpairs(lambda, BFtype)
    shift=l*functionaltable.dtable.grid;
    farray=functionaltable.functionalarray
    nm=functionaltable.normminus
    np=functionaltable.normplus
    dp=functionaltable.dphi

    spinfunctionalmatrix=Array{T}(undef, spinlength2d(l, functionaltable.dtable), length(indlist))

    for i=eachindex(indlist)
        spinfunctionalmatrix[:,i].=(T(0.5) *nlnorm(l , dp)) .* (nm[1:end-shift] ./ nm[shift+1:end] .* farray[1:end-shift, indlist[i][1], indlist[i][3], 1].*farray[shift+1:end, indlist[i][2], indlist[i][4], 2] .+ np[1:end-shift] ./ np[shift+1:end] .*farray[shift+1:end, indlist[i][1], indlist[i][3], 1].*farray[1:end-shift, indlist[i][2], indlist[i][4], 2] )
    end

    return spinfunctionalmatrix
end

function spingrid2d(l::Int, dtableinstance::DeltaTable{T}) where T <: Real
    [(deltal, l) for deltal in spindtable2d(l, dtableinstance)]
end


function Anorm(dphi, Δ, l, j, n)
    log_result = -n*log(16) +
        lgamma(0.5 + n) - lgamma(0.5) +
        lgamma(-1 + Δ + 2*n) - lgamma(-1 + Δ) +
        lgamma(0.5*(-j + Δ) + n) - lgamma(0.5*(-j + Δ)) +
        lgamma(0.5*(j + Δ) + n) - lgamma(0.5*(j + Δ)) +
        lgamma(0.5*(-1 - l + Δ) + n) - lgamma(0.5*(-1 - l + Δ)) +
        lgamma(0.5*(l + Δ) + n) - lgamma(0.5*(l + Δ)) -
        lgamma(n + 1) -
        lgamma(-0.5 + Δ + n) + lgamma(-0.5 + Δ) -
        lgamma(0.5*(-1 - j + Δ) + n) + lgamma(0.5*(-1 - j + Δ)) -
        lgamma(0.5*(-1 + j + Δ) + n) + lgamma(0.5*(-1 + j + Δ)) -
        lgamma(0.5*(-l + Δ) + n) + lgamma(0.5*(-l + Δ)) -
        lgamma(0.5*(1 + l + Δ) + n) + lgamma(0.5*(1 + l + Δ)) -
        lgamma(-1 + n + Δ + n) + lgamma(-1 + n + Δ)+log(2) + lgamma(0.5 * (1 - j + l)) + lgamma(0.5 *  (1 + j + l)) - 
             (lgamma(0.5 *  (2 - j + l)) + lgamma(0.5 *  (2 + j + l)))-log(pi)+2*log(Δ - l ) + 14*log(2 - l + Δ) +
        2*lgamma(5 - dphi - j/2 + n + Δ/2) +
        4*lgamma(1/2*(-l + Δ)) +
        2*lgamma(-j + 2*n + Δ) +
        2*lgamma(dphi + 1/2*(-l + Δ)) -
        2*log(-j + 2*n + Δ) - 14*log(2 - j + 2*n + Δ) -
        2*lgamma(dphi - j/2 + n + Δ/2) -
        2*lgamma(-l + Δ) -
        2*lgamma(1/2*(10 - 2*dphi - l + Δ)) -
        4*lgamma(1/2*(-j + 2*n + Δ))

    # get the actual result
    result = exp(log_result)
    if j==0
        result=result/2
    end
    
    return result*nlnorm(l, dphi)/nlnorm(j, dphi)
end




function shifted3dmatrix(l::Int, j::Int, n::Int, truncN::Int, table2d::HigherDimFunctionalTable{T}) where T
    
    if l==0
        return table2d[j][2 + div(l-j, 2)*table2d.grid + div(table2d.grid+1, 4)+n*table2d.grid: end-(truncN-n)*table2d.grid -  div(l-j, 2)*table2d.grid, :]
    end

    return table2d[j][1 + div(l-j, 2)*table2d.grid + div(table2d.grid+1, 2)+n*table2d.grid: end-(truncN-n)*table2d.grid -  div(l-j, 2)*table2d.grid, :]
end

# function shifted3ddtable(l::Int, j::Int, n::Int, truncN::Int, table2d::HigherDimFunctionalTable{T}) where T
#     return table2d.dtablelist[div(j, 2)+1][1 + div(l-j, 2)*table2d.grid + div(table2d.grid+1, 2)+n*table2d.grid: end-(truncN-n)*table2d.grid -  div(l-j, 2)*table2d.grid, :]|>vec
# end

function shifted3ddtable(l::Int, j::Int, n::Int, truncN::Int, table2d::HigherDimFunctionalTable{T}) where T
    if l==0
        return table2d.dtablelist[div(j, 2)+1][2 + div(l-j, 2)*table2d.grid + div(table2d.grid+1, 4)+n*table2d.grid: end-(truncN-n)*table2d.grid -  div(l-j, 2)*table2d.grid, :]|>vec
    end 
    return table2d.dtablelist[div(j, 2)+1][1 + div(l-j, 2)*table2d.grid + div(table2d.grid+1, 2)+n*table2d.grid: end-(truncN-n)*table2d.grid -  div(l-j, 2)*table2d.grid, :]|>vec
end

# function spindtable3d(l::Int, truncN::Int, dtableinstance::DeltaTable{T}) where T <: Real
#     return spindtable2d(l, dtableinstance)[1 + div(dtableinstance.grid+1, 2): end-truncN*dtableinstance.grid]
# end

function spindtable3d(l::Int, dtableinstance::DeltaTable{T}, truncN::Int) where T <: Real
    if l==0
        return spindtable2d(l, dtableinstance)[2 + div(dtableinstance.grid+1, 4): end-truncN*dtableinstance.grid]
    end
    return spindtable2d(l, dtableinstance)[1 + div(dtableinstance.grid+1, 2): end-truncN*dtableinstance.grid]
end

# function spinlength3d(l::Int, truncN::Int, dtableinstance::DeltaTable{T}) where T <: Real
#     dtableinstance.tablelength- l*dtableinstance.grid -truncN*dtableinstance.grid - div(dtableinstance.grid+1, 2)
# end

function spinlength3d(l::Int, dtableinstance::DeltaTable{T}, truncN::Int) where T <: Real
    if l==0
        return dtableinstance.tablelength- l*dtableinstance.grid -truncN*dtableinstance.grid - div(dtableinstance.grid+1, 4)-1
    end
    return dtableinstance.tablelength- l*dtableinstance.grid -truncN*dtableinstance.grid - div(dtableinstance.grid+1, 2)
end

# function spingrid3d(l::Int, truncN::Int, dtableinstance::DeltaTable{T}) where T <: Real
#     [(deltal, l) for deltal in spindtable3d(l, dtableinstance, truncN)]
# end

function spingrid3d(l::Int, truncN::Int, dtableinstance::DeltaTable{T}) where T <: Real
    [(deltal, l) for deltal in spindtable3d(l, dtableinstance, truncN)]
end


function spinjcontribution(l::Int, j::Int, truncN::Int, table2d::HigherDimFunctionalTable{T}) where T
    normlist=Vector{Vector{T}}(undef, truncN+1)
    Threads.@threads for i=0:truncN
        tempfunc=(x->Anorm(table2d.dphi, x, l, j, i))
        normlist[i+1]=tempfunc.(shifted3ddtable(l, j, 0, truncN, table2d))
    end
    normlist[1][1]=normlist[1][2]
    spinj=shifted3dmatrix(l, j, 0, truncN, table2d).*normlist[1]
    for i=1:truncN
        spinj+=shifted3dmatrix(l, j, i, truncN, table2d).*normlist[i+1]
    end
    return spinj
end

function spinfunctional3d(l::Int, lambda::Int, truncN::Int, table2d::HigherDimFunctionalTable{T}) where T <: Real
    
    indlist=functionalpairs(lambda, table2d.type)
    spinfunctionalmatrix=spinjcontribution(l, 0, truncN, table2d)
    for j=2:2:l
        spinfunctionalmatrix+=spinjcontribution(l, j, truncN, table2d)
    end
    return spinfunctionalmatrix[:, 1:length(indlist)]
end

function generateHigherDimFunctionalTable(llist::Vector{Int}, lambda::Int, truncN::Int, table2d::HigherDimFunctionalTable{T}, functionaltable::OneDimFunctionalTable{T} ) where T <: Real
    if lambda>table2d.order 
        error("lambda too high for the table")
    end

    # if functionaltable.type!="F" 
    #     error("Non fermionic functional not defined")
    # end

    spinnum=length(llist)
    dtablelist=Vector{Vector{T}}(undef, spinnum)
    gridlist=Vector{Vector{Tuple{T, Int64}}}(undef, spinnum)
    functionalarray=Vector{Matrix{T}}(undef, spinnum)

    # Threads.@threads for (i, l) in enumerate(llist)
    #     dtablelist[i]=spindtable2d(l, functionaltable.dtable)
    #     gridlist[i]=spingrid2d(l, functionaltable.dtable)
    #     functionalarray[i]=spinfunctional2d(l, lambda, functionaltable)
    # end

    for i in 1:length(llist)
        println(i)
        l=llist[i]
        dtablelist[i]=spindtable3d(l, functionaltable.dtable, truncN)
        gridlist[i]=spingrid3d(l, truncN, functionaltable.dtable)
        functionalarray[i]=spinfunctional3d(l, lambda, truncN, table2d)
    end

    return HigherDimFunctionalTable(3, table2d.dphi, llist, dtablelist, length.(dtablelist), gridlist, lambda, Dict(zip(llist, functionalarray)), -(table2d.functionalDic[0][1,:]), length(functionalpairs(lambda, table2d.type)), functionaltable.grid, table2d.type, nothing, functionalstringpairs(lambda, table2d.type))
end


# function pochhammer(x, n)
#     gamma(x + n) / gamma(x)
# end


# function Asansc(delta::T, l::Int, j::Int, n::Int ) where T
#     (j*gamma(l+1)*pochhammer(0, j)*pochhammer(0.5, (-j + l)/2.)*pochhammer(0.5,j + (-j + l)/2.)*pochhammer(0.5,n)*pochhammer(-1 + delta,2*n)*pochhammer((delta - j)/2.,n)*pochhammer((delta + j)/2.,n)*pochhammer((-1 + delta - l)/2.,n)*pochhammer((delta - l)/2.,n)*pochhammer((delta + l)/2.,n)*pochhammer((1 + delta + l)/2.,n))/(Power(16,n)*Factorial(j)*Factorial((-j + l)/2.)*Factorial(n)*pochhammer(0,1 + j + (-j + l)/2.)*pochhammer(1,l)*pochhammer(-0.5 + delta,n)*pochhammer((-1 + delta - j)/2.,n)*pochhammer((-1 + delta + j)/2.,n)*pochhammer(-1 + delta + n,n)) 
# end

# function Anorm(delta::T, l::Int, j::Int, n::Int ) where T
#     ((8 (2 + delta - l)^14 (-1 + delta - j + 2 n))/((delta - l) (delta - j + 2 n)^2 (2 + delta - j + 2 n)^14 pi^(3/2) )*(gamma(-0.5 + delta)*gamma((-1 + delta - j)/2.)*gamma((-1 + delta + j)/2.)*Power(gamma((2 + delta - l)/2.),3)*Power(gamma(delta/2. + dphi - l/2.),2)*gamma((1 + delta + l)/2.)*gamma((1 - j + l)/2.)*gamma((1 + j + l)/2.)*gamma(0.5 + n)*gamma(-1 + delta + n)*Power(gamma(5 + delta/2. - dphi - j/2. + n),2)*gamma((1 + delta - j + 2*n)/2.)*gamma((delta + j + 2*n)/2.)*gamma((-1 + delta - l + 2*n)/2.)*gamma((delta + l + 2*n)/2.))/(gamma(-1 + delta)*gamma((delta - j)/2.)*gamma((delta + j)/2.)*gamma((-1 + delta - l)/2.)*Power(gamma((1 + delta - l)/2.),2)*Power(gamma((10 + delta - 2*dphi - l)/2.),2)*gamma((delta + l)/2.)*gamma((2 - j + l)/2.)*gamma((2 + j + l)/2.)*gamma(1 + n)*gamma(-0.5 + delta + n)*gamma(delta/2. - j/2. + n)*Power(gamma(delta/2. + dphi - j/2. + n),2)*gamma(delta/2. - l/2. + n)*gamma((-1 + delta + j + 2*n)/2.)*gamma((1 + delta + l + 2*n)/2.))
# end




"""
    generateHigherDimFunctionalTable(llist::Vector{Int}, lambda::Int, functionaltable::OneDimFunctionalTable{T}, dim::Int=2)
generateHigherDimFunctionalTable(llist::Vector{Int64}, lambda::Int64, functionaltable::OneDimFunctionalTable{T}) where T<:Real
"""
function generateHigherDimFunctionalTable(llist::Vector{Int}, lambda::Int, BFtype::String, functionaltable::OneDimFunctionalTable{T}, dim::Int=2) where T <: Real
    @assert dim==2

    if lambda+1>functionaltable.order 
        error("lambda too high for the table")
    end

    # if functionaltable.type!="F" 
    #     error("Non fermionic functional not defined")
    # end

    spinnum=length(llist)
    dtablelist=Vector{Vector{T}}(undef, spinnum)
    gridlist=Vector{Vector{Tuple{T, Int64}}}(undef, spinnum)
    functionalarray=Vector{Matrix{T}}(undef, spinnum)

    # Threads.@threads for (i, l) in enumerate(llist)
    #     dtablelist[i]=spindtable2d(l, functionaltable.dtable)
    #     gridlist[i]=spingrid2d(l, functionaltable.dtable)
    #     functionalarray[i]=spinfunctional2d(l, lambda, functionaltable)
    # end

    Threads.@threads for i in 1:length(llist)
        l=llist[i]
        dtablelist[i]=spindtable2d(l, functionaltable.dtable)
        gridlist[i]=spingrid2d(l, functionaltable.dtable)
        functionalarray[i]=spinfunctional2d(l, lambda, BFtype, functionaltable)
    end

    return HigherDimFunctionalTable(dim, functionaltable.dphi, llist, dtablelist, length.(dtablelist), gridlist, lambda, Dict(zip(llist, functionalarray)), -(functionalarray[1][1, :]), length(functionalpairs(lambda, BFtype)), functionaltable.grid, BFtype, nothing, functionalstringpairs(lambda, BFtype))
end

