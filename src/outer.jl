
function gridin(op::Tuple{T, Int64}, gridVectors::Vector{Tuple{T, Int64}}, solVectors::Vector{Tuple{T, Int64}}, gridVectors1::Vector{Tuple{T, Int64}}) where T
    gridlist=filter!(x->x>T(1e-8), map(x->abs(x[1]-op[1]), vcat(filter(x->x[2]==op[2], gridVectors), filter(x->x[2]==op[2], solVectors), filter(x->x[2]==op[2], gridVectors1)) ))
    if length(gridlist)==0
        return T(1e-9)
    else
        return findmin(gridlist)[1]
    end  
end

function newgrid(op::Tuple{T, Int64}, gridVectors::Vector{Tuple{T, Int64}}, solVectors::Vector{Tuple{T, Int64}}, gridVectors1::Vector{Tuple{T, Int64}}) where T
    gridop=gridin(op, gridVectors, solVectors, gridVectors1)
    if gridop<T(2e-9)
        return Vector{Tuple{T, Int64}}()
    end
    if op[2]+gridop/3>op[1]
        return [(op[1]+gridop/3, op[2])]
    elseif op[1]>T(250.)
        return Vector{Tuple{T, Int64}}()
    else
        return [(op[1]-gridop/3, op[2]), (op[1]+gridop/3, op[2])]
    end
end

function newgridlist(oplist::Vector{Tuple{T, Int64}}, gridVectors::Vector{Tuple{T, Int64}}, solVectors::Vector{Tuple{T, Int64}}, gridVectors1::Vector{Tuple{T, Int64}}) where T
    return vcat(map(x-> newgrid(x, gridVectors, solVectors, gridVectors1), oplist)...)
end

function functionalop(op::Tuple{T, Int64}, table2d::HigherDimFunctionalTable{T}) where T
    map(f->f(op[1]), table2d.itp[findfirst(isequal(op[2]), table2d.ltable)])
end

function outer_approx(lp::LinearProgramCompact{T}, iter::Int, table2d::HigherDimFunctionalTable{T}) where T
    lp_outer = copy(lp)
    for i=1:iter
        lp_outer.gridVectors = newgridlist(lp_outer.gridsol, lp.gridVectors, lp_outer.gridsol, lp_outer.gridVectors);
        lp_outer.vecMatrix=reduce(hcat, map(x->functionalop(x, table2d), lp_outer.gridVectors))[1:lp.num_of_functional, :]    
        lp_outer.veccostlist=sparse(map(x->costf(x[2], x[1], lp_outer.gapchannel), lp_outer.gridVectors));
        lp_outer.rclist= Vector{Float64}(undef, length(lp_outer.gridVectors))
        # updateInverse!(lp_outer)
        lp_outer.invA=lp.invA
        println("SmallTable")
        iterateCompact!(lp_outer, 1000; init=false)
        # updateInverse!(lp)
        lp.invA=lp_outer.invA
        println("BigTable")
        iterateCompact!(lp, 1000; init=false)
    end
end