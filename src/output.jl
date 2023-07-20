

function string_subscript(n::Int)
    # Define a string with all subscript digits
    subscripts = "₀₁₂₃₄₅₆₇₈₉"
    
    # Convert the integer to a string, replace each digit with its subscript
    # equivalent, and print the result
    last(first((subscripts),n+1))
end

function functionaltype(plustype::Int, minustype::Int)
    # Define a string with all subscript digits
    functypestring = [["β+B", "α+B", "β+F", "α+F"], ["β-B", "α-B", "β-F", "α-F"]]
    # Convert the integer to a string, replace each digit with its subscript
    # equivalent, and print the result
    return (functypestring[1][plustype], functypestring[2][minustype])
end

function plotfunctional(l, lp, hightable)
    if l==lp.gapchannel
        x=hightable.dtablelist[div(l, 2)+1][1:2*hightable.grid]
        y=((x->costf(l, x, lp.gapchannel)).(hightable.dtablelist[div(l, 2)+1])-hightable[l]*lp.functional)             
        plot(x, log.(y[1:2*hightable.grid]))
    else
        x=hightable.dtablelist[div(l, 2)+1][1:2*hightable.grid]
        y=(-hightable[l]*lp.functional)
        plot(x, (y[1:2*hightable.grid]))
    end
end

function extractsol(solvec::Array{Tuple{Tuple{T, Int}, T},1}; cutoff= .2) where T
    test3=sort(solvec);
    protected=[];
    protectlist=Int64[];
    for (i,solvec1) in enumerate(test3)
        if solvec1[1][1]-floor(solvec1[1][1])<1e-40
            push!(protected,solvec1)
            push!(protectlist,i)
        end
    end
    for i in reverse(protectlist)
        deleteat!(test3,i)
    end

    temps=test3;
    ll=length(temps);
    op=[];
    for i=1:ll
        if i>1&&abs(temps[i][1][1]-temps[i-1][1][1])<=cutoff&&temps[i][1][2]==temps[i-1][1][2]
            continue
        end
        if i<ll&&abs(temps[i][1][1]-temps[i+1][1][1])<=cutoff&&temps[i][1][2]==temps[i+1][1][2]
            push!(op,(((temps[i][1][1]*temps[i][2]+temps[i+1][1][1]*temps[i+1][2])/(temps[i+1][2]+temps[i][2]),temps[i+1][1][2]),temps[i+1][2]+temps[i][2]))
        end
        if i==1&&(temps[i][1][2]!=temps[i+1][1][2]||abs(temps[i][1][1]-temps[i+1][1][1])>cutoff)
            push!(op,temps[1])
        end
        if i==ll
            push!(op,temps[i])
        end
        if i>1&&i<ll&&(abs(temps[i][1][1]-temps[i+1][1][1])>cutoff||temps[i][1][2]!=temps[i+1][1][2])
            push!(op,temps[i])
        end
    end
sort([op;protected])
end


function solutionprint(lp::LinearProgramCompact{T}, onetable::OneDimFunctionalTable{T}) where T
    
    sortedoplist=sort([(lp.gridsol[i], lp.coeffs[i]) for i=1: lp.num_of_functional],by=x->x[1][1])
    sortedoplist=extractsol(sortedoplist)
    sortedoplist=map(x->[x[1][1],x[1][2],x[2]* onetable.normplusitp(x[1][1]-x[1][2])* onetable.normminusitp(x[1][1]-x[1][2])*nlnorm(x[1][2], onetable.dphi)], sortedoplist)
    
    @printf("The solution with %d operators.\n", length(sortedoplist))
    println("|  L  |         Δ           |        OPE^2        |")
    println("|:---:|:-------------------:|:-------------------:|") 
    map(x-> @printf("| %3d | %18.14f | %18.14e |\n", x[2], x[1], x[3]),sortedoplist);
    println("\n")
end

function solution(lp::LinearProgramCompact{T}, onetable::OneDimFunctionalTable{T}) where T
    
    sortedoplist=sort([(lp.gridsol[i], lp.coeffs[i]) for i=1: lp.num_of_functional],by=x->x[1][1])
    sortedoplist=extractsol(sortedoplist)
    sortedoplist=map(x->((x[1][1],x[1][2]),x[2]* onetable.normplusitp(x[1][1]-x[1][2])* onetable.normminusitp(x[1][1]-x[1][2])*nlnorm(x[1][2], onetable.dphi) ), sortedoplist)
    return sortedoplist
end