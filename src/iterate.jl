using Printf

function addSparseTo!(v2::Vector{T}, v1::SparseVector{T, Int64},) where {T<: Real}
    for (i, ind) in enumerate(v1.nzind)
        v2[ind]+= v1.nzval[i]
    end
    return v2
end

function swap_columns!(m1, m2, ind1::Int, ind2::Int)
    @inbounds @simd for k=1:size(m1, 1)
        m1[k, ind1], m2[k, ind2] = m2[k, ind2], m1[k, ind1]
    end
end

function updateInverse!(lp::LinearProgramCompact{T}) where T
    copyto!(lp.invA.factors, lp.AMatrix)
    lp.invA=lu!(lp.invA.factors)
end

function updateFunctional!(lp::LinearProgramCompact{T}) where T
    ldiv!(lp.functional, (lp.invA)', lp.solcostlist)
end

function updateCoeffs!(lp::LinearProgramCompact{T}) where T
    ldiv!(lp.coeffs, lp.invA, lp.target)
end

function updateRClist!(lp::LinearProgramCompact{T}) where T
    addSparseTo!(mul!(lp.rclist, lp.vecMatrix', (-lp.functional)), lp.veccostlist)
end

function updateneginfo!(lp::LinearProgramCompact{T}) where T
    lp.negindlist=findall(<(0),lp.rclist)
    lp.negnum=length(lp.negindlist)
    lp.negrclist =lp.rclist[lp.negindlist]
end

function updatecostvar!(lp::LinearProgramCompact{T}) where T
    icolM=ldiv!(lp.invA, (lp.vecMatrix[:,lp.negindlist]))
    lp.costvar = Vector{T}(undef, lp.negnum)
    lp.bbaroutindlist = similar(lp.negindlist)
    @inbounds @simd for i in 1:lp.negnum
            minvals=convert(T, 1e256)
            minind=1
            for j in 1:lp.num_of_functional
                minvals1 = (icolM[j,i] > zero(T) ? lp.coeffs[j] / icolM[j,i] : convert(T, 1e256))
                if minvals1 < minvals
                    minind=j
                    minvals=minvals1
                end
        end
        lp.costvar[i]=minvals
        lp.bbaroutindlist[i]=minind
    end
    lp.costvar.=lp.negrclist.*lp.costvar
end

function swapbbar!(lp::LinearProgramCompact{T}) where T
    lp.costvarmax, negindbbar = findmin(lp.costvar)
    bbarinind=lp.negindlist[negindbbar]
    bbaroutind=lp.bbaroutindlist[negindbbar]
    @printf("%13.5e |", -lp.costvarmax)
   
    lp.bbarin=(bbarinind, lp.gridVectors[bbarinind])
    lp.bbarout=(bbaroutind, lp.gridsol[bbaroutind])
    @printf("%9.3f |", lp.bbarin[2][1])
    lp.bbarin[2][2]>=0 ? @printf("%3d |", lp.bbarin[2][2]) : @printf("AUX |")
    @printf("%9.3f |", lp.bbarout[2][1])
    lp.bbarout[2][2]>=0 ? @printf("%3d |", lp.bbarout[2][2]) : @printf("AUX |")
    swap_columns!(lp.vecMatrix, lp.AMatrix, bbarinind, bbaroutind)
    lp.gridVectors[bbarinind], lp.gridsol[bbaroutind] = lp.gridsol[bbaroutind], lp.gridVectors[bbarinind]
    lp.veccostlist[bbarinind], lp.solcostlist[bbaroutind] = lp.solcostlist[bbaroutind], lp.veccostlist[bbarinind]
end

function updatecost!(lp::LinearProgramCompact{T}) where T
    lp.cost=dot(lp.coeffs, lp.solcostlist)
end

function reverseswap!(lp::LinearProgramCompact{T}) where T
    bbarinind=lp.bbarin[1]
    bbaroutind=lp.bbarout[1]
    swap_columns!(lp.vecMatrix, lp.AMatrix, bbarinind, bbaroutind)
    lp.gridVectors[bbarinind], lp.gridsol[bbaroutind] = lp.gridsol[bbaroutind], lp.gridVectors[bbarinind]
    lp.veccostlist[bbarinind], lp.solcostlist[bbaroutind] = lp.solcostlist[bbaroutind], lp.veccostlist[bbarinind]
    updateInverse!(lp)
    updateCoeffs!(lp)
    if lp.bbarin[2][2]<0  
        lp.AUXnum -= 1
    end

    if lp.bbarout[2][2]<0 
         lp.AUXnum += 1
    end
    
    updatecost!(lp)
    updateFunctional!(lp) 
end

function checkstatus!(lp::LinearProgramCompact{T}, checknum::Int) where T
    err=abs((lp.costprevious-lp.cost + lp.costvarmax)/lp.costvarmax);
    # @printf("%2.2e |", err) 
    # checknum==1 ? @printf("%d |\n", lp.err) : @printf("%d |", lp.err)
    
    while (err>1e-1||lp.err>0)&&checknum<=lp.negnum-1
        reverseswap!(lp)
        
        negindbbar = partialsortperm(lp.costvar, 1+checknum )
        lp.costvarmax = lp.costvar[negindbbar]
        bbarinind=lp.negindlist[negindbbar]
        bbaroutind=lp.bbaroutindlist[negindbbar]
        # @printf("%13.5e |", -lp.costvarmax)
        lp.bbarin=(bbarinind, lp.gridVectors[bbarinind])
        lp.bbarout=(bbaroutind, lp.gridsol[bbaroutind])
        # @printf("%9.3f |", lp.bbarin[2][1])
        # lp.bbarin[2][2]>=0 ? @printf("%3d |", lp.bbarin[2][2]) : @printf("AUX |")
        # @printf("%9.3f |", lp.bbarout[2][1])
        # lp.bbarout[2][2]>=0 ? @printf("%3d |", lp.bbarout[2][2]) : @printf("AUX |")
        swap_columns!(lp.vecMatrix, lp.AMatrix, bbarinind, bbaroutind)
        lp.gridVectors[bbarinind], lp.gridsol[bbaroutind] = lp.gridsol[bbaroutind], lp.gridVectors[bbarinind]
        lp.veccostlist[bbarinind], lp.solcostlist[bbaroutind] = lp.solcostlist[bbaroutind], lp.veccostlist[bbarinind]    
        updateInverse!(lp)
        updateCoeffs!(lp)
        if lp.bbarin[2][2]<0  
            lp.AUXnum += 1
        end

        if lp.bbarout[2][2]<0 
             lp.AUXnum -= 1
        end

        if lp.bbarout[2][2]<0 && lp.AUXnum==0
            lp.status="Feasible"
        end

        lp.costprevious =lp.cost;
        updatecost!(lp)
        lp.err=count(<(0), lp.coeffs)
        err=abs((lp.costprevious-lp.cost + lp.costvarmax)/lp.costvarmax) 
        checknum+=1
    end

    if (err>1e-1||lp.err>0)
        checknum+=1
    end

    return checknum
end

function iterateCompact!(lp::LinearProgramCompact{T}, n::Int; init=true, quiet=false) where T
    starttime=time()
    println("Started at: $(Libc.strftime(time()))") 
    println("| Ite | negnum |   |costvar|  | sol_in_d | l  | sol_out_d| l  |     cost     | AUX |  time  |err| relerr|")
    println("|:---:|:------:|:------------:|:--------:|:--:|:--------:|:--:|:------------:|:---:|:------:|:-:|:-----:|") 
    for i=1:n
        ittime=time()
        @printf("|%4d |", i)
        updateFunctional!(lp)
        updateRClist!(lp)
        updateneginfo!(lp)
        @printf("%7d |", lp.negnum)
        quiet=false
        if lp.negnum==0
            if !quiet println("\n  Min cost achieved") end
            lp.status="Minimized"
             break
        end
        
        updatecostvar!(lp)
        swapbbar!(lp)
        updateInverse!(lp)
        updateCoeffs!(lp)
        if lp.bbarin[2][2]<0  
            lp.AUXnum += 1
        end

        if lp.bbarout[2][2]<0 
             lp.AUXnum -= 1
        end

        if lp.bbarout[2][2]<0 && lp.AUXnum==0
            lp.status="Feasible"
        end

        lp.costprevious =lp.cost;
        updatecost!(lp);
        lp.err=count(<(0), lp.coeffs);
        
        
        if i>1||!init
            if checkstatus!(lp, 1)>lp.negnum
                reverseswap!(lp)
                println("\nNeeds more precision to proceed.")
                break
            end
        end
        @printf("%13.5e |", lp.cost)
        @printf("%4d |", lp.AUXnum)
        @printf("%2.1e |", time()-ittime) 
        err=abs((lp.costprevious-lp.cost + lp.costvarmax)/lp.costvarmax);
        @printf("%2d | ", lp.err)
        @printf("%4.0e |\n", err) 
        
    end
    updateFunctional!(lp)
    quiet ? solution(lp) : nothing
end

