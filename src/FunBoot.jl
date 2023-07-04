module FunBoot

using SparseArrays, HDF5, LinearAlgebra, LoopVectorization, Printf, BSplineKit, PyPlot

export addSparseTo!, swap_columns!, lambdapairs, griddtable, OneDimFunctionalTable, DeltaTable, spindtable2d, spinlength2d, spingrid2d, spinfunctional2d, readsinglepoint

export iterateCompact!, status, functionalpairs, generateHigherDimFunctionalTable, filteroutop!, filterfunctional!, solution, solutionRaw

export generateHigherDimFunctionalTable, HigherDimFunctionalTable, getindex, convolveiterator, initializeIndices, initializeLP, LinearProgramCompact 

export updateFunctional!,  updateRClist!, updateneginfo!, updatecostvar!, swapbbar!, updateInverse!, updateCoeffs!, updatecost!, addonepoint!, functionalop

export str_subscript, filterBFind, filterBFwithString, generateintepolation, costf, outer_approx,  plotfunctional, extractsol, solutionprint, convert_type


include("OneDimTable.jl")

include("HigherDimTable.jl")

include("LP.jl")

include("iterate.jl")

include("output.jl")

include("outer.jl")

##################################################
end

