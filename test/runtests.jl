using FunBoot
#using AppleAccelerate
# if Intel CPU, use mkl
# if AMD CPU, delete it
singlepointfilename="2dIsing.h5"
newdtable, newfunctionaltable=readsinglepoint(singlepointfilename);
newfunctionaltable2d=generateHigherDimFunctionalTable(vcat(collect(0:2:120),[256]), 4, "F+B-", newfunctionaltable);
newfunctionaltable2d.itp=generateintepolation(newfunctionaltable2d, 10);
#newfunctionaltable2d.functionallist
newLP=initializeLP(newfunctionaltable2d, vcat(collect(0:2:120),[256]), 4, 50; spintoGap=0);
@elapsed iterateCompact!(newLP, 9000)
for i=1:10
    outer_approx(newLP, 20, newfunctionaltable2d);
end
solutionprint(newLP, newfunctionaltable)
