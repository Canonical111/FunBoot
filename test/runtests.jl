using FunBoot
#using AppleAccelerate
# if Intel CPU, using MKL
# if AMD CPU, delete it
singlepointfilename="2dIsing.h5"
newdtable, newfunctionaltable=readsinglepoint(singlepointfilename);
println("Test with 32 FB functional at 2d Ising");
newfunctionaltable2d=generateHigherDimFunctionalTable(vcat(collect(0:2:60),[120]), 3, "F+B-", newfunctionaltable);
newfunctionaltable2d.itp=generateintepolation(newfunctionaltable2d, 10);
#newfunctionaltable2d.functionallist
newLP=initializeLP(newfunctionaltable2d, vcat(collect(0:2:60),[120]), 3, 100; spintoGap=0);
@elapsed iterateCompact!(newLP, 9000)
for i=1:10
    outer_approx(newLP, 20, newfunctionaltable2d);
end
solutionprint(newLP, newfunctionaltable)
