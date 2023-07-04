using FunBoot
using AppleAccelerate
# if Intel CPU, use mkl
# if AMD CPU, delete it
singlepointfilename="phionedot2.h5"
newdtable, newfunctionaltable=readsinglepoint(singlepointfilename);
newfunctionaltable2d=generateHigherDimFunctionalTable(vcat(collect(0:2:80),[128]), 4, "F+F-", newfunctionaltable);
newfunctionaltable2d.itp=generateintepolation(newfunctionaltable2d, 8);
#newfunctionaltable2d.functionallist
newLP=initializeLP(newfunctionaltable2d, vcat(collect(0:2:80),[128]), 4, 100; spintoGap=0);
@elapsed iterateCompact!(newLP, 9000)
outer_approx(newLP, 20, newfunctionaltable2d);
outer_approx(newLP, 20, newfunctionaltable2d);
solutionprint(newLP, newfunctionaltable)
