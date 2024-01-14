using FunBoot
#using AppleAccelerate
# if Intel CPU, use MKL
# if AMD CPU, delete it

singlepointfilename="functionalIsing12.h5"
newdtable, newfunctionaltable=readsinglepoint(singlepointfilename);

newfunctionaltable2d=generateHigherDimFunctionalTable(collect(0:2:100), 2, "F+F-", newfunctionaltable);

for i=eachindex(newfunctionaltable2d.dtablelist)
    newfunctionaltable2d.dtablelist[i].+=1e-15
end

newfunctionaltable3d=generateHigherDimFunctionalTable(collect(0:2:12), 1, 300, newfunctionaltable2d, newfunctionaltable);

newfunctionaltable3d.itp=generateintepolation(newfunctionaltable3d, 10);

newLP3d=initializeLP(newfunctionaltable3d, collect(0:2:12), 1, 200; spintoGap=0);

@elapsed iterateCompact!(newLP3d, 3000)

solutionRaw(newLP3d)

for i=1:10
    outer_approx(newLP3d, 20, newfunctionaltable3d)
end

solutionprint(newLP3d, newfunctionaltable)
