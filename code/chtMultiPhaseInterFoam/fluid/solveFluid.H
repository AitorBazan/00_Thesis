if (finalIter)
{
    mesh.data::add("finalIteration", true);
}

#include "initCorrectPhi.H"

if (firstIter)
{
    thermol.correctMassSources(T);
    thermol.solve();
    rho = thermol.rho();

}

if (frozenFlow)
{
    #include "TEqnFluidPhase.H"
}
else
{
    
    #include "UEqnFluidPhase.H"
    #include "YEqnFluidPhase.H"      
    #include "TEqnFluidPhase.H"
        // --- PISO loop

    for (int corr=0; corr<nCorr; corr++)
    {
        #include "pEqnFluidPhase.H"
    }
    if (pimple.turbCorr())
    {
        turbulence.correct();
    }

}

if (finalIter)
{
    rho = thermol.rho();
    mesh.data::remove("finalIteration");
}
