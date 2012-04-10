
if (oCorr == 0)
{
    parcels.evolve();
    #include "chemistry.H"
    #include "rhoEqn.H"
}
#include "UEqn.H"
#include "YEqn.H"

// --- PISO loop
for (int corr=0; corr<nCorr; corr++)
{
    #include "hsEqn.H"
    #include "pEqn.H"
}
//turb.correct();

//rho = thermo.rho();

//if(runTime.write())
//{
//    chemistry.dQ()().write();
//}
