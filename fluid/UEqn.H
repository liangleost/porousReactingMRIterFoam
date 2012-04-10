    fvVectorMatrix UEqn
    (
//      pZones.ddt(rho,U)
      fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turb.divDevRhoReff(U)
     ==
        rho.dimensionedInternalField()*gFluid[i]
      + parcels.SU()
      + momentumSource.Su()
    );
    pZones.addResistance(UEqn);

    if (momentumPredictor)
    {
        uxres=solve(UEqn == -fvc::grad(p)).initialResidual();
    }

