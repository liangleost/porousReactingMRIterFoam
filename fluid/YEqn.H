
tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

//if(solveSpecies)
{
    label inertIndex = -1;
    volScalarField Yt = 0.0*Y[0];
    forAll(Y, j)
    {
        if (Y[j].name() != inertSpecie)
        {
            volScalarField& Yj = Y[j];
            volScalarField YRESJ = Y[j];
            solve
            (
                fvm::ddt(rho, Yj)
              + mvConvection->fvmDiv(phi, Yj)
              - fvm::laplacian(turb.muEff(), Yj)
              ==
                parcels.Srho(j)
               + kappa*chemistry.RR(j)().dimensionedInternalField()
              + massSource.Su(j),
                mesh.solver("Yj")
            );
            Yj.max(0.0);
            Yt += Yj;
            specieres=max(mag(Yj-YRESJ)/max(Yj,dimensionedScalar("1e-10",dimensionSet(0,0,0,0,0,0,0),1e-15))).value();
            YRESJ.clear();
        }
        else
        {
            inertIndex = j;
        }
        spres=max(spres,specieres);
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
