    // Initialise fluid field pointer lists
    PtrList<IOdictionary> chemistryProperties(fluidRegions.size());
    PtrList<rhoChemistryModel> pChemistryFluid(fluidRegions.size());
    PtrList<hsReactionThermo> thermoFluid(fluidRegions.size());
    PtrList<DimensionedField<scalar, volMesh> > kappaFluid(fluidRegions.size());
    PtrList<volScalarField> rhoFluid(fluidRegions.size());
    PtrList<volScalarField> KFluid(fluidRegions.size());
    PtrList<volVectorField> UFluid(fluidRegions.size());
    PtrList<volScalarField> GFluid(fluidRegions.size());
    PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
    PtrList<DimensionedField<scalar, volMesh> > chemistryShFluid(fluidRegions.size());
    PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
    PtrList<compressible::turbulenceModel> turbulence(fluidRegions.size());
    PtrList<icoPoly8ThermoReactingMultiphaseCloud> parcelFluid(fluidRegions.size());
    PtrList<porousZones> porousFluid(fluidRegions.size());
    PtrList<radiation::radiationModel> radiationFluid(fluidRegions.size());
    PtrList<scalarTimeActivatedExplicitSourceList> massSourceFluid(fluidRegions.size());
    PtrList<vectorTimeActivatedExplicitSourceList> momentumSourceFluid(fluidRegions.size());
    PtrList<scalarTimeActivatedExplicitSourceList> energySourceFluid(fluidRegions.size());

    List<scalar> initialMassFluid(fluidRegions.size());

    // Populate fluid field pointer lists
    forAll(fluidRegions, i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidRegions[i].name() << nl << endl;

        Info<< "    Adding to thermoFluid\n" << endl;
        pChemistryFluid.set
        (
            i,
            autoPtr<rhoChemistryModel>
            (
                rhoChemistryModel::New
                (
                    fluidRegions[i]
                )
            ).ptr()
        );
        thermoFluid.set
        (
            i,
            hsReactionThermo::New(fluidRegions[i]).ptr()
        );

        Info<< "*** Reading chemistry properties for region "
            << fluidRegions[i].name() << nl <<endl;
        chemistryProperties.set
        (
            i,
            new IOdictionary
            (
                IOobject
                (
                    "chemistryProperties",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            )
        );
        Info<< "    Adding to rhoFluid\n" << endl;
        rhoFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                thermoFluid[i].rho()
            )
        );

        Info<< "    Adding to KFluid\n" << endl;

        KFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "K",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermoFluid[i].Cp()*thermoFluid[i].alpha()
            )
        );

        Info<< "    Adding to UFluid\n" << endl;
        UFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to phiFluid\n" << endl;
        phiFluid.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(rhoFluid[i]*UFluid[i])
                    & fluidRegions[i].Sf()
            )
        );

        Info<< "    Adding to gFluid\n" << endl;
        gFluid.set
        (
            i,
            new uniformDimensionedVectorField
            (
                IOobject
                (
                    "g",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
        Info<< "    Adding to turbulence\n" << endl;
        turbulence.set
        (
            i,
            autoPtr<compressible::turbulenceModel>
            (
                compressible::turbulenceModel::New
                (
                    rhoFluid[i],
                    UFluid[i],
                    phiFluid[i],
                    thermoFluid[i]
                )
            ).ptr()
        );
        radiationFluid.set
        (
            i,
            autoPtr<radiation::radiationModel>
            (
                radiation::radiationModel::New
                (
                    thermoFluid[i].T()
                )
            ).ptr()
        );
        GFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "G",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );
        parcelFluid.set
        (
            i,
            new icoPoly8ThermoReactingMultiphaseCloud
            (
                "reactingCloud1",
                rhoFluid[i],
                UFluid[i],
                gFluid[i],
                thermoFluid[i]
            )
        );

        kappaFluid.set
        (
            i,
            new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "kappa",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("zero",dimless,0.0)
            )
        );


    Info<< "Creating mass source\n" << endl;
    massSourceFluid.set
    (
        i,
        new scalarTimeActivatedExplicitSourceList
        (
            "mass",
            fluidRegions[i],
            dimMass/dimTime/dimVolume,
            thermoFluid[i].composition().species()
        )
    );


    Info<< "Creating momentum source\n" << endl;
    momentumSourceFluid.set
    (
        i,
        new vectorTimeActivatedExplicitSourceList
        (
            "momentum",
            fluidRegions[i],
            dimMass*dimVelocity/dimTime/dimVolume,
            "U"
        )
    );

    Info<< "Creating energy source\n" << endl;
    energySourceFluid.set
    (
        i,
        new scalarTimeActivatedExplicitSourceList
        (
            "energy",
            fluidRegions[i],
            dimEnergy/dimTime/dimVolume,
            "h"
        )
    );



        porousFluid.set
        (
            i,
            new porousZones(fluidRegions[i])
        );

        chemistryShFluid.set
        (
            i,
            new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "chemistry::Sh",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("chemistrySh",dimEnergy/dimTime/dimVolume, 0.0)
            )
        );

        initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

    }
