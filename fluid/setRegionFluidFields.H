    const fvMesh& mesh = fluidRegions[i];
    #include "readAdditionalSolutionControls.H"
    #include "createPorousZones.H"
    rhoChemistryModel& chemistry = pChemistryFluid[i];
    hsReactionThermo& thermo = chemistry.thermo();
    Switch turbulentReaction(chemistryProperties[i].lookup("turbulentReaction"));
    dimensionedScalar Cmix("Cmix",dimless,1.0);
    if(turbulentReaction)
    {
        chemistryProperties[i].lookup("Cmix") >> Cmix;
    }
    volScalarField& G = GFluid[i];
    radiation::radiationModel& radiation = radiationFluid[i];
    DimensionedField<scalar,volMesh>& kappa = kappaFluid[i];
    icoPoly8ThermoReactingMultiphaseCloud& parcels = parcelFluid[i];
    basicMultiComponentMixture& composition = thermo.composition();
    word inertSpecie(thermo.lookup("inertSpecie"));
    if(!composition.contains(inertSpecie))
    {
        FatalErrorIn(args.executable())
            << "Specified inert specie '" << inertSpecie << "' not found in "
            << "species list. Available species:" << composition.species()
            <<exit(FatalError);
    }

    volScalarField& rho = rhoFluid[i];
    volScalarField& K = KFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];
    PtrList<volScalarField>& Y = composition.Y();
    compressible::turbulenceModel& turb = turbulence[i];
    #include "createExplicitSources.H"
    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();
    volScalarField& hs = thermo.hs();

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
    forAll(Y,speciei)
    {
        fields.add(Y[speciei]);
    }
    fields.add(hs);
    DimensionedField<scalar,volMesh>& chemistrySh = chemistryShFluid[i];

