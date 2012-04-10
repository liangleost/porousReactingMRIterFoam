    const dictionary& pimple = mesh.solutionDict().subDict("PIMPLE");

    int nCorr(readInt(pimple.lookup("nCorrectors")));

    int nNonOrthCorr =
        pimple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
    int debugLevel =
        pimple.lookupOrDefault<int>("debugLevel", 0);
    int infoLevel =
        pimple.lookupOrDefault<int>("infoLevel", 0);
    bool stepIterLevel =
        pimple.lookupOrDefault<Switch>("stepIterLevel", true);
    bool fileWrite =
        pimple.lookupOrDefault<Switch>("fileWrite", true);

    bool momentumPredictor =
        pimple.lookupOrDefault<Switch>("momentumPredictor", true);

    double pConvergence = 
        pimple.lookupOrDefault<double>("pConvergence",1e-8);
    double UxConvergence = 
        pimple.lookupOrDefault<double>("UxConvergence",1e-8);
    double rhoConvergence = 
        pimple.lookupOrDefault<double>("rhoConvergence",1e-8);
    double specieConvergence = 
        pimple.lookupOrDefault<double>("specieConvergence",1.0);
    double hsConvergence = 
        pimple.lookupOrDefault<double>("hsConvergence",1e-8);
    Info<<"pConvergence = "<<pConvergence<<"	UConvergence = "<<UxConvergence<<"	rhoConvergence = "<<rhoConvergence<<"	specieConvergence = "<<specieConvergence<<"	hsConvergence = "<<hsConvergence<<nl;