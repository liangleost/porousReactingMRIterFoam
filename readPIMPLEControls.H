    // We do not have a top-level mesh. Construct the fvSolution for
    // the runTime instead.
    fvSolution solutionDict(runTime);

    const dictionary& pimple = solutionDict.subDict("PIMPLE");

    int debugLevel(readInt(pimple.lookup("debugLevel")));
    int infoLevel(readInt(pimple.lookup("debugLevel")));
    bool stepIterLevel = pimple.lookupOrDefault<Switch>("stepIterLevel",true);
    bool fileWrite = pimple.lookupOrDefault<Switch>("fileWrite",true);
    int nOuterCorr(readInt(pimple.lookup("nOuterCorrectors")));
    double pConvergence = pimple.lookupOrDefault<double>("pConvergence",1e-8);
    double UxConvergence = pimple.lookupOrDefault<double>("UxConvergence",1e-8);
    double rhoConvergence = pimple.lookupOrDefault<double>("rhoConvergence",1e-8);
    double specieConvergence = pimple.lookupOrDefault<double>("specieConvergence",1.0);
    double hsConvergence = pimple.lookupOrDefault<double>("hsConvergence",1e-8);
    if(pimple.found("pConvergence")||pimple.found("UxConvergence")||pimple.found("rhoConvergence")||pimple.found("specieConvergence")||pimple.found("hsConvergence"))
    {
        Info<<"Convergence Criterion:\n	";
        if(pimple.found("pConvergence"))
        {
            Info<<"p - "<<pConvergence;
        }
        if(pimple.found("UxConvergence"))
        {
            Info<<"	Ux - "<<UxConvergence;
        }
        if(pimple.found("rhoConvergence"))
        {
            Info<<"	rho - "<<rhoConvergence;
        }
        if(pimple.found("specieConvergence"))
        {
            Info<<"	Species - "<<specieConvergence;
        }
        if(pimple.found("hsConvergence"))
        {
            Info<<"	T - "<<hsConvergence;
        }
        Info<<nl;
    }