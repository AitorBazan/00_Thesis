    // We do not have a top-level mesh. Construct the fvSolution for
    // the runTime instead.
    fvSolution solutionDict(runTime);

    const dictionary& pimpleCHT = solutionDict.subDict("PIMPLE");

    const int nOuterCorr =
        pimpleCHT.getOrDefault<int>("nOuterCorrectors", 1);
        