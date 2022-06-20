const dictionary& pimpleCHT = mesh.solutionDict().subDict("PIMPLE");

int nNonOrthCorr =
    pimpleCHT.getOrDefault<int>("nNonOrthogonalCorrectors", 0);
