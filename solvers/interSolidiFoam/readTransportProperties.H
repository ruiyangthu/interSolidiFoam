    // set up properties for fluid - gas mixture
    immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

    volScalarField& alpha1(mixture.alpha1());
    volScalarField& alpha2(mixture.alpha2());

    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();
    const dimensionedScalar& cp1 = mixture.cp1();
    const dimensionedScalar& cp2 = mixture.cp2();
    const dimensionedScalar& lambda1 = mixture.lambda1();
    const dimensionedScalar& lambda2 = mixture.lambda2();
    const dimensionedScalar& beta1 = mixture.beta1();
    const dimensionedScalar& beta2 = mixture.beta2();
    
    const dictionary& meltingDict = mesh.solutionDict().subDict("MELTING");
    int minTCorr(readLabel(meltingDict.lookup("minTempCorrector")));
    int maxTCorr(readLabel(meltingDict.lookup("maxTempCorrector")));
    scalar alphaTol(readScalar(meltingDict.lookup("alphaTolerance")));
    scalar alphaRel(readScalar(meltingDict.lookup("alphaRelaxatio")));
   
 
    volScalarField nu = mixture.nu();

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Reading latent heat of fusion
    dimensionedScalar L(transportProperties.lookup("L"));

    // Reading melting temperature range
    dimensionedScalar Ts(transportProperties.lookup("Ts"));
    dimensionedScalar Tl(transportProperties.lookup("Tl"));

    dimensionedScalar Tmelt
    (
        "Tmelt",
        (Tl+Ts)/scalar(2.0)
    );

    // Reading large D'arcy-type source term constant DCl
    dimensionedScalar DCl(transportProperties.lookup("DCl"));
    dimensionedScalar DCs(transportProperties.lookup("DCs"));

/*
    // set up properties for solidification

    // Reading density rho
    dimensionedScalar rhos(transportProperties.lookup("rhos"));

    // Reading specific heat capacity
    dimensionedScalar cps(transportProperties.lookup("cps"));

    // Reading thermal conductivity
    dimensionedScalar lambdas(transportProperties.lookup("lambdas"));
*/

