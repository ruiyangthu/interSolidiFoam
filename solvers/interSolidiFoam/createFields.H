Info<< "Reading field T\n" << endl;
volScalarField T
(
 IOobject
 (
  "T",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
 IOobject
 (
  "p_rgh",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );

Info<< "Reading field U\n" << endl;
volVectorField U
(
 IOobject
 (
  "U",
  runTime.timeName(),
  mesh,
  IOobject::MUST_READ,
  IOobject::AUTO_WRITE
  ),
 mesh
 );





#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
#include "readTransportProperties.H"
//myimmiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

Info<< "Calculating phase-fraction\n" << endl;
//Calculating phase-fraction f of melt 
//volScalarField& alpha1(mixture.alpha1());
//volScalarField& alpha2(mixture.alpha2());

//const dimensionedScalar& rho1 = mixture.rho1();
//const dimensionedScalar& rho2 = mixture.rho2();
volScalarField fl ("fl", max(min((T-Ts)/(Tl-Ts),scalar(1)),scalar(0)));
volScalarField fs ("fs", scalar(1) - fl);

Info<< "Liquid Melt volume fraction fl = "
<< fl.weightedAverage(mesh.Vsc()).value()
<< "  Min(fl) = " << min(fl).value()
<< "  Max(fl) = " << max(fl).value()
<< endl;

Info<< "Solid Melt volume fraction fs = "
<< fs.weightedAverage(mesh.Vsc()).value()
<< "  Min(fs) = " << min(fs).value()
<< "  Max(fs) =\n" << max(fs).value()
<< endl;


volScalarField alphas 
(
 IOobject
 (
  "alphas",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 max(min(alpha1*fs,scalar(1)),scalar(0))
 );
alphas.write();
// Need to store rho for ddt(rho, U)
volScalarField rho
(
 IOobject
 (
  "rho",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT
  ),
 alpha1*rho1 + alpha2*rho2,
 alpha1.boundaryField().types()
 );
rho.oldTime();

// Mass flux
surfaceScalarField rhoPhi
(
 IOobject
 (
  "rhoPhi",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 fvc::interpolate(rho)*phi
 );

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
 incompressible::turbulenceModel::New(U, phi, mixture)
 );

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p
(
 IOobject
 (
  "p",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 p_rgh + rho*gh
 );

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
 p,
 p_rgh,
 pimple.dict(),
 pRefCell,
 pRefValue
 );

if (p_rgh.needReference())
  {
    p += dimensionedScalar
      (
       "p",
       p.dimensions(),
       pRefValue - getRefCellValue(p, pRefCell)
       );
    p_rgh = p - rho*gh;
  }
mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES flux from previous time-step
surfaceScalarField alphaPhi
(
 IOobject
 (
  "alphaPhi",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
 phi*fvc::interpolate(alpha1)
 );

// MULES Correction
tmp<surfaceScalarField> talphaPhiCorr0;

//Creating fields for thermophysical properties
Info<< "Creating fields for thermophysical properties\n" << endl;

volScalarField DC
(
 IOobject
 (
  "DC",
  runTime.timeName(),
  mesh
  ),
 alpha1*DCl*pow(scalar(1)-fl,2)/(pow(fl,3)+DCs),
 alpha1.boundaryField().types()
 );
DC.oldTime();

volScalarField beta
(
 IOobject
 (
  "beta",
  runTime.timeName(),
  mesh
  ),
 alpha1*beta1 + alpha2*beta2,
 alpha1.boundaryField().types()
 );
beta.oldTime();

volScalarField lambda
(
 IOobject
 (
  "lambda",
  runTime.timeName(),
  mesh
  ),
 alpha1*lambda1 + alpha2*lambda2,
 alpha1.boundaryField().types()
 );
lambda.oldTime();

volScalarField cp
(
 IOobject
 (
  "cp",
  runTime.timeName(),
  mesh
  ),
 alpha1*cp1 + alpha2*cp2,
 alpha1.boundaryField().types()
 );
cp.oldTime();

volScalarField Tcorr("Tcorr", (Tl-Ts)*fl + Ts);

dimensionedScalar dt = runTime.deltaT();
volScalarField Sp ("Sp", -rho*L*alpha1/dt/(Tl-Ts));
volScalarField Sc ("Sc", -Sp*Tcorr);


// Need to store rhoCp for ddt(rhoCp, U)
volScalarField rhoCp
(
 IOobject
 (
  "rhoCp",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT
  ),
 alpha1*rho1*cp1 + alpha2*rho2*cp2,
 alpha1.boundaryField().types()
 );
rhoCp.oldTime();

surfaceScalarField phiCp
(
 IOobject
 (
  "phiCp",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 phi*fvc::interpolate(cp)
 );

volScalarField rho1Alpha1
(
 IOobject
 (
  "rho1Alpha1",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 rho1*alpha1,
 alpha1.boundaryField().types()
 );
rho1Alpha1.oldTime();


//Calculating enthalpy and corrected temperature
Info<< "Calculating enthalpy and corrected temperature\n" << endl;

volScalarField h
(
 IOobject
 (
  "h",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 cp*(T-Tmelt) + alpha1*fl*L
 );
h.write();

