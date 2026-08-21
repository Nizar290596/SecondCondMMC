# The cloudProperties File

The mmcDropletSprayFoam settings are controlled with the fuelCloudProperties
file. In the following an example file with the different sub-model options
are listed.

```C++
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2312                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      binary;
    arch        "LSB;label=32;scalar=64";
    class       dictionary;
    location    "constant";
    object      sprayCloudProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// The solution dictionary follows mostly the default OpenFOAM settings
// Note: integrationSchemes for T is different than the default. 
solution
{
    active          true;
    coupled         yes;
    cellValueSourceCorrection no;
    maxCo           0.3;
    calcFrequency   1;
    maxTrackTime    1;
    sourceTerms
    {
        schemes
        {
            U               semiImplicit 1;
        }
        resetOnStartup  false;
    }
    interpolationSchemes
    {
        rho             cell;
        U               cellPoint;
        p               cell;
        thermo:mu       cell;
        T               cell;
        Cp              cell;
        kappa           cell;
    }
    integrationSchemes
    {
        U               Euler;			// Euler Implicit method
        T               CrankNicolson;		// Choose here either
        					// CrankNicolson
        					// Euler Explicit
    }
}

subModels
{
    // Sub-models for the droplet parcel
    // =================================
    
    
    // Standard sub-model of the KinematicParcel
    // - Particle forces model (drag)
    // - injector models
    // - Dispersion model
    // - Patch interaction model
    // - Stochastic collision model
    // - Surface film model
    // - Packing model
    // - Damping model
    // - Exchange model (isotropy model)
    particleForces
    {
        sphereDrag      ;
    }
    // Injection model can be omitted 
    //    injectionModels     ;
    
    dispersionModel             none;
    patchInteractionModel       standardWallInteraction;
    stochasticCollisionModel    none;
    surfaceFilmModel            none;
    packingModel                none;
    dampingModel                none;
    isotropyModel               none;
    
    standardWallInteractionCoeffs
    {
        type            rebound;
    }
    
    // Sub-models for the ThermoSprayDroplet
    // =====================================
    
    // Liquid properties
    // - UDF: User defined functions
    // - OFLiquidProperties: Use OpenFOAMs liquid properties 
    LiquidPropertiesModel UDF;
    
    // How to couple the mmc stochastic particle to the source terms
    // of the droplet.
    dropletToMMCModel nearestNeighbor;
    
    
    // When the OpenFOAM liquid properties is selected, select
    // here the correct liquid name.
    // See also: ${WM_PROJECT_DIR}/src/thermophysicalProperties/liquidProperties
    OFLiquidPropertiesCoeffs
    {
        fluidName   C2H5OH;
    }
    
    
    // Use user defined functions for the liquid properties.
    // See also the liquidPropertiesModel class for more information!
    UDFCoeffs
    {
        saturationPressure
        {
            type            coded;
            code            #{
               std::vector<scalar> a {59.796, -6595.0, -5.0474, 6.3e-07, 2.0};  // C2H5OH
               return ( std::exp(a[0]+a[1]/x+a[2]*std::log(x)+a[3]*std::pow(x,a[4])) );
           #};
        }
        saturationTemperature
        {
            type            coded;
            code            #{
                std::vector<scalar> a {3.5154e+02, 2.5342e+01, 2.2231e+00, 2.0784e-01};  // C2H5OH
                scalar P1 = log(x/101325.0);
                return ( a[0]+a[1]*P1+a[2]*P1*P1+a[3]*P1*P1*P1 );
           #};
        }
        binaryDiffusionCoefficient 5.74e-05;
        latentHeat
        {
            type            coded;
            code            #{
           scalar Tc = 516.25;
           scalar T1 = x/Tc;
           std::vector<scalar> a {958345.091059, -0.4134, 0.75362, 0.0};  // C2H5OH
           return ( a[0]*pow(1.0-T1,a[1]+a[2]*T1+a[3]*T1*T1) );
           #};
        }
        liquidHeatCapacity
        {
            type            coded;
            code            #{
               scalar a_in = 3086.27;
               scalar b_in = -13.9017;
               scalar c_in = 0.0449525;
               scalar d_in = -1.89534803995077e-05;
               scalar e_in = 0;
               scalar f_in = 0;
               return ((((f_in*x+e_in)*x+d_in)*x +c_in)*x+b_in)*x+a_in;           #};
        }
        liquidAbsoluteEnthalpy
        {
            type            coded;
            code            #{
            std::vector<scalar> a {-6752827.25039109, 2052.57331394213, -0.060995463326749, -0.00238048724015426, 1.30130890620591e-05};  // C2H5OH
            return ( a[0]+a[1]*x+a[2]*x*x+a[3]*x*x*x+a[4]*x*x*x*x );
           #};
        }
        rhoL
        {
            type            coded;
            code            #{
                std::vector<scalar> a {70.1308387, 0.26395, 516.25, 0.2367, 0.0};  // C2H5OH
                const scalar T = std::min(a[2],x);
                return ( a[0]/pow(a[1],1.0+pow(1.0-T/a[2],a[3])) );
           #};
        }
    }
    
    // Coefficitions for the coupling of the stochastic particles
    // with the droplets
    nearestNeighborCoeffs
    {
        useSPProps      false;
        lambda_x        true;
        lambda_f        false;
        lambda_T        true;
        ri              0.00025856;
        fm              0.0072;
        Tm              40;
    }
    
    // Required to set the correct species index
    // Note: The dropletComponentName must match the name 
    //       in the chemistry model. This must not be identical to 
    //       the OpenFOAM liquidProperties name.
    dropletComponent C2H5OH;


    // Distribute Source Term
    // ======================

    // This is only possible for the DNS droplet spray simulation

    sourceDistribution  PSI;        // Default is particle-source-in-cell model (PSI)
                                    // (Optional) activate distributed

    // If distSourceModel is active the sourceDistrbutionCoeff need 
    // to provide the averaging volume
    distributedCoeffs 
    {
        // Type of the averaging volume, either box or sphere        
        type       box;
        // Potentially add additional entries required for the averaging
        // volume
    }
}

```




