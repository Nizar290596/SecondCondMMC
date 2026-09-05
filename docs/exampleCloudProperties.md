# The cloudProperties File

Most of the settings in mmcFoam are controlled via the `cloudProperties`
file located in the constant/ directory of the OpenFOAM case. 

To give an overview of possible options and settings in the cloudProperties
an example file is provided below

```C++
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.0                                   |
|   \\  /    A nd           | Web:      http://www.openfoam.org               |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

FoamFile
{
    version         2.0;
    format          ascii;
   
    root            "";
    case            "";
    instance        "";
    local           "";
   
    class           dictionary;
    object          cloudProperties;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solution
{

    interpolationSchemes
    {
        U           cellPoint;
        DEff        cellPoint;
        rho         cellPoint;
        gradRho     cellPoint;
        f           cellPoint;
        p           cellPoint;
    }
}



// Controls the particle number within a super cell 

particleManagement
{
    numCtlOn  on;       // Switch to turn on/off particle management

    Nlo  10;            // Minimum number of particles in a supercell
    Npc  12;            // Target number of particles in a supercell
    Nhi  18;            // Maximum number of particles in a supercell
}


// In the following dictionary all submodels of mmcFoam
// are selected and controlled. 
subModels
{

    // (Optional) Specify thermodynamic model to use
    // Specify here the name of the thermophysical properties dictionary
    // to read. By defaul it selects the thermophysicalProperties dictionary
    // in the constant folder. 
    thermoModel     thermophysicalProperties;

// Select all submodels

    // Selects how particles enter the domain
    // Options:
    // 1. FreeStream
    // 2. NoInflow
    InflowBoundaryModel               FreeStream;


    // Composition model of the particles
    // Options:
    // 1. singlePhaseMixture
    compositionModel          singlePhaseMixture;


    // Radiation models - typically selected with soot
    // Options:
    // 1. OpticallyThinRadiation
    // 2. none
    radiationModel        OpticallyThinRadiation;

    // Select how particle particle pairs are mixed
    // Options:
    // 1. none
    // 2. MMCcurl
    // 3. sootMMCCurl -- special MMC curl model for soot
    // 4. ransMMCcurl -- for RANS, dense particle method
    // 5. Curl        -- for dense particle method
    mixingModel                          MMCcurl;


    // Aerosol sub-models
    // Options:
    // 1. none
    // 2. synthesisModel
    synthesisModel                          none;


    // Activate soot modeling
    // Options:
    // 1. none
    // 2. twoEquationModel
    sootModel                               none;


    // Thermophysical coupling of particles with the Eulerian fields
    // Options:
    // 1. none
    // 2. ParticleInCell
    // 3. KernelEstimation
    // 4. FlameletCurves
    thermoPhysicalCouplingModel KernelEstimation;
    

    // Reaction models
    // Options:
    // 1. none
    // 2. finiteRateNewtonLinODEChemistry
    // 3. finiteRateParticleChemistry
    // 4. flameSheetParticleChemistry
    // 5. sootParticleChemistry
    reactionModel finiteRateParticleChemistry;

    // Dynamic load balancing (finiteRateChem Only)
    // Note: Only works if it is run in parallel otherwise a FatalError occurs
    balanceReactionLoad     true;


// Define all settings of the submodels via subdictionary
// The notation is the name of the submodel and 'Coeffs' appended to it
// E.g., myNewModelCoeffs -- or general <modelName>Coeffs
// For detailed information about the options of a model see the *.H and *.C file 
// of the selected model

    
    KernelEstimationCoeffs
    {
        fLow                0.015;  // lower flamability limit
        fHigh               0.985;  // upper flamability limit
        fm                  0.030;  // same as in mixing model

        dfMax               0.015;  // maximum distance in f space allowed to do extrapolation
        C2                  0.1;    // used to control resolution

        nNearest            20;     // number of nearest neighbours per kernel (default 20)
        rMax                1.0e9;  // cap on the physical kernel radius

        // Conditioning axis for the kernel estimator. Defaults to
        // thermophysicalCoupling/condVariable if omitted. Set it to
        // phiModified to condition on the modified progress variable; any
        // other value must name a registered coupling variable.
        //
        // Do NOT put phiModified in thermophysicalCoupling/condVariable:
        // that entry is also handed to the reaction model as the name of the
        // mixture-fraction coupling variable and must stay z.
        condVariable        phiModified;
    }


    OpticallyThinRadiationCoeffs
    {
    }

    MMCcurlCoeffs
    {
        //- Parameters r_i and f_m as per Eq.(29) Cleary & Klimenko, Phys Fluids 23, 115102, 2011
        // (0.5*(3120*2.42e-09/6/(0.000175)**(2-2.36)*1/0.03)**(1/2.36))/(3**0.5)
        r_i                     0.00193;

        // One normalisation factor per pairing axis, keyed "<axis>_m".
        // The axes of the first conditioning are the reference variables
        // declared in mmcVariablesDefinitions, so with shadow positions:
        Xim_i
        {
            sPx_m                  0.03;
            sPy_m                  0.03;
            sPz_m                  0.03;
        }

        aISO                       true;    // Select either aISO false or true -- if false it
                                            // defaults to C&K model
                                            // NOTE: the C&K time scale needs a reference
                                            // variable with a transported Eulerian field
                                            // (referenceType interpolated). With shadow
                                            // positions its |grad Xi|^2 is identically zero
                                            // and no pair would mix.
        
        meanTimeScale              true;    // Use the mean of the two particle times
                                            // if set to false it uses the min()

        pairingMethod              global;  // Define the region where particles are mixed
                                            // Options:
                                            // 1. global     -- all particles can mix
                                            // 2. local      -- only particles on the 
                                            //                  processor can mix
                                            // 3. subVolumes -- create subVolumes in which particles
                                            //                  are mixed
    }


    finiteRateParticleChemistryCoeffs
    {
        zLower              0.01;
        zUpper              0.99;
    }

}


// ===========================================================================
//                      Second conditioning (optional)
// ===========================================================================
// A second, independent mixing stage that acts on a random subset of the
// cloud and pairs it in its own reference space.
//
// Each stage is selected by its own keyword and configured by its own
// coefficients sub-dictionary:
//
//     mixingModel             MMCcurl;
//     secondCondMixingModel   secondCondMMCcurl;
//
// While secondConditioning/enabled is true, MMCcurl restricts the first
// conditioning to the progress variable phi, which every particle carries; the
// composition of the flagged subset is then mixed by secondCondMMCcurl in its
// own reference space. Mixing the composition in both stages would mix those
// scalars twice, so the first stage stands down from it.
//
// Consequence: with the gate enabled the composition is not mixed by the first
// stage at all, so particles outside the subset keep their injection
// composition. They also skip the chemistry and are excluded from the
// thermophysical coupling and from the T/Y Eulerian statistics. Set
// enabled false for the usual single-level behaviour.

secondConditioning
{
    enabled             true;   // constructs the second mixing model

    R                    0.1;   // probability that a new particle joins the
                                // subset (secondCondFlag = 1)

    Tu                 300.0;   // phi = (T - Tu)/(Tb - Tu) at injection,
    Tb                2000.0;   // clipped to [0,1]

    beta                 1.0;   // phi_deg = phi * exp(beta * omega_OU)
    tauOU              1.0e-3;  // OU relaxation time [s]

    A_phi                0.0;   // W(phi) = A*(1-phi)*exp[Z*(phi-1)]
    Z_phi                0.0;
}

// ... and inside subModels, alongside the usual MMCcurlCoeffs:
//
//     secondCondMMCcurlCoeffs
//     {
//         r_i             1.0e-3;
//
//         // The axes this stage pairs on, in order. Any scalar registered on
//         // the particle may be used. Defaults to phiModified followed by the
//         // mmcVarSet reference variables.
//         referenceAxes   ( phiModified sPx sPy sPz );
//
//         // One entry per axis above, keyed "<axis>_m". A tight
//         // phiModified_m makes the k-d tree split preferentially on the
//         // progress variable.
//         Xim_i
//         {
//             phiModified_m   0.01;
//             sPx_m           0.05;
//             sPy_m           0.05;
//             sPz_m           0.05;
//         }
//
//         particleFilter  secondCondFlag;   // or: none
//         pairingMethod   local;            // local | subVolumes | global
//
//         meanTimeScale   true;
//         CL              0.5;
//         CE              0.1;
//
//         // Diagnostic pair log. Samples this many paired couples (globally,
//         // not per rank) at each write time and writes them to
//         //   postProcessing/secondCondPairs/secondCondPairs_<time>.dat
//         // Columns: proc, both secondCondFlags, both phiModified values and
//         // their separation, the separation over the remaining reference
//         // axes, the physical separation, and both temperatures before and
//         // after the mixing. 0 (default) disables it.
//         nPairSamples    1000;
//     }


// This is the subdictionary of the thermophysical model
// Note: This is a bug and in the future might be moved 
//       to e.g. KernelEstimationCoeffs
thermophysicalCoupling
{
    CH3OCH3;
    CH3OCH2O2;
    CH3OCHO;
    C2H2;
    C2H4;
    C2H6;
    CH2O;
    CH3;
    CH4;
    O;
    O2;
    N2;
    H2;
    H2O;
    CO;
    CO2;
    OH;
    AR;

    // The coupling variable passed to the reaction model as the mixture
    // fraction. Must name a registered coupling variable - it is looked up as
    // XiC(condVariable) in ReactingPopeParticle::calc(). To condition a
    // coupling model on something else, use that model's own condVariable
    // entry (see KernelEstimationCoeffs).
    condVariable         z;

    // The blending of the relaxation timescale allows 
    // for a smooth initial transient simulation period.
    // After one domain flush tauRelaxStart should be 
    // set to 200 and tauRelax to the commonly used 
    // factor of 10.
    tauRelaxBlending false;
    tauRelax            20;
    tauRelaxStart      200;
    tauRelaxDelta      0.5;

    tauUnits      timestep;
}



// ===========================================================================
//                      Statistic Tools of mmcFoam
// ===========================================================================
// It is possible to collect statistics of mmcFoam during the run time 
// of the simulation. Either by averaging the particle data to the Eulerian
// field or by sampling the particle data in a defined probe volume. 


eulerianStatistics
{
    enabled           true;

    sampleInterval      10;     // Every n time steps

    axisymmetric      true;     // If true a plane is constructed for averaging

    // Entries only required for axisymmetric true:
        rotAxis        (0 0 1);

        planeType       pointAndNormal;
        triangulate     false;
        pointAndNormalDict
        {
            point   (0 0 0);
            normal  (0.0001 1.0 0); // Small offset so cells are cut better
        }

    // List here all variables that you want to sample
    z;
    T;
    df;
    dx;
    U;
}

particleSampling
{
    enabled                  true;

    sampleInterval             10;
    
    condVariable                z;  // Define the conditioning variable 
                                    // required for lowerBound and upperBound

    lowerBound               0.00;  // Exclude particles at or below this 
                                    // threshold of the reference variable
    upperBound               1.00;  // Exclude particles at or above this 
                                    // threshold of the reference variable.

    fields                  (H T);  // Optional list of fields to sample.
                                    // if fields is left blank all fields
                                    // are sampled

    collated                 true;  // Write the output in a collated binary format
                                    // Default is: true

    // Define the particle sampling probes. 
    // Currently only particleProbeBox is supported, which spans up a box
    // with the two points p1 and p2.

    particleSampleProbes
    {
        zD05    // Name of the probe 
        {
            type    particleProbeBox;
            p1      (-1 -1 0.03675);
            p2      ( 1 1  0.03775);
        }
        zD075   // Name of the probe
        {
            type    particleProbeBox;
            p1      (-1 -1 0.055825);
            p2      ( 1 1  0.055925);
        }
    }

}

// ************************************************************************* //

```




