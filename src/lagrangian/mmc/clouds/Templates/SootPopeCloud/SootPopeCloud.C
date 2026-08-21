/*---------------------------------------------------------------------------*\
                                       8888888888                              
                                       888                                     
                                       888                                     
  88888b.d88b.  88888b.d88b.   .d8888b 8888888  .d88b.   8888b.  88888b.d88b.  
  888 "888 "88b 888 "888 "88b d88P"    888     d88""88b     "88b 888 "888 "88b 
  888  888  888 888  888  888 888      888     888  888 .d888888 888  888  888 
  888  888  888 888  888  888 Y88b.    888     Y88..88P 888  888 888  888  888 
  888  888  888 888  888  888  "Y8888P 888      "Y88P"  "Y888888 888  888  888 
------------------------------------------------------------------------------- 

License
    This file is part of mmcFoam.

    mmcFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    mmcFoam is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with mmcFoam. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
#include "SootPopeCloud.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::SootPopeCloud<CloudType>::setModels()
{
    sootModel_.reset
    (
        SootModel<SootPopeCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    radiationModel_.reset
    (
        CloudRadiationModel<SootPopeCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


template<class CloudType>
void Foam::SootPopeCloud<CloudType>::cloudReset(SootPopeCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    sootModel_ = c.sootModel_;

    radiation_ = c.radiation_;
}


template<class CloudType>
void Foam::SootPopeCloud<CloudType>::setSootData()
{
    //- Read in information on soot species
    if (!this->cloudProperties().found("sootData"))
        return;

    const dictionary sootDict(this->cloudProperties().subDict("sootData"));

    const Switch sootSp(sootDict.lookup("sootingFlame"));

    sootingFlame_ = sootSp;

    if (sootingFlame_)
    {
        const scalar sootSD(readScalar(sootDict.lookup("sootSolidDensity")));

        sootSolidDensity_ = sootSD;

        List<scalar> sootMMWlist(sootDict.lookup("sootMinMolWt"));

        sootMinMolWt_ = sootMMWlist;

        numfVsoot_ = sootMinMolWt_.size();

        const scalar sootMKn(readScalar(sootDict.lookup("MinKnScaledMix")));

        sootMinKn_ = sootMKn;

        const scalar InterfV(readScalar(sootDict.lookup("IntermittencyfV")));

        IntermittencyfV_ = InterfV * 1e-9;

        Info << "Soot species included" << nl;

        Info << "Soot species minimum molcular weight is " << sootMinMolWt_ << " [kg/kmol]" << nl;

        Info << "Soot solid density is " << sootSolidDensity_ << " [kg/m^3]" << nl;

        Info << "Minimal Kn for non-scaled mixing is " << sootMinKn_ << nl << endl;

        Info << "Minimal fV of intermittency is " << IntermittencyfV_ << nl << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SootPopeCloud<CloudType>::SootPopeCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    const volVectorField& U,
    const volScalarField& DEff,
    const volScalarField& rho,
    const volVectorField& gradRho,
    const mmcVarSet& Xi,
    const Switch initAtCnstr,
    bool readFields
)
    :
    CloudType
    (
        cloudName,
        mesh,
        U,
        DEff,
        rho,
        gradRho,
        Xi,
        false,
        true
    ),

    sootPopeCloud(),

    sootModel_(nullptr),

    radiationModel_(nullptr),
    
    radiation_(false),

    sootingFlame_(false),

    sootSolidDensity_(),

    sootMinMolWt_(),

    numfVsoot_(1),

    sootMinKn_(),

    IntermittencyfV_(),

    indexO2N2_(2, 0)
{
    Info << nl << "Creating soot Pope Particle Cloud." << nl << endl;

    setModels();

    setSootData();

    setEulerianStatistics();

    if(readFields)
    {
        if(this->size()>0)
        {
            Info << nl << "Reading Soot Pope particle cloud data from file." << endl;
            particleType::sootParticleIOType::readFields(*this, this->composition());
        }
        else
        {
            if(initAtCnstr)
                this->initReleaseParticles();
        }
    }

    // To obtain the index of O2 and N2 in the mech
    indexO2N2_[0] = this->composition().slgThermo().carrier().species()["O2"];
    indexO2N2_[1] = this->composition().slgThermo().carrier().species()["N2"];
}


template<class CloudType>
Foam::SootPopeCloud<CloudType>::SootPopeCloud
(
    const word& cloudName,
    const fvMesh& mesh,
    const mmcVarSet& Xi,
    const Switch initAtCnstr,
    bool readFields
)
    :
    CloudType
    (
        cloudName,
        mesh,
        Xi,
        false,      // Only top level cloud calls initAtCnstr
        readFields
    ),

    sootPopeCloud(),

    sootModel_(nullptr),

    radiationModel_(nullptr),
    
    radiation_(false),

    sootingFlame_(false),

    sootSolidDensity_(),

    sootMinMolWt_(),

    numfVsoot_(1),

    sootMinKn_(),

    IntermittencyfV_(),

    indexO2N2_(2, 0)
{
    Info << nl << "Creating soot Pope Particle Cloud." << nl << endl;

    setModels();

    setSootData();

    setEulerianStatistics();

    if(readFields)
    {
        if(this->size()>0)
        {
            Info << nl << "Reading Soot Pope particle cloud data from file." << endl;
            particleType::sootParticleIOType::readFields(*this, this->composition());
        }
        else
        {
            if(initAtCnstr)
                this->initReleaseParticles();
        }
    }

    // To obtain the index of O2 and N2 in the mech
    indexO2N2_[0] = this->composition().slgThermo().carrier().species()["O2"];
    indexO2N2_[1] = this->composition().slgThermo().carrier().species()["N2"];
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SootPopeCloud<CloudType>::~SootPopeCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::SootPopeCloud<CloudType>::setParticleProperties
(
    particleType& particle,
    const scalar& mass,
    const scalar& wt,
    const scalar& patchI,
    const scalar& patchFace,
    const bool& iniRls
)
{
    CloudType::setParticleProperties(particle,mass,wt,patchI,patchFace,iniRls);

    particle.fVsoot() = calcfVsoot(particle.Y(),particle.T(),particle.pc());

    particle.Intermittency() = calcIntermittency(particle.fVsoot());

    particle.sootVf() = 0.;

    particle.sootIm() = 1.;

    particle.ySoot() = 0.;

    particle.nSoot() = 0.;

    particle.sootContribs() = {0., 0., 0., 0., 0., 0.};

    particle.XO2XN2() = 0.;

    particle.initStatisticalSampling();
}


template<class CloudType>
void Foam::SootPopeCloud<CloudType>::setEulerianStatistics()
{

    if (this->eulerianStatsDict().found("fVsoot"))
    {
        const dimensionSet dim = dimless;
        
        for(label ii = 0; ii< numfVsoot(); ii++)
        {
            this->eulerianStats().newProperty(word("fVsoot"+std::to_string(ii)),dim);
            Info << "created field = " << word("fVsoot"+std::to_string(ii)) << endl;
        }

    }

    if (this->eulerianStatsDict().found("Intermittency"))
    {
        const dimensionSet dim = dimless;
        
        for(label ii = 0; ii< numfVsoot(); ii++)
        {
            this->eulerianStats().newProperty(word("Intermittency"+std::to_string(ii)),dim);
        }

    }

    if (this->eulerianStatsDict().found("sootVf"))
    {
        const dimensionSet dim = dimless;

        this->eulerianStats().newProperty("sootVf", dim);
    }

    if (this->eulerianStatsDict().found("sootIm"))
    {
        const dimensionSet dim = dimless;

        this->eulerianStats().newProperty("sootIm", dim);
    }

    List<string> sootContribsTag = {"yInC2H2", "ySgC2H2", "yOxO2", "yOxOH", "nInC2H2", "nCg"};
    if(this->eulerianStatsDict().found("sootContribs"))
    {
        const dimensionSet dim = dimless;

        for(int i = 0; i < sootContribsTag.size(); i++)
        {
            this->eulerianStats().newProperty(word("sootContribs" + sootContribsTag[i]), dim);
        }
    }

    if (this->eulerianStatsDict().found("XO2XN2"))
    {
        const dimensionSet dim = dimless;

        this->eulerianStats().newProperty("XO2XN2", dim);
    }
}


template<class CloudType>
void Foam::SootPopeCloud<CloudType>::updateEulerianStatistics()
{
    CloudType::updateEulerianStatistics();

    forAllIters(*this, iter)
    {
        this->eulerianStats().findCell(iter().position());

        if (this->eulerianStatsDict().found("fVsoot"))
        {
            for(label ii = 0; ii< numfVsoot(); ii++)
            {
              this->eulerianStats().calculate("fVsoot"+std::to_string(ii),iter().wt(),iter().fVsoot()[ii]);
            }
        }

        if (this->eulerianStatsDict().found("Intermittency"))
        {
            for(label ii = 0; ii< numfVsoot(); ii++)
            {
              this->eulerianStats().calculate("Intermittency"+std::to_string(ii),iter().wt(),iter().Intermittency()[ii]);
            }
        }

        if(this->eulerianStatsDict().found("sootVf"))
        {
            this->eulerianStats().calculate("sootVf", iter().wt(), iter().sootVf());
        }

        if(this->eulerianStatsDict().found("sootIm"))
        {
            this->eulerianStats().calculate("sootIm", iter().wt(), iter().sootIm());
        }

        List<string> sootContribsTag = {"yInC2H2", "ySgC2H2", "yOxO2", "yOxOH", "nInC2H2", "nCg"};
        if(this->eulerianStatsDict().found("sootContribs"))
        {
            for(int i = 0; i < 6; i++)
            {
                this->eulerianStats().calculate("sootContribs" + sootContribsTag[i], iter().wt(), iter().sootContribs()[i]);
            }
        }

        if(this->eulerianStatsDict().found("XO2XN2"))
        {
            this->eulerianStats().calculate("XO2XN2", iter().wt(), iter().XO2XN2());
        }
    }
}


template<class CloudType>
Foam::scalarField Foam::SootPopeCloud<CloudType>::calcfVsoot
(
    const scalarField& Y,
    const scalar& T,
    const scalar& pc
)
{
    scalarField fVsoot(numfVsoot(),0.0);

    if (sootingFlame_)
    {
        const scalar rho = this->composition().particleMixture(Y).rho(pc,T);

        forAll(fVsoot,ii)
        {
            forAll(Y,ns)
            {
                if (this->composition().molWt(ns) >= sootMinMolWt_[ii]) 
                {
                    fVsoot[ii] += Y[ns] * rho / sootSolidDensity_;
                }
            }
        }
    }

    return fVsoot;
}

template<class CloudType>
Foam::scalarField Foam::SootPopeCloud<CloudType>::calcIntermittency
(
    const scalarField& fVsoot
)
{
    scalarField Intermittency(numfVsoot(),0.0);

    if (sootingFlame_)
    {
        forAll(Intermittency,ii)
        {
                if (fVsoot[ii] >= IntermittencyfV_) 
                {
                    Intermittency[ii] = 1.0;
                }
        }
    }

    return Intermittency;
}


template<class CloudType>
void Foam::SootPopeCloud<CloudType>::calcXO2XN2
(
    const scalarField& Y, 
    scalar& XO2N2
) const
{
    const double WO2 = 32.0;
    const double WN2 = 28.0;

    const scalar& YO2 = Y[this->indexO2N2()[0]];
    const scalar& YN2 = Y[this->indexO2N2()[1]];

    if(YN2 > SMALL)
    {
        XO2N2 = YO2*WN2/YN2/WO2;
    }
    else
    {
        XO2N2 = YO2*WN2/SMALL/WO2;
    }
}


template<class CloudType>
void Foam::SootPopeCloud<CloudType>::writeFields() const
{
    CloudType::writeFields();
    
    if (this->size())
    {
        particleType::sootParticleIOType::writeFields(*this, this->composition());
    }
}

