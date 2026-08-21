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
#include "eulerianStatistics.H"
#include "IFstream.H"

namespace Foam
{
defineTypeNameAndDebug(eulerianStatistics, 0);

const word Foam::eulerianStatistics::ExtAv = "Av";
const word Foam::eulerianStatistics::ExtSqrAv = "SqrAv";
const word Foam::eulerianStatistics::ExtVar = "Var";
const word Foam::eulerianStatistics::ExtStdDev = "StdDev";
const word Foam::eulerianStatistics::ExtZwtAccum = "wt";
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct given name and mesh
Foam::eulerianStatistics::eulerianStatistics
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    regIOobject
    (
        IOobject
        (
            "eulerianStatistics",
            mesh.time().constantPath(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    ),
    
    enabled_(readBool(dict.lookup("enabled"))),

    axisymmetric_(dict.lookupOrDefault<Switch>("axisymmetric",false)),
        
    mesh_(mesh),
    
    sampleInterval_(dict.lookupOrDefault("sampleInterval",1)),
       
    ZAv_(),
    
    ZSqrAv_(),
 
    ZVar_(),
        
    ZStdDev_(),

    ZwtAccum_()
{
    if (enabled_)
        Info << "Eulerian Statistics:"
             << dict << endl;
    if (axisymmetric_ && enabled_)
        plane_.reset(new mmcSupport::samplePlane(mesh,dict));
}


// * * * * * * * * * * * * *  Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * *  Member Functions  * * * * * * * * * * * //

void Foam::eulerianStatistics::newProperty
(
    const word& name,
    const dimensionSet& dim
)
{
    if (!enabled_)
        return;
    
    names_.append(name);

    if (axisymmetric_)
        addFieldsAxisymetric(name);
    else
        addFields(name,dim);
}


void Foam::eulerianStatistics::addFieldsAxisymetric(const word& name)
{
    const label fieldSize = plane_->size();

    ZAv_.append
    (
        new Field<scalar>(fieldSize,0)
    );
    
    ZSqrAv_.append
    (
        new Field<scalar>(fieldSize,0)
    );
        
    ZVar_.append
    (
        new Field<scalar>(fieldSize,0)
    );
    
    ZStdDev_.append
    (
        new Field<scalar>(fieldSize,0)
    );
    
    ZwtAccum_.append
    (
        new Field<scalar>(fieldSize,0)
    );
    
    // Attempt to read fields if present 
    // Required for restart runs
    axisymmetricRead(name,ZAv_.size()-1);
};


void Foam::eulerianStatistics::axisymmetricRead
(
    const word& name,
    const label i
)
{
    // Check if field exist 
    fileName samplePath
    (
        mesh_.time().globalPath()/"mmcStatistics/eulerianStatistics/"
      + mesh_.time().timeName()
    );
    
    if(Foam::exists(samplePath/name+ExtAv))
    {
        // Only the master processor reads the data
        if (Pstream::parRun() && !Pstream::master())
            return;
        
        Info << "Read eulerianStatistics data for "<<name 
             << " from "<<samplePath<<endl;
        
        IFstream ifsZAv(samplePath/name+ExtAv);
        ifsZAv >> ZAv_[i];
        
        IFstream ifsZSqrAv(samplePath/name+ExtSqrAv);
        ifsZSqrAv >> ZSqrAv_[i];
        
        IFstream ifsZVar(samplePath/name+ExtVar);
        ifsZVar >> ZVar_[i];
        
        IFstream ifsZStdDev(samplePath/name+ExtStdDev);
        ifsZStdDev >> ZStdDev_[i];
        
        IFstream ifsZwtAccum(samplePath/name+ExtZwtAccum);
        ifsZwtAccum >> ZwtAccum_[i];
    }
}


void Foam::eulerianStatistics::addFields
(
    const word& name,
    const dimensionSet& dim
)
{

    // Note: We have to create a volScalarField and a normal field
    //       Using a reference to the internal volScalarField does not work
    //       as this leads to a double free corruption error during the 
    //       destruction of the PtrLists.

    volZAv_.append
    (
        new volScalarField
        (
            IOobject
            (
                name+ExtAv,
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                name+ExtAv,
                dim,
                0
            )
        )
    );
    
    ZAv_.append(new Field<scalar>(volZAv_.last().field()));
    
    volZSqrAv_.append
    (
        new volScalarField
        (
            IOobject
            (
                name+ExtSqrAv,
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                name+ExtSqrAv,
                dim*dim,
                0
            )
        )
    );
    
    ZSqrAv_.append(new Field<scalar>(volZSqrAv_.last().field()));
        
    volZVar_.append
    (
        new volScalarField
        (
            IOobject
            (
                name+ExtVar,
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                name+ExtVar,
                dim*dim,
                0
            )
        )
    );
        
    ZVar_.append(new Field<scalar>(volZVar_.last().field()));
        
    volZStdDev_.append
    (
        new volScalarField
        (
            IOobject
            (
                name+ExtStdDev,
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                name+ExtStdDev,
                dim,
                0
            )
        )
    );
        
    ZStdDev_.append(new Field<scalar>(volZStdDev_.last().field()));
        
    volZwtAccum_.append
    (
        new volScalarField
        (
            IOobject
            (
                name+ExtZwtAccum,
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                name+ExtZwtAccum,
                dimless,
                0
            )
        )
    );
    
    ZwtAccum_.append(new Field<scalar>(volZwtAccum_.last().field()));
}


void Foam::eulerianStatistics::calculate
(
    const word& name,
    const scalar wt,
    const scalar Z
)
{
    if (!enabled_)
        return;

    const label timeStep = mesh_.time().timeIndex();
    
    if (timeStep % sampleInterval_)
        return;
    
    //- Find correct field
    label iField = -1;
    
    forAll(names_,i)
    {
        if (name == names_[i])
            iField = i;
    }
    
    if (iField == -1 || ind_ == -1)
        return;
    
    //- Calculate averages
    const scalar deltaT = mesh_.time().deltaT().value();

    ZwtAccum_[iField][ind_] += wt * deltaT;

    ZAv_[iField][ind_] += (wt * deltaT) /ZwtAccum_[iField][ind_] * (Z - ZAv_[iField][ind_]);
    
    ZSqrAv_[iField][ind_] += (wt * deltaT) /ZwtAccum_[iField][ind_] * (sqr(Z) - ZSqrAv_[iField][ind_]);

    ZVar_[iField][ind_] = max(ZSqrAv_[iField][ind_] - sqr(ZAv_[iField][ind_]), 0.0);

    ZStdDev_[iField][ind_] = sqrt(ZVar_[iField][ind_]);

    // Only required for non-axisymmetric averaging 
    if (!axisymmetric_ && mesh_.time().writeTime())
    {
        // Copy the field data into the volScalarFields
        volZwtAccum_[iField][ind_] = ZwtAccum_[iField][ind_];
        volZAv_[iField][ind_] = ZAv_[iField][ind_];
        volZSqrAv_[iField][ind_] = ZSqrAv_[iField][ind_];
        volZVar_[iField][ind_] = ZVar_[iField][ind_];
        volZStdDev_[iField][ind_] = ZStdDev_[iField][ind_];
    }
}


void Foam::eulerianStatistics::findCell(const point& p)
{
    if (!enabled_)
        return;

    if (axisymmetric_)
        ind_ = plane_->cellIndex(p);
    else
        ind_ = mesh_.findCell(p);
}


void Foam::eulerianStatistics::axisymmetricWrite() const
{
    if (Pstream::parRun())
    {
        // Combine first all results to one processor
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        if (!Pstream::master())
        {
            UOPstream toBuffer(Pstream::masterNo(),pBufs);
            forAll(names_,i)
            {
                toBuffer << ZAv_[i] 
                         << ZSqrAv_[i] 
                         << ZwtAccum_[i];
            }
        }

        pBufs.finishedSends();


        if (Pstream::master())
        {
            // Only the master has to read
            forAll(names_,i)
            {
                // Need to create local copies as function has to be const
                Field<scalar> ZAv = ZAv_[i];
                Field<scalar> ZSqrAv = ZSqrAv_[i];
                Field<scalar> ZwtAccum = ZwtAccum_[i];

                for (label procI=0; procI < Pstream::nProcs(); procI++)
                {
                    if (procI == Pstream::masterNo())
                        continue;

                    UIPstream fromBuffer(procI,pBufs);
                    Field<scalar> tempZAv;
                    Field<scalar> tempZSqrAv;
                    Field<scalar> tempZwtAccum;

                    fromBuffer >> tempZAv;
                    fromBuffer >> tempZSqrAv;
                    fromBuffer >> tempZwtAccum;

                    // It can happen that not all cells in the plane have a 
                    // value assigned to them if e.g. no particle is located
                    // there. Then this would lead to a division by zero. 
                    // --> Requires manual for loop
                    forAll(ZAv,cellI)
                    {
                        if ((ZwtAccum[cellI]+tempZwtAccum[cellI]) > SMALL)
                        {
                            ZAv[cellI] = 
                                (
                                    (ZAv[cellI]* ZwtAccum[cellI]) 
                                  + (tempZAv[cellI]*tempZwtAccum[cellI])
                                )/ (ZwtAccum[cellI]+tempZwtAccum[cellI]);
                            
                            ZSqrAv[cellI] = 
                                (
                                    (ZSqrAv[cellI]* ZwtAccum[cellI]) 
                                  + (tempZSqrAv[cellI]*tempZwtAccum[cellI])
                                )/ (ZwtAccum[cellI]+tempZwtAccum[cellI]);
                        }
                    }
                    
                    ZwtAccum += tempZwtAccum;
                }

                Field<scalar> ZVar = max(ZSqrAv - sqr(ZAv),0.0);
                Field<scalar> ZStdDev = sqrt(ZVar);
                
                plane_->writeField(names_[i]+ExtAv,ZAv);
                plane_->writeField(names_[i]+ExtSqrAv,ZSqrAv);
                plane_->writeField(names_[i]+ExtVar,ZVar);
                plane_->writeField(names_[i]+ExtStdDev,ZStdDev);
                plane_->writeField(names_[i]+ExtZwtAccum,ZwtAccum);
            }
        }
    }
    else
    {
        // write the fields as vtp files
        forAll(names_,i)
        {
            plane_->writeField(names_[i]+ExtAv,ZAv_[i]);
            plane_->writeField(names_[i]+ExtSqrAv,ZSqrAv_[i]);
            plane_->writeField(names_[i]+ExtVar,ZVar_[i]);
            plane_->writeField(names_[i]+ExtStdDev,ZStdDev_[i]);
            plane_->writeField(names_[i]+ExtZwtAccum,ZwtAccum_[i]);
        }
    }
}


bool Foam::eulerianStatistics::writeData(Ostream& os) const 
{
    // If axisymmetric is switched off, the fields are written with the
    // object registry write, as volScalarFields are registered with the 
    // object registry
    
    
    if (axisymmetric_)
    {
        os << "Write axisymmetric file" << endl;
        axisymmetricWrite();
    }
    else
    {
        os << "Writing eulerian fields"<<endl;
    }
   
    os.check(FUNCTION_NAME);
    return os.good();
}



// ************************************************************************* //
