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

Class
    eulerianStatistics

Description   
    Class for calculating eulerian statistics (means, variance) 
    from Lagrangian particle data
     
    Example Usage:
     
        controlled with the constant/cloudProperties dictionary:
        
        eulerianStatistics
        {
            enabled           true;

            sampleInterval      10; // sampling interval

            axisymmetric      true; // Circumferential average (default off)

            // Following settings are only required for axisymmetric=true
            rotAxis         (0 0 1);      

            planeType       pointAndNormal;
            triangulate     false;
            pointAndNormalDict
            {   
                point   (0 0 0);
                normal  (0 1.0 0);
            }

            // Fields to sample
            z;
            T;
            df;
            dx;
            U;
            rw;
        }
    
SourceFiles
    eulerianStatistics.H
    eulerianStatistics.C
  
\*---------------------------------------------------------------------------*/

#ifndef eulerianStatistics_H
#define eulerianStatistics_H

#include "volFields.H"
#include "IOdictionary.H"
#include "samplePlane.H"
#include "regIOobject.H"
#include "PstreamBuffers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Forward declaration of classes
class fvMesh;
  
/*---------------------------------------------------------------------------*\
                    Class eulerianStatistics Declaration
*----------------------------------------------------------------------------*/

class eulerianStatistics
:
    public regIOobject
{
    // Private data
        
        //- Enabled
        const Switch enabled_;
        
        //- Axissymmetric averaging (circular averaging)
        const Switch axisymmetric_;
        
        //- References to the mesh
        const fvMesh& mesh_;
        
        //- Sampling timestep interval
        const label sampleInterval_;
       
        //- sample plane required for axisymmetric averaging
        //  Not set if axisymmetric is switched off
        autoPtr<mmcSupport::samplePlane> plane_;
        
        //- Name list of properties to average
        DynamicList<word> names_;

        //- Current index of the field for averaging
        label ind_;
        
        //- Averages
        PtrList<Field<scalar>> ZAv_;
        PtrList<Field<scalar>> ZSqrAv_;
        
        //- Variance and standard deviation
        PtrList<Field<scalar>> ZVar_;
        PtrList<Field<scalar>> ZStdDev_;
        
        //- Weighting
        PtrList<Field<scalar>> ZwtAccum_;
        
        // For non axisymmetric averaging volScalarFields have to be created
        PtrList<volScalarField> volZAv_;
        PtrList<volScalarField> volZSqrAv_;
        PtrList<volScalarField> volZVar_;
        PtrList<volScalarField> volZStdDev_;
        PtrList<volScalarField> volZwtAccum_;
        
        
        //- Name extensions
        static const word ExtAv;
        static const word ExtSqrAv;
        static const word ExtVar;
        static const word ExtStdDev;
        static const word ExtZwtAccum;
        

    // Private Member Functions
        
        //- Write out field data if axisymmetric is switched on
        void axisymmetricWrite() const;
        
        //- Read field data if axisymmetric is switched on
        void axisymmetricRead
        (
            const word& name,
            const label i
        );
        
        //- Add fields for new property for axisymmetric averaging
        void addFieldsAxisymetric(const word& name);
        
        //- ADd fields for new property for non-axisymmetric averaging
        void addFields
        (
            const word& name,
            const dimensionSet& dim
        );

public:

    TypeName("eulerianStatistics");

    // Constructors

        //- Construct given name and mesh
        eulerianStatistics
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    // Member Functions
        
        inline Switch enabled() const {return enabled_;}
        
        //- Add new property
        void newProperty
        (
            const word& name,
            const dimensionSet& dim
        );

        //- Calculate mean and variance
        void calculate
        (
            const word& name,
            const scalar wt,
            const scalar Z
        );
        
        //- Set the cellIndex, required to be called before calculate
        void findCell(const point& pos);
        
        //- Write the eulerian statistics
        virtual bool writeData(Ostream& os) const;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
