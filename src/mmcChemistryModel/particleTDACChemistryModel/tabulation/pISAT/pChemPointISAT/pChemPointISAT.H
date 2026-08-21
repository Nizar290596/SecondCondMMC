/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::pChemPointISAT

Description
    Adaption necessary for particleTDACChemistryModel

\*---------------------------------------------------------------------------*/

#ifndef pChemPointISAT_H
#define pChemPointISAT_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
class pBinaryNode;

template<class ReactionThermo, class ThermoType>
class particleTDACChemistryModel;


/*---------------------------------------------------------------------------*\
                       Class pChemPointISAT Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class pChemPointISAT
{
    // Private data

        //- Pointer to the chemistryModel object
        particleTDACChemistryModel<ReactionThermo, ThermoType>& chemistry_;

        //- Vector storing the composition, temperature and pressure
        //  and deltaT if a variable time step is set on
        scalarField phi_;

        //- Vector storing the mapping of the composition phi
        scalarField Rphi_;

        //- LT the transpose of the L matrix describing the Ellipsoid Of
        //  Accuracy use List of Lists to be able to change size if DAC is used
        scalarSquareMatrix LT_;

        //- A the mapping gradient matrix
        scalarSquareMatrix A_;

        //- Vector storing the scale factor
        scalarField scaleFactor_;

        //- Reference to the node in the binary tree holding this chemPoint
        pBinaryNode<ReactionThermo, ThermoType>* node_;

        //- The size of the composition space (size of the vector phi)
        label completeSpaceSize_;

        //- Number of times the element has been grown
        label nGrowth_;

        //- Tolerance for the Ellipsoid of accuracy
        static scalar tolerance_;

        //- Number of active species stored in the chemPoint
        label nActiveSpecies_;

        //- Vectors that store the index conversion between the simplified
        //  and the complete chemical mechanism
        List<label> simplifiedToCompleteIndex_;

        label timeTag_;

        //for testing
        label lastTimeUsed_;
        label TimeNotUsed_;
        static label leafCount_;


        bool toRemove_;

        label maxNumNewDim_;

        Switch printProportion_;

        //- Variable to store the number of retrieves the chemPoint
        //  will generate at each time step
        label numRetrieve_;

        //- Variable to store the number of time steps the chempoint is allowed
        //   to still live according to the maxChPLifeTime_ parameter
        label nLifeTime_;

        List<label> completeToSimplifiedIndex_;

        //- Number of equations in addition to the species eqs.
        label nAdditionalEqns_;

        label idT_;
        label idp_;
        label iddeltaT_;

        //- QR decomposition of a matrix
        //  Input : nCols cols number
        //  R the matrix to decompose
        //  QT an empty matrix that stores the transpose of the Q matrix
        //  R is returned in the given R matrix
        //  which is used to store the ellipsoid of accuracy
        void qrDecompose
        (
            const label nCols,
            scalarSquareMatrix& R
        );

        //- QR update of the matrix A
        void qrUpdate
        (
            scalarSquareMatrix& R,
            const label n,
            const scalarField& u,
            const scalarField& v
        );

        void rotate
        (
            scalarSquareMatrix& R,
            const label i,
            const scalar a,
            const scalar b,
            label n
        );


public:

    // Constructors

        //- Construct from components
        pChemPointISAT
        (
            particleTDACChemistryModel<ReactionThermo, ThermoType>& chemistry,
            const scalarField& phi,
            const scalarField& Rphi,
            const scalarSquareMatrix& A,
            const scalarField& scaleFactor,
            const scalar& tolerance,
            const label& completeSpaceSize,
            const dictionary& coeffsDict,
            pBinaryNode<ReactionThermo, ThermoType>* node = nullptr
        );

        //- Construct from another chemPoint and reference to a binary node
        pChemPointISAT
        (
            const pChemPointISAT<ReactionThermo, ThermoType>& p,
            pBinaryNode<ReactionThermo, ThermoType>* node
        );

        //- Construct from another chemPoint
        pChemPointISAT
        (
            pChemPointISAT<ReactionThermo, ThermoType>& p
        );

//for testing 
        scalar checkSolutionError_;
        label leafID_;


    // Member functions

        //- Access to the particleTDACChemistryModel
        inline particleTDACChemistryModel<ReactionThermo, ThermoType>& chemistry()
        {
            return chemistry_;
        }

        inline label nGrowth()
        {
            return nGrowth_;
        }

        inline label& completeSpaceSize()
        {
            return completeSpaceSize_;
        }

        inline const scalarField& phi() const
        {
            return phi_;
        }

        inline const scalarField& Rphi() const
        {
            return Rphi_;
        }

        inline const scalarField& scaleFactor()
        {
            return scaleFactor_;
        }

        inline const scalar& tolerance()
        {
            return tolerance_;
        }

        static void changeTolerance(scalar newTol)
        {
            tolerance_ = newTol;
        }

        inline pBinaryNode<ReactionThermo, ThermoType>*& node()
        {
            return node_;
        }

        inline const scalarSquareMatrix& A() const
        {
            return A_;
        }

        inline scalarSquareMatrix& A()
        {
            return A_;
        }

        inline const scalarSquareMatrix& LT() const
        {
            return LT_;
        }

        inline scalarSquareMatrix& LT()
        {
            return LT_;
        }

        inline label nActiveSpecies()
        {
            return nActiveSpecies_;
        }

        inline List<label>& completeToSimplifiedIndex()
        {
            return completeToSimplifiedIndex_;
        }

        inline List<label>& simplifiedToCompleteIndex()
        {
            return simplifiedToCompleteIndex_;
        }

        //- Increases the number of retrieves the chempoint has generated
        void increaseNumRetrieve();

        //- Resets the number of retrieves at each time step
        void resetNumRetrieve();

        //- Increases the "counter" of the chP life
        void increaseNLifeTime();

        label simplifiedToCompleteIndex(const label i);

        inline const label& timeTag()
        {
            return timeTag_;
        }

        inline label& lastTimeUsed()
        {
            return lastTimeUsed_;
        }

        inline label& TimeNotUsed()
        {
            return TimeNotUsed_;
        }

        inline bool& toRemove()
        {
            return toRemove_;
        }

        inline label& maxNumNewDim()
        {
            return maxNumNewDim_;
        }

        inline const label& numRetrieve()
        {
            return numRetrieve_;
        }

        inline const label& nLifeTime()
        {
            return nLifeTime_;
        }

        inline bool variableTimeStep() const
        {
            return chemistry_.variableTimeStep();
        }

        // ISAT functions

            //- To RETRIEVE the mapping from the stored chemPoint phi, the query
            // point phiq has to be in the EOA of phi.
            // To test if phiq is in the ellipsoid:
            // ||L^T.dphi|| <= 1
            bool inEOA(const scalarField& phiq);

            //- More details about the minimum-volume ellipsoid covering an
            //  ellispoid E and a point p are found in [1].
            //  Here is the main steps to obtain the
            //  modified matrix L' describind the new ellipsoid.
            //  1) calculate the point p' in the transformed space :
            //  p' = L^T.(p-phi)
            //  2) compute the rank-one decomposition:
            //  G = I + gamma.p'.p'^T,
            //  with gamma = (1/|p'|-1)*1/|p'|^2
            //  3) compute L':
            //  L'L'^T = (L.G)(L.G)^T,
            //  L'^T is then obtained by QR decomposition of (L.G)^T = G^T.L^T
            //  [1] Stephen B. Pope, "Algorithms for ellipsoids", FDA 08-01,
            //  Cornell University, 2008
            bool grow(const scalarField& phiq);

            //- If phiq is not in the EOA, then the mapping is computed.
            //  But as the EOA is a conservative approximation of the region of
            //  accuracy surrounding the point phi, we could expand it by
            //  comparing thecomputed results with the one obtained by linear
            //  interpolation.  The error eps is calculated:
            //  eps = ||B.(dR - dRl)||,
            //  with dR = Rphiq - Rphi, dRl = A.dphi and B the diagonal scale
            //  factor matrix.
            //  If eps <= tolerance, the EOA is too conservative and a GROW is
            //  performed,
            //  otherwise, the newly computed mapping is associated to the
            //  initial composition and added to the tree.
            bool checkSolution
            (
                const scalarField& phiq,
                const scalarField& Rphiq
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pChemPointISAT.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
