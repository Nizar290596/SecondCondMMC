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

#include "particleNumberController.H"

// * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * * * 

Foam::particleNumberController::particleNumberController
(
    const fvMesh& mesh,
    const dictionary& particleManagementDict
)
:
    mesh_(mesh),
    particleManagementDict_(particleManagementDict),
    numCtrlOn_(particleManagementDict_.lookupOrDefault<Switch>("numCtlOn", true)),
    Npc_(readScalar(particleManagementDict_.lookup("Npc"))),
    Nlo_(readScalar(particleManagementDict_.lookup("Nlo"))),
    Nhi_(readScalar(particleManagementDict_.lookup("Nhi"))),
    superMesh_
    (
        constructSuperMesh(mesh)
    )
{
    nSuperCells_ = superMesh_->cells().size();
    
    // Set sizes
    superCellForCell_.resize(mesh.cells().size());
    
    cellsInSuperCell_.resize(nSuperCells_);
    
    // Reserve the space for each super cell
    forAll(cellsInSuperCell_,superCellI)
    {
        cellsInSuperCell_[superCellI].reserve(mesh.nCells()/nSuperCells_);
    }
    
    // Counter of how many cells could not be found with findCell()
    label numCellsNotFoundInSuperMesh = 0;
    
    forAll(mesh.cells(),celli)
    {
        point p = mesh.cellCentres()[celli];
        
        label superCellI = superMesh_->findCell(p);
        
        if (superCellI == -1)
        {
            superCellI = superMesh_->findNearestCell(p);
            if (numCellsNotFoundInSuperMesh < 20)
                WarningInFunction() << "Cell "<<p
                    <<" not found in super cell with fvMesh::findCell()"<<endl;
            else if (numCellsNotFoundInSuperMesh == 20)
                WarningInFunction() << "Warning messages are truncated"<<endl;
        }
        
        if (superCellI == -1)
            FatalError 
                << "Cell "<<p<<" not found in super mesh"
                <<exit(FatalError);

        superCellForCell_[celli] = superCellI;
        
        cellsInSuperCell_[superCellI].append(celli);
    }

    constructSuperCellToProcID();

    Info << "Particle Number Controller Settings:" << nl
         << "\tactive: " << numCtrlOn_ << nl
         << "\tNlo:    " << this->Nlo_ << nl
         << "\tNpc:    " << this->Npc_ << nl
         << "\tNhi:    " << this->Nhi_ << nl << endl;
};

// * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * * * 

Foam::fvMesh* 
Foam::particleNumberController::constructSuperMesh
(
    const fvMesh& mesh
)
{
    Info << "Read superMesh..."<<endl;
    // First generate a new Time instance to set the correct 
    // dbDir() path!
    time_.reset(
        new Time
        (
            "controlDict",
            mesh.time().rootPath(),
            mesh.time().globalCaseName(),
            false,
            false
        )
    );
    
    fvMesh* meshPtr = new fvMesh
    (
        IOobject
        (
            "super",
            time_->caseConstant(),
            time_(),
            IOobject::MUST_READ
        )
    );
    
    return meshPtr;
}


void Foam::particleNumberController::constructSuperCellToProcID()
{
    nSuperCellsOnProcessor_ = 0;
    
    if (Pstream::parRun())
    {
        // each processor creates a list of super cells included in its local mesh
        superCellOnProcessor_.resize(nSuperCells_,false);
        superCellToProcID_.resize(nSuperCells_);
        
        forAll(cellsInSuperCell_,superCellI)
        {
            if (cellsInSuperCell_[superCellI].size() > 0)
            {
                superCellOnProcessor_[superCellI] = true;
                nSuperCellsOnProcessor_++;
            }
        }
        
        List<List<bool>> superCellsInProc_AllProcessor(Pstream::nProcs());
        superCellsInProc_AllProcessor[Pstream::myProcNo()] = superCellOnProcessor_;
        
        // Gather and scatter List
        Pstream::gatherList(superCellsInProc_AllProcessor);
        Pstream::scatterList(superCellsInProc_AllProcessor);

        forAll(superCellsInProc_AllProcessor,procI)
        {
            const List<bool>& superCellInProc = superCellsInProc_AllProcessor[procI];
            forAll(superCellInProc,superCellI)
            {
                if (superCellInProc[superCellI])
                    superCellToProcID_[superCellI].append(procI);
            }
        }
    }
    else
    {
        superCellOnProcessor_.resize(nSuperCells_,true);
        nSuperCellsOnProcessor_=nSuperCells_;
    }
}


Foam::List<Foam::DynamicList<Foam::point>>
Foam::particleNumberController::getCellPositionsInSuperCell()
{    
    return fieldDataInSuperCell(mesh_.C());
}



Foam::List<Foam::particleNumberController::pContainer> 
Foam::particleNumberController::initParticlePosBasedOnField
(
    const volScalarField& rho
)
{
    // =========================================================================
    // The ItoPopeParticles are initialized based on their mass, 
    // e.g. the total particle mass should represent the mass in the domain
    // =========================================================================
    
    // First collect the mass value in each super cell
    volScalarField::Internal mass
    (
        "Mass",
        mesh_.V() * rho
    );
    
    List<DynamicList<scalar>> massInSuperCell = fieldDataInSuperCell(mass);
    
    // The position of cells in each super cell
    List<DynamicList<point>> cellPositionsInSuperCell = getCellPositionsInSuperCell();
    
    //- Calculate the sum of mass in each super cell
    List<scalar> superCellMass(nSuperCells_,0.0);
    
    forAll(massInSuperCell, superCellI)
    {
        for (scalar& m : massInSuperCell[superCellI])
            superCellMass[superCellI] += m;
    }

    // Populate each superCell with Npc Pope particles at random locations with mass m.

    DynamicList<pContainer> listParticlePosMass(Npc_*nSuperCellsOnProcessor());

    // Note to ranom number generator
    // All processors should return the same results for the same data set,
    // --> All processors use the same seed for their random number generator
    List<int> seedsForSuperCell = synchronizeRandomSeed();

    
    for (label superCellI = 0; superCellI < nSuperCells_; ++superCellI)
    {
        // if the super cell is not on the processor continue
        if (!superCellOnProcessor_[superCellI] || superCellMass[superCellI] == 0)
            continue;
        
        // get random points in the supercell
        List<point> points = randomPointsInSuperCell
            (
                superCellI,
                Npc_,
                seedsForSuperCell[superCellI]
            );
        
        StochasticLib1 rndGen(seedsForSuperCell[superCellI]);
        
        DynamicList<pContainer> newParticles;

        // Loop over all random points in the super cell and check if they 
        // are on this processor
        for (point& p : points)
        {
            label icell = -1;
            label itetface = -1;
            label itetpt = -1;
                
            mesh_.findCellFacePt(p,icell,itetface,itetpt);
  
            // Check if point is in the local mesh
            if (icell == -1)
                continue;
           
            // Create the tetIndices for this point
            tetIndices tetInd(icell,itetface,itetpt);
           
           
            // Convert the point into barycentric coordinates
            barycentric b = transformPoint(p,tetInd.tet(mesh_));
            pContainer newParticle;

            newParticle.mass     = 0.0;
            newParticle.location = b;
            newParticle.tetFace  = itetface;
            newParticle.tetPt    = itetpt;
            newParticle.cell     = icell;
             
            newParticles.append(newParticle);
        }
        
        
        //- Create list of particle masses and points
        for (pContainer& p : newParticles)
            p.mass = superCellMass[superCellI] / points.size();
        
        listParticlePosMass.append(newParticles);
    }

    return listParticlePosMass;
}


Foam::List<int>
Foam::particleNumberController::synchronizeRandomSeed() const
{
    List<int> seedsForSuperCell(nSuperCells_,0);
    
    if (Pstream::myProcNo() == 0)
    {
        // The root processor generates the seed and scatters it to all 
        // processors
        StochasticLib1 rndGen(time(0));
        forAll(seedsForSuperCell,superCellI)
        {
            seedsForSuperCell[superCellI] = rndGen.IRandom(-247483640,247483640);
        }
    }
    
    Pstream::broadcast(seedsForSuperCell);

    return seedsForSuperCell;
}


Foam::List<Foam::point>
Foam::particleNumberController::randomPointsInSuperCell
(
    const label superCellI,
    const label nPoints,
    const int seed
) const
{
    // Decompose the super mesh into tetrahedar
    List<tetIndices> superCellTets = polyMeshTetDecomposition::cellTetIndices
    (
        superMesh(),
        superCellI
    );
    
    const scalar superCellVolume = superMesh().V()[superCellI];
    
    Random ofRand(seed);
    
    DynamicList<point> points;
    points.reserve(nPoints);
    
    forAll(superCellTets,superTetI)
    {
        const tetIndices& superCellTetIs = superCellTets[superTetI];

        tetPointRef superTet = superCellTetIs.tet(superMesh());

        // Number of particles inserted per tet is proportional to volume.
        // Assumes constant density in superCell.
        scalar superTetVolume = superTet.mag();
        
        label nPointsInSuperTet = std::ceil
            (superTetVolume / superCellVolume * nPoints);
        
        // Randomly scatter the points
        for (label k = 0; k < nPointsInSuperTet; k++)
        {
            barycentric b = barycentric01(ofRand);
            points.append(superTet.barycentricToPoint(b));
        }
    }

    // if size of points is significantly larger than nPoints,
    // shuffle the list to avoid initializing particles always in the frist 
    // tetrahedars of the cell
    if (points.size() > 1.2*nPoints)
    {
        ofRand.shuffle(points);
    }

    points.resize(nPoints);
    
    return points;
}


Foam::barycentric
Foam::particleNumberController::transformPoint
(
    const point& p, 
    const tetPointRef tet
)
{
    // Algorithm taken from: 
    // https://www.cdsimpson.net/2014/10/barycentric-coordinates.html
    
    // get the three vectors of the tet
    const point& a = tet.a();
    const point& b = tet.b();
    const point& c = tet.c();
    const point& d = tet.d();
    
    point vap = p - a;
    point vbp = p - b;

    point vab = b - a;
    point vac = c - a;
    point vad = d - a;

    point vbc = c - b;
    point vbd = d - b;
    // ScTP computes the scalar triple product
    auto ScTP = [](point& a, point& b, point& c) -> double
    {
        return a & (b^c);
    };
    
    
    double va6 = ScTP(vbp, vbd, vbc);
    double vb6 = ScTP(vap, vac, vad);
    double vc6 = ScTP(vap, vad, vab);
    double vd6 = ScTP(vap, vab, vac);
    double v6 = 1.0 / ScTP(vab, vac, vad);
    return barycentric(va6*v6, vb6*v6, vc6*v6, vd6*v6);
}

// * * * * * * * * * * *  Public Member Functions * * * * * * * * * * * * * * * 

Foam::label 
Foam::particleNumberController::getSuperCellID(const label cellID) const
{
    return superCellForCell_[cellID];
}
