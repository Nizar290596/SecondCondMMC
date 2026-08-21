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

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>

\*---------------------------------------------------------------------------*/

#include "samplePlane.H"
#include "sampledPlane.H"
#include "vtkSurfaceWriter.H"
#include <unordered_map>

void Foam::mmcSupport::samplePlane::recombinePointsAndFaces()
{
    // Gather all points and faces
    List<List<face>> allFaces(Pstream::nProcs());
    allFaces[Pstream::myProcNo()] = faces_;
    Pstream::gatherList(allFaces);
    Pstream::scatterList(allFaces);

    List<List<point>> allPoints(Pstream::nProcs());
    allPoints[Pstream::myProcNo()] = points_;
    Pstream::gatherList(allPoints);
    Pstream::scatterList(allPoints);

    // Recombine by finding the right number of points and faces
    DynamicList<face> recombinedFaces;
    DynamicList<point> recombinedPoints;

    forAll(allFaces,procI)
    {
        const auto& faces = allFaces[procI];
        const auto& points = allPoints[procI];
        std::unordered_map<label,label> pointIndexMap;
        for (const auto& f : faces)
        {
            face newFace(f.size());
            // Loop over the points
            forAll(f,i)
            {
                auto it = pointIndexMap.find(f[i]);
                if (it != pointIndexMap.end())
                {
                    newFace[i] = it->second;
                }
                else
                {
                    recombinedPoints.append(points[f[i]]);
                    pointIndexMap.insert
                    (
                        std::pair<label,label>(f[i],recombinedPoints.size()-1)
                    );
                    newFace[i] = recombinedPoints.size()-1;
                }
            }

            recombinedFaces.append(newFace);
        }
    }

    faces_ = recombinedFaces;
    points_ = recombinedPoints;
}


Foam::vector2D 
Foam::mmcSupport::samplePlane::transform(const vector& p, const bool sign) const
{
    // axial component
    const scalar ax = (p & rotAxis_);
    const point pRad = p - rotAxis_*ax;
    scalar r = mag(pRad);

    if (sign && (pRad & rad_) < 0)
        r *= -1.0;

    // third component is always zero 
    return vector2D(ax,r);
}


// * * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * *

Foam::mmcSupport::samplePlane::samplePlane
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    rotAxis_(dict.get<vector>("rotAxis"))
{
    // make sure that rotAxis_ is normalized
    rotAxis_ = rotAxis_/mag(rotAxis_);
    
    bool readData = read();

    if (!readData)
    {
        // Create the sample plane
        autoPtr<sampledPlane> plane(new sampledPlane("plane",mesh, dict));
        plane->update();
        
        faces_ = plane->faces();
        points_= plane->points();
        normal_= plane->normal();

        // recombine points and faces
        recombinePointsAndFaces();
    }

    rad_ = normal_ ^ rotAxis_;
    
    // Normalize rad
    rad_ /= mag(rad_);
    
    coordinates_.resize(faces_.size());
    forAll(coordinates_,i)
    {
        coordinates_[i] = transform(faces_[i].centre(points_),true);
    }
    
    // Set weights to 1.0
    List<scalar> weights(2,1.0);

    // Construct a kdTree to lookup the closest cell quickly
    treePtr_.reset
    (
        new kdTree<vector2D>
        (
            coordinates_,
            weights,
            true,           // median based sorting
            false           // no noramlize with min/max of dimension 
        )
    );
    
    Info << "Sample Plane Information:" <<nl
         << "\trotAxis: "<<rotAxis_<< nl
         << "\tnormal:  "<<normal_ << nl
         << "\trad:     "<<rad_    << endl;

    // Write out faces and points for potential restart
    if (Pstream::master() && !readData)
        write();
    
}

// * * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * *


Foam::label
Foam::mmcSupport::samplePlane::cellIndex
(
    const point& p
) const
{
    // Calculate the radius of the point
    auto projectedPoint = transform(p);
    
    // returns a foundParticle type
    // use idx to access the index
    auto e = treePtr_->nNearest(projectedPoint,1);
    return e[0].idx;
}


void Foam::mmcSupport::samplePlane::writeField
(
    const word& fieldName,
    const Field<scalar>& field
) const
{
    fileName samplePath
    (
        mesh_.time().globalPath()/"mmcStatistics/eulerianStatistics/"
      + mesh_.time().timeName()
    );
    
    samplePath.expand();
    
    fileName filePath = (samplePath/fieldName);
    
    #ifdef FULLDEBUG
    Info << "Write " << fieldName <<" to file "<< filePath<<endl;
    #endif
    Foam::mkDir(samplePath);

    // Write as Field
    OFstream os(filePath);
    os << field << coordinates_;

    surfaceWriters::vtkWriter vtk(points_,faces_,filePath,false);
    vtk.write(fieldName,field);
}


void Foam::mmcSupport::samplePlane::write() const
{
    // Write points, faces and normal to file
    fileName samplePath
    (
        mesh_.time().globalPath()/"mmcStatistics/eulerianStatistics/"
    );
    
    samplePath.expand();
    
    fileName filePath = (samplePath/"samplePlane.dat");

    Foam::mkDir(samplePath);

    OFstream os(filePath);

    os << points_ << faces_ << normal_;
}


bool Foam::mmcSupport::samplePlane::read()
{
    // Write points, faces and normal to file
    fileName samplePath
    (
        mesh_.time().globalPath()/"mmcStatistics/eulerianStatistics/"
    );
    
    samplePath.expand();
    
    fileName filePath = (samplePath/"samplePlane.dat");

    if (Foam::exists(filePath))
    {
        IFstream is(filePath);
        is >> points_ >> faces_ >> normal_;
        return true;
    }
    return false;
}


// ************************************************************************* //
