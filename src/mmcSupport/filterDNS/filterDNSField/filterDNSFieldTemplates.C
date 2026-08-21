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


template<class fieldType>
Foam::tmp<fieldType>
Foam::filterDNSField::mean
(
    const fieldType& field
) const
{
    // Create a density field 
    volScalarField rho
    (
        IOobject
        (
            "filteredField_rhoTMP",
            mesh().time().timeName(),
            mesh().time(),
            IOobject::NO_READ
        ),
        this->mesh(),
        dimensioned<scalar>("0",dimDensity,1.0)
    );

    return FavreAverage(field,rho);
}


template<class fieldType>
Foam::tmp<fieldType>
Foam::filterDNSField::FavreAverage
(
    const fieldType& field,
    const volScalarField& rho
) const
{
    typedef typename fieldType::value_type Type;

    // Create a new temporary field
    tmp<fieldType> filteredField
    (
        new fieldType
        (
            IOobject
            (
                "filteredField_"+field.name(),
                filterMesh().time().timeName(),
                filterMesh().time(),
                IOobject::NO_READ
            ),
            this->filterMesh(),
            dimensioned<Type>("0",field.dimensions(),pTraits<Type>::zero)
        )
    );

    fieldType& filteredFieldRef = filteredField.ref();

    updateField(field,rho,filteredFieldRef);
    
    return filteredField;
}


template<class fieldType>
void Foam::filterDNSField::updateField
(
    const fieldType& field,
    fieldType& filteredField
) const
{
    // Create a density field 
    volScalarField rho
    (
        IOobject
        (
            "filteredField_rhoTMP",
            mesh().time().timeName(),
            mesh().time(),
            IOobject::NO_READ
        ),
        this->mesh(),
        dimensioned<scalar>("0",dimDensity,1.0)
    );

    return updateField(field,rho,filteredField);
}


template<class fieldType>
void Foam::filterDNSField::updateField
(
    const fieldType& field,
    const volScalarField& rho,
    fieldType& filteredField
) const
{
    typedef typename fieldType::cmptType cmptType;

    typedef typename fieldType::value_type Type;

    // Create a new field for the filtered density field
    volScalarField rhoFiltered
    (
        IOobject
        (
            "filteredField_rho",
            filterMesh().time().timeName(),
            filterMesh().time(),
            IOobject::NO_READ
        ),
        this->filterMesh(),
        dimensioned<scalar>("0",dimDensity,0)
    );


    // Initialize the filteredField to zero
    forAll(filteredField,j)
    {
        filteredField[j] = pTraits<Type>::zero;
    }

    
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);


    // In the buffer we stream the remote cells that have to be send
    // The information is coded as [filterCellID, [Value, volFraction*rho]]
    // The second entry is the weighted density
    //
    // Unfortunately we have to use lists here and cannot directly stream to the
    // buffer due to an unkown readRaw error which appears with the Tuple2 type.
    typedef Tuple2<label,Tuple2<Type,scalar>> cellInfoType;
    List<DynamicList<cellInfoType>> toSend(Pstream::nProcs());

    // Loop over the DNS field for local cells
    // =======================================
    forAll(field,celli)
    {
        // Loop over all local cells
        forAll(cellToFilterCellProcessorID_[celli],j)
        {
            // Check if it is a local processor
            if (cellToFilterCellProcessorID_[celli][j] == Pstream::myProcNo())
            {
                for (direction cmpt=0; cmpt < pTraits<Type>::nComponents; ++cmpt)
                {
                    cmptType& val = setComponent(filteredField[cellToFilterCell_[celli][j]],cmpt);
                    val += cellToFilterCellVolRatio_[celli][j] * rho[celli]
                      * component(field[celli],cmpt);
                }
                rhoFiltered[cellToFilterCell_[celli][j]] += 
                    cellToFilterCellVolRatio_[celli][j]*rho[celli];
            }
            else
            {
                const label proci = cellToFilterCellProcessorID_[celli][j];
                cellInfoType cellInfo;
                cellInfo.first() = cellToFilterCell_[celli][j];
                cellInfo.second().first() = field[celli];
                cellInfo.second().second() = cellToFilterCellVolRatio_[celli][j]*rho[celli];
                toSend[proci].append(cellInfo);
            }
        }
    }
    
    // Stream list to buffer
    forAll(toSend,procI)
    {
        if (procI == Pstream::myProcNo())
            continue;
        UOPstream pStream(procI,pBufs);
        pStream << toSend[procI];
    }

    // Correct parallel run
    // ====================
    pBufs.finishedSends();

    for (label procI=0; procI < Pstream::nProcs(); procI++)
    {
        if (procI != Pstream::myProcNo())
        {
            UIPstream pStream(procI, pBufs);
            DynamicList<cellInfoType> receivedCellInfo;
            pStream >> receivedCellInfo;
            for (auto& cellInfo : receivedCellInfo)
            {
                // Get the local filter cell ID
                const label filterCellID = cellInfo.first();
                const Type& receivedValue = cellInfo.second().first();
                const scalar& receivedScaledRho = cellInfo.second().second();

                for (direction cmpt=0; cmpt < pTraits<Type>::nComponents; ++cmpt)
                {
                    cmptType& val = setComponent(filteredField[filterCellID],cmpt);
                    val += receivedScaledRho*component(receivedValue,cmpt);
                }
                rhoFiltered[filterCellID] += receivedScaledRho;
            }
        }
    }

    forAll(filteredField,cellj)
    {
        if (nCellsInFilterCell_[cellj] != 0)
            filteredField[cellj] /= rhoFiltered[cellj];
    }

    // Update the boundary
    updateBoundary(filteredField);
    filteredField.correctBoundaryConditions();
}


template<class fieldType>
void Foam::filterDNSField::updateBoundary
(
    fieldType& field
) const
{
    
    auto& bField = field.boundaryFieldRef();
    
    // HPC PATCH BOUNDARY
    // Construct the patch internal boundary field manually due to an OpenFOAM
    // bug documented in:
    // https://develop.openfoam.com/Development/openfoam/-/issues/3277
    auto bInternalFieldTmp = constructInternalBoundaryField(bField);
    

    // Get const access to the boundary internal field
    auto& bInternalField = bInternalFieldTmp();

    // Loop over all patches
    forAll(bField,patchi)
    {
        auto& bPatchField = bField[patchi];
        auto& bPatchInternalField = bInternalField[patchi];
        // Loop over all entries in the patch field
        forAll(bPatchField,i)
        {
            bPatchField[i] = bPatchInternalField[i];
        }
    }
}


template<class bFieldType>
Foam::tmp<bFieldType> Foam::filterDNSField::constructInternalBoundaryField
(
    bFieldType& bField
) const
{
    auto tresult = tmp<bFieldType>::New
    (
        bFieldType::Internal::null(),
        bField
    );

    auto& result = tresult.ref();

    forAll(result, patchi)
    {
        result[patchi] == bField[patchi].patchInternalField();
    }

    return tresult;
}


template<class fieldType>
Foam::tmp<fieldType> 
Foam::filterDNSField::mapFilterFieldToDNS
(
    const fieldType& filteredField
) const
{
    // Create a DNS field 
    typedef typename fieldType::value_type Type;

    // Create a new temporary field
    tmp<fieldType> fieldRef
    (
        new fieldType
        (
            IOobject
            (
                "mappedField_"+filteredField.name(),
                mesh().time().timeName(),
                mesh().time(),
                IOobject::NO_READ
            ),
            this->mesh(),
            dimensioned<Type>("0",filteredField.dimensions(),pTraits<Type>::zero)
        )
    );
    fieldType& field = fieldRef.ref();

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Store the information of the filter cells to send to other processors
    // Stored as a Tuple of the cellID of the DNS mesh on the remote processor
    // and the filter cell value:
    // first():  filter cell value
    // second().first():  cellID of the DNS cell on the remote processor
    // second().second(): volume overlap of the DNS cell in the filter cell
    List<DynamicList<Tuple2<Type,Tuple2<label,scalar>>>> filterCellsToSendToProcList(Pstream::nProcs());

    // Loop over the filteredField and map all local ones
    forAll(filteredField,i)
    {
        // Get all DNS cells
        const List<label>& cellID = filterCellToCell_[i];
        forAll(cellID,j)
        {
            field[cellID[j]] = filteredField[i];
        }

        // Insert all requested cells from the given processor
        const auto& remoteCells = filterCellToRemoteCells_[i];
        forAll(remoteCells,j)
        {
            const label DNSCellID = remoteCells[j].first();
            const label toProc = remoteCells[j].second().first();
            const scalar volOverlap = remoteCells[j].second().second();

            auto& filterCellsToSendToProc = filterCellsToSendToProcList[toProc];
            filterCellsToSendToProc.append
            (
                Tuple2<Type,Tuple2<label,scalar>>
                (
                    filteredField[i],
                    Tuple2<label,scalar>
                    (
                        DNSCellID,
                        volOverlap
                    )
                )
            );
        }
    }

    // Stream list to buffer
    forAll(filterCellsToSendToProcList,procI)
    {
        if (procI == Pstream::myProcNo())
            continue;

        UOPstream pStream(procI,pBufs);
        pStream << filterCellsToSendToProcList[procI];
    }

    // Correct parallel run
    // ====================
    pBufs.finishedSends();

    for (label procI=0; procI < Pstream::nProcs(); procI++)
    {
        if (procI == Pstream::myProcNo())
            continue;

        UIPstream pStream(procI,pBufs);
        DynamicList<Tuple2<Type,Tuple2<label,scalar>>> receivedFilteredCells;
        pStream >> receivedFilteredCells;
        for (auto& tuple : receivedFilteredCells)
        {
            field[tuple.second().first()] = tuple.second().second()*tuple.first();
        }
    
    
    }
    // Update the boundary
    updateBoundary(field);
    field.correctBoundaryConditions();

    updateBoundary(field);

    return fieldRef;
}


template<class fieldType>
Foam::tmp<fieldType> 
Foam::filterDNSField::mapSourceFilterFieldToDNS
(
    const fieldType& filteredField
) const
{
    // Create a DNS field 
    typedef typename fieldType::value_type Type;

    // Create a new temporary field
    tmp<fieldType> fieldRef
    (
        new fieldType
        (
            IOobject
            (
                "mappedSourceField_"+filteredField.name(),
                mesh().time().timeName(),
                mesh().time(),
                IOobject::NO_READ
            ),
            this->mesh(),
            dimensioned<Type>("0",filteredField.dimensions(),pTraits<Type>::zero)
        )
    );
    fieldType& field = fieldRef.ref();

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Store the information of the filter cells to send to other processors
    // Stored as a Tuple of the cellID of the DNS mesh on the remote processor
    // and the filter cell value:
    // first():  filter cell value
    // second().first():  cellID of the DNS cell on the remote processor
    // second().second(): volume overlap of the DNS cell in the filter cell
    List<DynamicList<Tuple2<Type,Tuple2<label,scalar>>>> filterCellsToSendToProcList(Pstream::nProcs());

    // Loop over the filteredField and map all local ones
    forAll(filteredField,i)
    {
        // Get all DNS cells
        const List<label>& cellID = filterCellToCell_[i];
        forAll(cellID,j)
        {
            field[cellID[j]] = filteredField[i]/nCellsInFilterCell_[i];
        }

        // Insert all requested cells from the given processor
        const auto& remoteCells = filterCellToRemoteCells_[i];
        forAll(remoteCells,j)
        {
            const label DNSCellID = remoteCells[j].first();
            const label toProc = remoteCells[j].second().first();
            const scalar volOverlap = remoteCells[j].second().second();

            auto& filterCellsToSendToProc = filterCellsToSendToProcList[toProc];
            filterCellsToSendToProc.append
            (
                Tuple2<Type,Tuple2<label,scalar>>
                (
                    filteredField[i]/nCellsInFilterCell_[i],
                    Tuple2<label,scalar>
                    (
                        DNSCellID,
                        volOverlap
                    )
                )
            );
        }
    }

    // Stream list to buffer
    forAll(filterCellsToSendToProcList,procI)
    {
        if (procI == Pstream::myProcNo())
            continue;

        UOPstream pStream(procI,pBufs);
        pStream << filterCellsToSendToProcList[procI];
    }

    // Correct parallel run
    // ====================
    pBufs.finishedSends();

    for (label procI=0; procI < Pstream::nProcs(); procI++)
    {
        if (procI == Pstream::myProcNo())
            continue;

        UIPstream pStream(procI,pBufs);
        DynamicList<Tuple2<Type,Tuple2<label,scalar>>> receivedFilteredCells;
        pStream >> receivedFilteredCells;
        for (auto& tuple : receivedFilteredCells)
        {
            field[tuple.second().first()] = tuple.second().second()*tuple.first();
        }
    }

    return fieldRef;
}

// ************************************************************************* //

