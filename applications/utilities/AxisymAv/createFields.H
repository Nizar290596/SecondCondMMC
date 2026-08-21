    Info << "Creating fields" << nl << endl;


    IOdictionary PostProcDict
    (
        IOobject
        (
            "PostProcess",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    //find field names need circumferential average
    dictionary AxisymAvDict
    (
        PostProcDict.subDict("AxisymAvDict")
    );

    List<word> AxisymAvFieldNames = AxisymAvDict.toc();

    Info << "Fields for circumferential average: " << AxisymAvFieldNames << endl;

    PtrList<volScalarField> AvFields;

    PtrList<volScalarField> CircAvFields;

    forAll(AxisymAvFieldNames,Yi)
    {

        Info << "Reading AvField " << AxisymAvFieldNames[Yi] << nl << endl;

        volScalarField* Av;
            
        Av = new volScalarField
        (
            IOobject
            (
                AxisymAvFieldNames[Yi],
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        AvFields.append(Av);


        Info << "creating CircAvField for " << AxisymAvFieldNames[Yi] << nl << endl;

        volScalarField* AvCirc;
            
        AvCirc = new volScalarField
        (
            IOobject
            (
                AxisymAvFieldNames[Yi]+"CircAv",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            AvFields[Yi]
        );

        CircAvFields.append(AvCirc);

    }
/*
    volScalarField TAv
    (
        IOobject
        (
            "TAv",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );


    volScalarField TAvCircAv
    (
        IOobject
        (
            "TAvCircAv",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        TAv*0
    );
*/



//    Info<< "TAvCircAv.internalField()[100] = " << TAvCircAv.internalField()[100] << endl;
//    Info<< "TAvCircAv.internalField()[1000] = " << TAvCircAv.internalField()[1000] << endl;
//    Info<< "TAvCircAv.internalField()[10000] = " << TAvCircAv.internalField()[10000] << endl;




