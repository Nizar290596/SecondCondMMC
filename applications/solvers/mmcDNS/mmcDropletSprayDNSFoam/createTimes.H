Foam::Info<< "Create DNS time\n" << Foam::endl;
 
Foam::Time runTime(Foam::Time::controlDictName, args);

// It is important to create a new time object to construct a second,
// independent object registry for the filtered field. Otherwise a lot of 
// the internal lookup functions of the mmc cloud may not work properly!

Foam::Info<< "Create filter time\n" << Foam::endl;
// Create an own time and object registry for the filter mesh
Time filterTime(Time::controlDictName,args);    
