# Add a Spray Parcel

A new thermo parcel for sprays is added. 

 - derived: Typename definition of the parcels and calling the make
   pre-processor macros of the include/ directory
 - include: Pre-processor macros to create the different sub models for the
   clouds  
 - Templates: The template class, sprayThermoParcel, is defined

## Adding a new parcel

1. Add a new class to the Templates/ 
2. Create a new basic directory in derived/ and create the typedefs as well as
   the make folder
3. Done

## Adding a new submodule

In the derived/make<Parcel> file the pre processor macros are called for the
submodels. These macros are defined in the spray/parcel/include/ folder. In
this folder new sub-model macros have to be added, e.g. for makeParcelForces.H
add makeParticleForceModelType(VirtualMassForce, CloudType); to add the 
VirtualMassForce sub model 
