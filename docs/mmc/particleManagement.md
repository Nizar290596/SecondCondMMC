# Particle Management

Active particle management is often required to ensure a sufficient density of 
stochastic mmc particles in all areas of the simulation domain. This is due to 
particles within the hot reaction zone moving much faster and further than 
particles in the cold areas. Hence, without particle management, the particle 
density within the reaction zones would quickly decrease and potentially be 
insufficient to represent the underlying statistics correctly. Therefore, 
active particle management is required to ensure enough particles exist within 
the reaction zone. Determining a smooth particle density field on the fly is 
difficult, so the mesh is divided into control zones. A defined particle number
range is set and actively controlled in each zone or cell. As these zones 
represent a second more coarse mesh than the underlying LES mesh, it is called
 **superMesh** in the code.  

## Usage

The particle management is controlled over a dictionary in the 
`constant/cloudProperties` file. The dictionary has to be present, even if the 
particle management is switched off.

```C++
particleManagement
{
    numCtlOn          false; // switch management on or off
    Nlo                  30; // minimum value of particles
    Npc                  50; // target value of particles in each super cell
    Nhi                  80; // maximum value of particles
}
```

If the particle management is switched on it will start adding particles once
less than `Nlo` or more than `Nhi` particles are present within one super cell. 
If too many particles are present the n particles with the lowest weight are 
deleted, such that the target value `Npc` is reached. E.g. in the given example
above, if 81 particles are present in the super cell, the 31 particles with the 
lowest weight will be deleted. To avoid a bias in the process, double the number
of particles is selected and then deleted with a 50% chance. In the exmaple this
would mean, selecting 61 particles and deleting 30 with the lowest weight.


## Implementation

The particle management is implemented in the particleNumberController class
located in the mmc, clouds, ItoPopeCloud folder. This class also is responsible 
for managing the superMesh and all related operations. 

### Parallel Scalability

In previous mmcFoam versions, the superMesh had to be decomposed with the 
general LES mesh, effectively limiting the maximum decomposition of the domain 
to the cell number of the superMesh. To solve this limitation, the current 
particleNumberController always reads the complete superMesh from the case 
constant folder. This means that all processors know the entire superMesh, 
which, due to the typically small size of the superMesh is not a significant 
memory or performance factor. The particleNumberController class then collects 
the information to which processors the super cells are distributed and stores 
the information in a bi-directional list, allowing the association of LES cells 
to super cell and super cell to LES cells. Further details about the 
implementation are found in the code comments of the class itself. 


