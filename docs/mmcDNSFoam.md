# mmcDNSFoam 

The `mmcDNSFoam` solver, part of the `mmcFoam` framework, is designed for coupled Direct Numerical Simulations (DNS) while ensuring independence from stochastic particle models. Unlike standard `mmcFoam` solvers, `mmcDNSFoam` solves the governing equations directly on the Eulerian grid, with stochastic particles treated as passive tracers that represent Large Eddy Simulation (LES) solutions.

## Control Settings

In the system/filterMesh/fvOptions file the filtering of the DNS fields can be switched on or off. Here the default value is on, which means that the DNS fields are filtered to the provided filter mesh and then used as an input for the stochastic particles. 

In the system/filterMesh/fvSolution:
```
PDFMETHOD
{
    filterDNSFields true;   // Switch filtering on or off
}
```

If the filtering is switched on the cloudProperties file and all other required
thermodynamicProperties files etc. must be located in the constant/filterMesh 
directory. 



