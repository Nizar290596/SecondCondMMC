# Code Structure

The code basis is structured in:

 - **applications:** All solvers and utilities 
 - **docs:** Documentation files in Markdown and PDF format
 - **Doxygen:** Doxygen of the code - (optional) coupled with OpenFOAM
 - **etc:** bash environments and bash support functions
 - **python:** Additional python libraries to read and post-process data
 - **src:**
    - **mmcSupport:** Contains all support classes for mmcFoam
    - **mmcChemistryModel:** Particle chemistry models
    - **lagrangian:**
        - **mmc:** Core classes of mmcFoam
        - **spray:** Additional classes for spray coupling with mmcFoam
 - **tests:** unitTests with Catch2
 - **tutorials:** Example case files for different flames and mmcFoam solver
        
