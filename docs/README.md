# mmcFoam 

mmcFoam is a turbulent combustion code for OpenFOAM. It was originally
developed in 2013 as a turbulent non-premixed combustion code by Yipeng Ge and
atthew Cleary. It has since evolved into a collaborative code development and
sharing project. Participing institutions are the University of Sydney,
University of Queensland, University of Stuttgart, University of New South
Wales, Karlsruhe Institute of Technology, Indian Institute
of Technology (Kanpur), Physikalisch-Technische Bundesanstalt
(Braunschweig) and National University of Singapore.

From 2013 to 2016 a number of updates were added including options for sparse
and dense mixing models, alternative chemical integrators and efficient
particle mixing pair selection algorithms. In May 2017 mmcFoam-4.x was released
following a major restructure by Sebastian Galindo-Lopez which generalised the
code expanding its capability beyond non-premixed combustion. At the same time
additional code blocks contributed from the University of Sydney and the
University of Stuttgart were added for multiphase combustion, premixed and
stratified combustion and particle nucleation and growth. mmcFoam-5.x was
released in May 2018 with major input from Zhijie Huo and Sebastian
Galindo-Lopez at the University of Sydney and Eshan Sharma at IIT Kanpur which
included a restructure of the spray code.

In 2023 the mmcFoam code basis was changed to the ESI-OpenCFD branch of
OpenFOAM with the release of mmcFoam-v2312. This version is also compatible
with OpenFOAM v2306 and v2212. 

## CONTRIBUTORS TO CODE DEVELOPMENT:

(Current) Main developers of the current mmc code basis are:

 * Jan Gärtner (University of Stuttgart) - porting to ESI-OpenCFD branch v2306
 * Matthew Cleary (University of Sydney)
 * Sebastian Galindo (University of Sydney)
 

mmcFoam is a team effort and all contributors can be found in [Contributors List](./docs/contributors.md) 

## PUBLISHING WITH mmcFoam:

The mmcFoam collaboration is designed to support a collegiate research
community and to foster sharing of state-of-the-art developments from all
groups. If you publish mmcFoam simulation results please cite the following key
paper:
<cite>
S. Galindo-Lopez, F. Salehi, M.J. Cleary, A.R. Masri, G. Neuber, O.T.
Stein, A. Kronenburg, A. Varna, E.R. Hawkes, B. Sundaram, A.Y. Klimenko and Y.
Ge, A stochastic  multiple  mapping conditioning  computational  model  in
OpenFOAM for turbulent combustion, Computers and Fluids (2018),
<a href="https://doi.org/10.1016/j.compﬂuid.2018.03.083">https://doi.org/10.1016/j.compﬂuid.2018.03.083</a>
</cite>

A full list of mmcFoam publications may be found in the publication [log](./docs/publication_log.md).

From time to time source code and solvers will be released prior to
publication. A list of such solvers is listed below. You are welcome to use
these but please refrain from publishing them yourself before the main
developers have had that opportunity. Please contact Matthew Cleary
(m.cleary@sydney.edu.au) if you intend to use them.  
Solvers under development:

 * All solvers published and up to date as of 21 January 2025


## DOWNLOADING AND INSTALLING THE SOURCE:

This section describes how to download and compile the very latest source code of
mmcFoam-v2406. The code is continually updated and it is recommended that users
regularly pull the latest code on the master branch of the git repository.

To distinguish different mmcFoam releases and to use different release versions
simultaneously, mmcFoam appends the current git branch or release to the solver
name. E.g. the solver `mmcFoam` compiles with the name `mmcFoam-v2406` if the 
version 2406 is checked out. Solvers on the master branch do not have any 
appendix to their name!

Download,install and compile mmcFoam as follows.

1. `git clone --recursive git@bitbucket.org:mcleary/mmcfoam-esi.git mmcFoam` 
   or alternative 
   `git clone --recursive https://github.tik.uni-stuttgart.de/ITV/mmcFoam mmcFoam`

2. (Optional) if the `--recursive` flag was not used during cloning execute:
    * `cd mmcFoam` (if not already in the mmcFoam folder)
    * `git submodule update --init --recursive`

3. `cd mmcFoam`

4. `./Allwmake` or use `./Allwmake -j` to compile on multiple cores.
   Use `./Allwmake -h` for more information.

Notes: The environment variables of mmcFoam are stored in the etc/ folder of
mmcFoam. To source the environment variables execute in a shell:
```bash
source <path_to_mmcFoam>/mmcFoam/etc/bashrc
```
This file has to be sourced if parts of mmcFoam are compiled using `wmake` and not the `Allwmake` script.

## CHECK OUT SPECIFIC VERSION

mmcFoam is using a similar branching model for the different OpenFOAM versions 
as openfoam.com does. Hence, each mmcFoam release for an OpenFOAM version has 
its own branch, labeled with the version of OpenFOAM. For example, branch v2306
 contains the development for the mmcFoam-v2306 version. 
However, as bug fixes may occur, the released branches are not static. To allow
an easy reference to a release version, a tagged version is created with each 
release, titled, e.g., mmcFoam-v2306. Therefore, for each new mmcFoam release, 
a corresponding branch with, e.g., v2306 and a corresponding named 
tag mmcFoam-v2306 exists. Further, the annotated tag message contains the 
release notes of this mmcFoam release, summarizing the most important changes. 


## USE MULTIPLE VERSIONS SIMULTANEOUSLY

mmcFoam solvers are compiled with an extension if they are on a git branch or
checked out tag. E.g. compiling mmcFoam with the v2012 branch/tag compiles
mmcFoam into the name mmcFoam-v2012. If it is compiled on the master branch no
extension is given. 


## GETTING STARTED:

Execute a tutorial file of mmcFoam in your run/ folder:

1. `cd $WM_PROJECT_USER_DIR`

2. `mkdir run`

3. `cd run`

4. Source the mmcFoam environment variables with
  ```bash
source <path_to_mmcFoam>/mmcFoam/etc/bashrc
  ```

5. `cp -r $MMC_PROJECT_USER_DIR/tutorials/ mmcFoamTutorials`

6. `cd mmcFoamTutorials/oneStepFlameSheet`

7. `./Allrun`

## Post Processing with Python

It is possible to collect particle and eulerian statistics with the built in 
tools of mmcFoam. To post-process the data a python library is provided to 
make reading and processing the data easier.

For further information about the python library see: [README-Python](./python/README.md)

---

## DOCUMENTATION

mmcFoam is documented on several levels:

1. Code documentation via Doxygen:
    * Go to `Doxygen/` directory of mmcFoam and execute `makeDoxygen` to generate
     a Doxygen manual of the code basis
    * Access the code documentation with your browser by opening the `mainPage.html` file
2. Documentation of some solvers, utilities, and libraries as a presentation or 
   text document in PDF. These files are stored in [pdfDocumentations](./docs/pdfDocumentations)
3. A more detailed report of the upgrade to mmcFoam-v2306 can be found [here](./docs/pdfDocumentations/mmcFoam-v2306-FullDocumentation.pdf) 
4. mmcFoam is also documented in the published literature. See the 
   [publication](./publication_log.md) file for more information.
5. Work in progress: A mkdocs web page can be created -- but unfinished and not documented

**Note:** An example `cloudProperties` file can be found in [exampleCloudProperties](./docs/exampleCloudProperties.md).

---

## PUBLICATION

Published work with mmcFoam can be found in this list [here](./docs/publication_log.md)

## LICENSE

mmcFoam is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.  See the file COPYING in this directory or
http://www.gnu.org/licenses/, for a
description of the GNU General Public License terms under which you
may redistribute files.
