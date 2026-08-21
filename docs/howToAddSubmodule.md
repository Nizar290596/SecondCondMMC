# Adding a New Submodule

To add a new sub-model, it is recommended to start from an existing sub-model, e.g., the MMCCurl sub-model, copy and rename the file as a starting point for the new sub-model. Once the new sub-model is defined the macro to generate this sub-model has to be extended. These macros are found in the `include/` folder in the `popeParticles` directory. Their notation follows following syntax, `make<ParticleType><SubModelName>Models.H`, e.g., for the mixing model it is, `makeMixingPopeParticleMixingModels.H`. However, this does not mean that the sub-model is only working with the in the name specified particle type, but shows on which particle properties it depends. 
If the sub-model has been extending exiting models the process ends here, however, if a completely new sub-model class is added the calls of the macros in the `basic<ParticleName>/makeBasic<ParticleName>Submodels.C` has to be adjusted. 

The general workflow to add a new sub-model is then:

 - Write new sub-model class in the submodels directory. Preferable by copying and renaming existing models to keep a coherent structure. 
 - Add the new sub-model to the respective `make<ParticleType><SubModelName>Models.H` in the `include/` directory of the popeParticles.
 - If it is a new sub-model class that did not exist before also adjust the calls of the macros in the `basic<ParticleName>/makeBasic<ParticleName>Submodels.C` located in the `popeParticles/derived` directory.  
