Info << "Moving Fuel Particles" << endl;

// Evolve fuel cloud
fuelCloud.evolve();

// Gather particle statistics
pStatsFuelCloud.particleSampling(fuelCloud);

