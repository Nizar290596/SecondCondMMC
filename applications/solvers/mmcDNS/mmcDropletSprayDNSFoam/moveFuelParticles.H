Info << "Moving Fuel Particles" << endl;

// Evolve fuel cloud
fuelCloud.evolve();

// Distribute the sources if necessary
fuelCloud.distributeSources();

