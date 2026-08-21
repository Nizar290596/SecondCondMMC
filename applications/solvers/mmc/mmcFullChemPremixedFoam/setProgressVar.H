// Progress variable on the DNS field is calculated base of CO2 mass fraction
label meshSize = mesh.C().size();
forAll(YEqvE, specieI) {
    if (YEqvE[specieI].name() == "CO2") {
        for(int i=0; i < meshSize; i++){
	     Xi.Vars("cLES").field()[i] = YEqvE[specieI][i]/YCO2eq;
             if (Xi.Vars("cLES").field()[i] > 1.0) Xi.Vars("cLES").field()[i] = 1.0;
        }
    }
}

Xi.Vars("cLES").field().correctBoundaryConditions();
