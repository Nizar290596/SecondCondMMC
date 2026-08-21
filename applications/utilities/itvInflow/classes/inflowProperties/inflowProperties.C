/*
 * inflowProperties.C
 *
 *  Created on: Aug 15, 2011
 *      Author: gregor
 */

#include "inflowProperties.H"

namespace Foam
{
    inflowProperties::inflowProperties(const Time& runTime, const IOdictionary inflowPropDict)
    :
    intL_(inflowPropDict.lookupOrDefault<scalar>("integralLength", 0.0)),
    var_(inflowPropDict.lookup("var")),
    diffTime_(runTime.endTime().value()-runTime.startTime().value()),
    ubulk_(inflowPropDict.lookupOrDefault<scalar>("Ubulk", 0.0)),
    patchName_(inflowPropDict.lookup("patchName")),
    divFree_(inflowPropDict.lookup("divFree")),
    lambda_(inflowPropDict.lookupOrDefault<scalar>("lambda", 0.0))
    {
        // Kempf formula to estimate a time step that diffusion takes 20 time steps
        diffDeltaT_ = sqr(intL_)/(Foam::constant::mathematical::pi*2*20*1);
        if(lambda_ == 0.0){
            Info << "Error: Please set lambda correctly" << endl;
        }
    }
}
