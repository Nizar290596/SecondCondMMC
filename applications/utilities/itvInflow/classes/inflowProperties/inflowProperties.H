/*
 * inflowProperties.H
 *
 *  Created on: Aug 15, 2011
 *      Author: gregor
 *
 *    A simple class to hold inflow properties to
 *    assure data consistency
 */

#ifndef    INFPROP_H_
#define INFPROP_H_

#include "volFields.H"
#include "OFstream.H"
#include "fvCFD.H"

namespace Foam
{

class inflowProperties {

private:

    const scalar intL_;    // desired intergral length
    const tensor var_;    // desired variance tensor
    const scalar diffTime_;    //  
    const scalar ubulk_;    // 
    const word patchName_;    //
    scalar diffDeltaT_;
    bool divFree_;
    const scalar lambda_;
    
public:
    inflowProperties(const Time& runTime,const IOdictionary inflowPropDict);    //
    
    scalar intL(){return intL_;}
    
    tensor var(){return var_;}

    scalar diffTime(){return diffTime_;}

    scalar ubulk(){return ubulk_;}

    scalar deltaT(){return diffDeltaT_;}

    word patchName(){return patchName_;}

    bool divFree(){return divFree_;}

    scalar lambda(){return lambda_;}
};
}
#endif /* RADIALDATA_H_ */
