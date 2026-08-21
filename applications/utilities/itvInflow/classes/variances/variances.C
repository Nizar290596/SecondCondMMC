/*
 * variances.C
 *
 *  Created on: Aug 15, 2011
 *      Author: gregor
 */

#include "variances.H"

namespace Foam
{

variances::variances(volVectorField &U_)
:
U(U_),
Mesh(U_.mesh()),
variance(vector::zero)
{
    forAll(Mesh.V(),Iter){
        MeshV += Mesh.V()[Iter];
    }
    update();
}

void variances::update(){
    // first zero mean check

    vector meanOffs(vector::zero);
    forAll(U,Iter){
        meanOffs += U[Iter];
    }

    meanOffs /=  Mesh.V().size();

    forAll(U,Iter){
        U[Iter] -= meanOffs;         // mean = 0
    }

    scalar alpha(0),beta(0),n(0);
    variance = vector::zero;

    forAll(U,iter){
            n=iter+1;
            alpha = (n-1)/n;
            beta = 1/n;
            for(int i=0;i<3;i++){
                variance[i] = alpha*variance[i] + beta*sqr(U[iter][i]);                     // variance first step. mean is assumed to be zero
            }
        }

        // Info << "mean: " << mean << endl;
        // Info << "variance: " << variance << endl;
}

void variances::normalise(){
    update();

    forAll(U, celli)                                        // normalize with to obtain a avg(u*u) = 1 field
    {
            for(int i=0;i<3;i++){
                U[celli][i] =  U[celli][i]/Foam::sqrt(variance[i]);
            }
    }
}

void variances::normalise2(){                                // keep a unity mean variance
    update();
    scalar meanVarSqrt = Foam::sqrt((variance[0]+variance[1]+variance[2])/3.0);
    scalar invMeanVarSqrt = 1/meanVarSqrt;

    Info << "meanVarSqrt" << meanVarSqrt << endl;
    forAll(U, celli)                                        // normalize with to obtain a avg(u*u) = 1 field
    {
            for(int i=0;i<3;i++){
                U[celli][i] =  U[celli][i]*invMeanVarSqrt;
            }
    }
}

void variances::scale(tensor rms){
    update();
    Info << "scaling variances" << endl;
    scalar a11,a21,a22,a31,a32,a33;
    a11 = Foam::sqrt(rms.xx());
    a21 = rms.xy()/a11;
    a22 = Foam::sqrt(rms.yy()-(a21*a21));
    a31 = rms.xz()/a11;//rms.xz()*a11;
    a32 = (rms.yz()-a21*a31)/a22;
    a33 = Foam::sqrt(rms.zz()-(a31*a31)-(a32*a32));
    vector Uunscaled(vector::zero);

    forAll(U, celli)                                        // normalize with rms to obtain a avg(u*u) = 1 field
        {
            Uunscaled     = U[celli];
            U[celli][0] = a11 * Uunscaled[0];
            U[celli][1] = a21 * Uunscaled[0] + a22 * Uunscaled[1];
            U[celli][2] = a31 * Uunscaled[0] + a32 * Uunscaled[1] + a33 * Uunscaled[2];
        }
}

void variances::corTensor(){                // calculate the corellation tensor
    tensor corr(tensor::zero);
    scalar alpha(0),beta(0),n(0);

    forAll(U, celli){
        n=celli+1;
        alpha = (n-1)/n;
        beta = 1/n;

        corr.xx() = alpha  * corr.xx() + beta * (U[celli][0]*U[celli][0]);
        corr.yy() = alpha  * corr.yy() + beta * (U[celli][1]*U[celli][1]);
        corr.zz() = alpha  * corr.zz() + beta * (U[celli][2]*U[celli][2]);
        corr.xy() = alpha  * corr.xy() + beta * (U[celli][0]*U[celli][1]);
        corr.xz() = alpha  * corr.xz() + beta * (U[celli][0]*U[celli][2]);
        corr.yz() = alpha  * corr.yz() + beta * (U[celli][1]*U[celli][2]);
    }

    Info << "avg R tensor" << corr << endl;
}

} // End namespace Foam


