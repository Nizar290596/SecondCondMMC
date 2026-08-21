/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

    This file is not part of OpenFOAM but an original routine developed 
    in the OpenFOAM techonology.

    You can redistribute it and/or modify it under the terms of the GNU General
    Public License as published by the Free Software Foundation; either version
    3 of the License, or (at your option) any later version. See GNU General
    Public License at <http://www.gnu.org/licenses/gpl.html>

    OpenFOAM is a trademark of OpenCFD.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.


Class
    Foam::syntheticTurbulence


Authors
    Andrea Montorfano, Federico Piscaglia
    Dipartimento di Energia, Politecnico di Milano
    via Lambruschini 4
    I-20156 Milano (MI)
    ITALY

Contact
    andrea.montorfano@polimi.it , ph. +39 02 2399 3909
    federico.piscaglia@polimi.it, ph. +39 02 2399 8620    


        
\*---------------------------------------------------------------------------*/

#include "syntheticTurbulence.H"

void Foam::syntheticTurbulence::update()
{

    if
    (
        curTimeIndex_ != this->runTime_.timeIndex() 
    )
    {
        //Pout << "Generation of the fluctuations" <<endl;
            
        scalarField phiAng(nModes_,0.0);
        scalarField alpha(nModes_,0.0);
        scalarField psi(nModes_,0.0);
        scalarField theta(nModes_,0.0);

        vectorField ut(coords_.size(),pTraits<vector>::zero);

        Random randNum(runTime_.timeIndex());                
    
        // generate angles

        if(Pstream::master())
        {
            for (label n=0;n<nModes_;n++)
            {
                
                phiAng[n] = 2.0*constant::mathematical::pi*randNum.sample01<scalar>();
                alpha[n] = 2.0*constant::mathematical::pi*randNum.sample01<scalar>();
                psi[n]   = 2.0*constant::mathematical::pi*randNum.sample01<scalar>();

                theta[n] = 0.0;

                bool flag(false);

                while (!flag)
                {

                    scalar x = constant::mathematical::pi*randNum.sample01<scalar>();
                    scalar y = 0.5*randNum.sample01<scalar>();

                    scalar fx = 0.5*sin(x);

                    if ( y <= fx )
                    {
                        flag = true;
                        theta[n] = x;
                    }
                }

            }

            if(Pstream::parRun())
            {
                for (label slave : Pstream::subProcs())
                {

                    if (debug)
                    {
                        Pout << "sending angles to slave: " << slave << endl;
                    }

                    OPstream toSlave(Pstream::commsTypes::blocking, slave);
                    

                    toSlave << phiAng
                        << alpha
                        << theta
                        << psi; 

                    
                }
            }


        }

        else
        {
            IPstream fromMaster(Pstream::commsTypes::blocking, Pstream::masterNo());

            if (debug)
            {
                Pout << "receiving angles on slave: " 
                    << Pstream::myProcNo() << endl;
            }
         
            fromMaster >> phiAng;
            fromMaster >> alpha;
            fromMaster >> theta;
            fromMaster >> psi;            
        }


        if (debug)
        {
            Pout << "Number of values calculated (out of "
                << coords_.size() << ")" << endl;
        }

        scalar timeStart=runTime_.elapsedClockTime();
        scalar newTime=0;

        scalarField deltak(nModes_,0.0);

        forAll(coords_,i)
        {

    
            // --------------------------------------------------------------
            // generate wavenumbers

            if ((i%1000==0)&&(i!=0) && debug) 
            {
                newTime=runTime_.elapsedClockTime()-timeStart;
                scalar percent = static_cast<scalar>
                (
                    static_cast<label>
                    (
                        static_cast<scalar>(i)/static_cast<scalar>
                        (
                            k_.size()
                        )*10000.0
                    )
                )/100.0;
                scalar timeRem = 1.0/(percent+SMALL)*100*newTime - newTime;
                scalar remH = static_cast<label>(timeRem/3600.0);
                scalar remMin = static_cast<label>((timeRem-remH*3600.0)/60.0);
                scalar remS = static_cast<label>
                (
                    timeRem-remH*3600.0-remMin*60.0
                );
                Pout << i <<"\t" << percent 
                    << "%,\tremaining time (estimated): " 
                    <<  remH << "h " 
                    << remMin << "min " 
                    << remS << "s" << endl;
            }



               if (k_[i] > minK_)
            {

                kw_[0] = kwZero_[i];


            if (distType_== "uniform" || distType_==string::null) 
                {
                    uniDist_(kw_);
                }
                else if (distType_ == "logarithmic")
                {
                    logDist_(kw_);                
                }
                else
                {
                    WarningIn("Foam::syntheticTurbulence::update()") << endl
                        << "distribution type not specified. "
                        << "Assuming uniform distribution" << endl;

                    uniDist_(kw_);
                }


                 deltak[0] = kw_[1] - kw_[0];                
                   deltak[nModes_-1] = kw_[nModes_-1] - kw_[nModes_-2];
                
                for(label n=1;n<nModes_-1;n++)
                {
                       scalar dkp;
                      scalar dkn;

                    dkp = (kw_[n] - kw_[n-1]) / 2.0;
                    dkn = (kw_[n+1] - kw_[n]) / 2.0;
                    deltak[n] = dkp + dkn;
                }

                   // kn: wavenumber vector, kappa^n_j, 

                scalarField k0 = kw_*sin(theta)*cos(phiAng);
                scalarField k1 = kw_*sin(theta)*sin(phiAng);
                scalarField k2 = kw_*cos(theta);

                vectorField kn(nModes_,pTraits<vector>::zero);

                kn.replace(0,k0);
                kn.replace(1,k1);
                kn.replace(2,k2);

                scalarField s0 = cos(alpha)*cos(theta)*cos(phiAng) -
                    sin(alpha)*sin(phiAng);

                scalarField s1 = cos(alpha)*cos(theta)*sin(phiAng) 
                    + sin(alpha)*cos(phiAng);

                scalarField s2 = -cos(alpha)*sin(theta);

                vectorField sigma(nModes_,pTraits<vector>::zero);

                sigma.replace(0,s0);
                sigma.replace(1,s1);
                sigma.replace(2,s2);

                // -------------------------------------------
                // calculate spectrum Ek(k)

                scalarField Ek = fSpect_*MVKS_(i, ke_[i], kw_, kkol_[i]); 

                if ( min(Ek*deltak) < 0.0)
                {
                    FatalErrorIn("Foam::syntheticTurbulence::update()")
                        << "Negative energy spectrum: " << endl
                        << exit(FatalError);
                }

                scalarField uHat = 1*sqrt(deltak*Ek);

                ut[i] += 2.0*sum(uHat*cos((kn & coords_[i]) + psi)*sigma);

            }


            
        }    // end of loop on faces

        // --------------------------------------------------------------
        // Billson filter                                 

        //original milano
        filter(ut);
               curTimeIndex_ = runTime_.timeIndex();

        //variante freiberg
        //curTimeIndex_ = runTime_.timeIndex();
        //filter(ut);

    } //endif curTimeIndex
    
} // end of function syntheticTurbulence::update()
