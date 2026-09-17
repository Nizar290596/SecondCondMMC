#!/usr/bin/env python3
"""Exercise actual cloud reference updates with lightweight particle fixtures."""
from pathlib import Path
import subprocess
import tempfile
ROOT=Path(__file__).resolve().parents[2]
source=(ROOT/'src/lagrangian/mmc/clouds/Templates/MixingPopeCloud/MixingPopeCloud.C').read_text()
def block(marker):
    start=source.index(marker); end=source.index('{',start)+1; depth=1
    while depth:
        depth+=(source[end]=='{')-(source[end]=='}');end+=1
    return source[start:end]
methods='\n'.join('template<class CloudType>\n'+block('void Foam::MixingPopeCloud<CloudType>::'+name) for name in ['updatePhiReaction','updateOUProcess','initializeSecondConditioningState'])
prefix=r'''
#include "SecondConditioningNumerics.H"
#include <vector>
#include <random>
#include <cassert>
#include <iostream>
#define forAllIters(c,i) for(auto i=(c).begin();i!=(c).end();++i)
#define FatalErrorInFunction Foam::errorStream
namespace Foam {
using scalar=double; using std::min; using std::max;
const int FatalError=1;
inline int exit(int){throw std::runtime_error("OpenFOAM fatal error");}
struct ErrorStream{template<class T>ErrorStream& operator<<(const T&){return *this;}}errorStream;
inline double OUStateUpdate(double w,double dt,double tau,double g){return secondConditioningNumerics::ouUpdate(w,dt,tau,g);}
struct Particle {
    double phi_=0,modified_=0,omega_=0,age_=0,temperature_=298;int flag_=0;
    double& phi(){return phi_;} double& phiModified(){return modified_;}
    double& omegaOU(){return omega_;} double& burnedAge(){return age_;}
    double& T(){return temperature_;} int& secondCondFlag(){return flag_;}
};
struct Random {
    std::mt19937 gen{1234};
    double Normal(int,int){return std::normal_distribution<double>{}(gen);}
    double RandomValue(){return std::uniform_real_distribution<double>{}(gen);}
};
struct Generator:Random{double Random(){return RandomValue();}};
template<class CloudType>struct MixingPopeCloud {
    using particleType=Particle;
    std::vector<Particle> particles;
    struct Iter {
        std::vector<Particle>::iterator i;
        bool operator!=(const Iter& b)const{return i!=b.i;}
        Iter& operator++(){++i;return *this;}
        Particle& operator()(){return *i;}
    };
    Iter begin(){return {particles.begin()};}Iter end(){return {particles.end()};}
    bool enabled=true;
    bool secondCondMixingEnabled()const{return enabled;}
    double secondCondAPhi_=9060,secondCondZPhi_=20,secondCondMaxSourceStep_=.05;
    double secondCondBeta_=.05,secondCondTauOU_=33e-6;
    double secondCondAgeRate_=0,secondCondAgeThreshold_=.99;
    double secondCondR_=.25,secondCondTu_=298,secondCondTb_=2248;
    Generator secondCondRandom_;
    void updatePhiReaction(scalar);
    void updateOUProcess(scalar);
    void initializeSecondConditioningState(particleType&);
};
}
'''
suffix=r'''
int main(){
    Foam::MixingPopeCloud<int> cloud; cloud.particles.resize(10000);
    int flagged=0;
    for(auto& p:cloud.particles){p.T()=1273;cloud.initializeSecondConditioningState(p);flagged+=p.secondCondFlag();assert(p.phi()==.5);}
    assert(flagged>2300&&flagged<2700);
    // Regression: disabled OU must still refresh the modified reference.
    cloud.secondCondBeta_=0; cloud.secondCondTauOU_=0;
    cloud.updatePhiReaction(1e-5);cloud.updateOUProcess(1e-5);
    for(auto& p:cloud.particles){assert(p.phi()>.5);assert(p.phiModified()==p.phi());assert(p.omegaOU()==0);}
    // A=0 must not prevent refreshing a phi changed by dense mixing.
    cloud.secondCondAPhi_=0;
    for(auto& p:cloud.particles)p.phi()=.8;
    cloud.updatePhiReaction(1e-5);cloud.updateOUProcess(1e-5);
    for(auto& p:cloud.particles)assert(p.phiModified()==.8);
    cloud.secondCondAgeRate_=2;
    for(auto& p:cloud.particles)p.phi()=1;
    cloud.updateOUProcess(.1);
    for(auto& p:cloud.particles){assert(std::abs(p.burnedAge()-.1)<1e-14);assert(std::abs(p.phiModified()-1.2)<1e-14);}
    cloud.enabled=false;Foam::Particle legacy;legacy.flag_=1;legacy.age_=2;
    cloud.initializeSecondConditioningState(legacy);
    assert(legacy.secondCondFlag()==0&&legacy.phi()==0&&legacy.burnedAge()==0);
    std::cout<<"Cloud reference: subset initialization, reaction on every particle, beta/tau zero refresh, A=0 refresh, optional age and disabled initialization passed\n";
}
'''
with tempfile.TemporaryDirectory() as tmp:
    cpp=Path(tmp)/'reference.cpp';exe=Path(tmp)/'reference'
    cpp.write_text(prefix+methods+suffix)
    subprocess.run(['g++','-std=c++17','-O1','-g','-D_GLIBCXX_ASSERTIONS','-fsanitize=undefined','-I'+str(ROOT/'src/lagrangian/mmc/submodels/Mixing'),str(cpp),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)
