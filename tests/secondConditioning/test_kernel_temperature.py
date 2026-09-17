#!/usr/bin/env python3
"""Compile production reconstruction with aliased target/carrier temperature."""
from pathlib import Path
import os
import re
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
source = (root/'src/lagrangian/mmc/submodels/Thermodynamic/ThermoPhysicalCouplingModel/KernelEstimation/KernelEstimation.C').read_text()
start = source.index('void Foam::KernelEstimation<CloudType>::reconstructTemperatureTargets')
end = source.index('{', start) + 1
depth = 1
while depth:
    depth += (source[end] == '{') - (source[end] == '}')
    end += 1
method = 'template<class CloudType>\n' + source[start:end]
initialization = re.search(r'    const scalarField initialTemperature[^\n]+\n    scalarField thermalTarget[^\n]+', source)[0]
accumulation = re.search(r'            thermalTarget\[celli\]\s*= sumWtT/sumWt;', source)[0]
inert_sum = re.search(r'            if \(YEqvETarget\[specieI\].name\(\) != inertSpecie_\)\n\s*Yt \+= YEqvETarget\[specieI\];', source)[0]
# Prevent reintroducing enthalpy writes to T before reconstruction.
compute = source[source.index('void Foam::KernelEstimation<CloudType>::computeTargets'):source.index('// * * * * * * * * * * * * * * * * Constructors')]
assert not re.search(r'\bTEqvETarget\[celli\]\s*=', compute)

prefix = r'''
#include <vector>
#include <cmath>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>
#define forAll(v,i) for (std::size_t i=0; i<(v).size(); ++i)
#define FatalErrorInFunction Foam::errors
namespace Foam {
using scalar=double;
using scalarField=std::vector<double>;
struct volScalarField:scalarField {
    using scalarField::vector;
    std::string name_;
    const std::string& name()const{return name_;}
    const scalarField& primitiveField() const {return *this;}
};
scalarField& operator+=(scalarField& a,const volScalarField& b){
    forAll(a,i) a[i]+=b[i];
    return a;
}
template<class T>using PtrList=std::vector<T>;
constexpr int FatalError=1;
inline int exit(int){throw std::runtime_error("fatal inversion input");}
struct Error {template<class T>Error& operator<<(const T&){return *this;}} errors;
struct Mixture {
    const volScalarField* carrier;
    double THa(double h,double p,double t0) const {
        if (!(t0>0)) throw std::runtime_error("Negative initial temperature T0");
        assert(p>0);
        // A thermo callback must never observe temporary enthalpy in T.
        for (double t:*carrier) assert(t>0);
        return (h+600000)/1000;
    }
};
struct Composition {
    const volScalarField* carrier;
    Mixture particleMixture(const scalarField& y) const {
        assert(y.size()==2 && std::abs(y[0]+y[1]-1)<1e-12);
        return {carrier};
    }
};
struct Owner {
    volScalarField temp{300,1200,500};
    volScalarField pressure{101325,101325,101325};
    const volScalarField& T()const{return temp;}
    const volScalarField& p()const{return pressure;}
    Composition composition()const{return {&temp};}
};
template<class CloudType>struct KernelEstimation {
    Owner data;
    bool coupleEnthalpy_=true;
    scalarField indicator{1,1,0};
    Owner& owner(){return data;}
    const scalarField& Indicator()const{return indicator;}
    void reconstructTemperatureTargets(const scalarField&,const scalarField&,
        const PtrList<volScalarField>&,volScalarField&);
    void exercise(bool alias) {
'''
middle = r'''
        volScalarField separate=this->owner().T();
        volScalarField& target=alias ? data.temp : separate;
        PtrList<volScalarField> species{{0.2,0.3,0.1},{0.8,0.7,0.9}};
        const double sumWt=2;
        const double moments[]={-273141.667194,700000,-999999};
        for (std::size_t celli=0;celli<2;++celli) {
            const double sumWtT=sumWt*moments[celli];
'''
suffix = r'''
        }
        reconstructTemperatureTargets(thermalTarget,initialTemperature,species,target);
        assert(std::abs(target[0]-326.858332806)<1e-8);
        assert(std::abs(target[1]-1300)<1e-12);
        assert(target[2]==500); // Unsupported cell retains its temperature.
        if (!alias) assert(data.temp[0]==300 && data.temp[1]==1200);
        coupleEnthalpy_=false;
        reconstructTemperatureTargets({400,1500,-1},initialTemperature,species,target);
        assert(target[0]==400 && target[1]==1500 && target[2]==500);
        coupleEnthalpy_=true;
        bool rejected=false;
        try {reconstructTemperatureTargets(thermalTarget,{-1,1200,500},species,target);}
        catch(const std::runtime_error&){rejected=true;}
        assert(rejected);
    }
};
}
'''
main = r'''
void checkInertClosure(){
    Foam::PtrList<Foam::volScalarField> YEqvETarget{{0.2,0.3},{0.8,0.7}};
    YEqvETarget[0].name_="fuel"; YEqvETarget[1].name_="N2";
    const std::string inertSpecie_="N2";
    Foam::scalarField Yt(2,0);
    forAll(YEqvETarget,specieI){
''' + inert_sum + r'''
    }
    assert(std::abs(1-Yt[0]-0.8)<1e-12);
    assert(std::abs(1-Yt[1]-0.7)<1e-12);
}
int main(){
    checkInertClosure();
    // Reproduce the former enthalpy-as-T0 error, without an OpenFOAM build.
    Foam::Owner owner;
    owner.temp[0]=-273141.667194;
    bool reproduced=false;
    try {Foam::Mixture{&owner.temp}.THa(owner.temp[0],101325,owner.T()[0]);}
    catch(const std::runtime_error&){reproduced=true;}
    assert(reproduced);
    Foam::KernelEstimation<int>{}.exercise(true);
    Foam::KernelEstimation<int>{}.exercise(false);
    std::cout << "Kernel thermal targets: negative absolute h, shared/separate T, unsupported cells and temperature mode passed\n";
}
'''
with tempfile.TemporaryDirectory() as temp:
    src=Path(temp)/'kernel_temperature.cpp'
    exe=Path(temp)/'kernel_temperature'
    src.write_text(prefix+initialization+middle+accumulation+suffix+method+main)
    subprocess.run([os.environ.get('CXX','g++'),'-std=c++17','-O2','-Wall','-Wextra','-Werror',
                    '-fsanitize=undefined',str(src),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)
