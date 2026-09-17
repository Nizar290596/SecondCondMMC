#!/usr/bin/env python3
"""Compile the actual pairing and triplet-exchange methods using minimal containers.
This verifies the algorithms, not OpenFOAM integration or MPI communication.
"""
from pathlib import Path
import subprocess
import tempfile
ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'src/lagrangian/mmc/submodels/Mixing/MixingModel/ParticleInteractionModels/sparseParticleModel'

def block(text, marker):
    start=text.index(marker); end=text.index('{', start)+1; depth=1
    while depth:
        depth += (text[end]=='{')-(text[end]=='}'); end+=1
    return text[start:end]
source=(BASE/'mixParticleModel.C').read_text()
methods='\n'.join('template<class CloudType>\n'+block(source,'void Foam::mixParticleModel<CloudType>::'+name) for name in ['findPairs','KkdTreeLikeSearch','SmixList'])
compare=block((BASE/'mixParticleModel.H').read_text(),'class lessArg').replace('FatalError << "Comparison operator did not work properly"<<exit(FatalError);','throw std::runtime_error("Invalid comparator");')
prefix=r'''
#include <vector>
#include <array>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>
#include <stdexcept>
#define forAll(a,i) for(int i=0;i<int((a).size());++i)
namespace Foam {
using label=int; using scalar=double; using std::max; using std::min;
constexpr double GREAT=1e30;
using std::ceil;
inline double mag(double x){return std::abs(x);}
inline double sqr(double x){return x*x;}
struct vector: std::array<double,3> {
    vector(double x=0,double y=0,double z=0):std::array<double,3>{{x,y,z}}{}
    vector operator-(const vector& b) const {return {(*this)[0]-b[0],(*this)[1]-b[1],(*this)[2]-b[2]};}
};
inline double magSqr(const vector& v){return v[0]*v[0]+v[1]*v[1]+v[2]*v[2];}
template<class T> struct List:std::vector<T> {
    using std::vector<T>::vector;
    using std::vector<T>::operator=;
    void setSize(size_t n){this->resize(n);}
    List& operator=(const T& x){std::fill(this->begin(),this->end(),x);return *this;}
    void append(const T& x){this->push_back(x);}
    void append(T&& x){this->push_back(std::move(x));}
};
template<class T> using DynamicList=List<T>;
struct Clock {
    double value()const{return .2;}
    const Clock& mesh()const{return *this;} const Clock& time()const{return *this;}
    const Clock& deltaT()const{return *this;}
};
template<class CloudType> struct mixParticleModel {
    struct eulerianFieldData {
        vector p; List<double> r=List<double>(4,0.0); int index=0,rank=0;
        const auto& position()const{return p;} const auto& XiR()const{return r;}
        int particleIndex()const{return index;}
        bool local()const{return rank==0;}
    };
    struct Particle {double value=0,weight=1,exposure=0;};
    List<double> Xii_=List<double>(4,1.0);
    bool physicalLocalization_=true;
    double ri_=1,retainedPairFraction_=1,maxPairDistance_=0;
    mutable List<int> splitAxisHistogram_;
    DynamicList<eulerianFieldData> eulerianFieldDataList_;
    DynamicList<List<int>> particlePairs_;
    List<Particle> particleList_;
    Clock clock;
    const Clock& owner()const{return clock;}
    void mixpair(Particle& p,const eulerianFieldData&,Particle& q,const eulerianFieldData&,double& dt){
        const double mean=(p.weight*p.value+q.weight*q.value)/(p.weight+q.weight);
        const double a=1-std::exp(-dt);
        p.value+=a*(mean-p.value);q.value+=a*(mean-q.value);
        p.exposure+=dt;q.exposure+=dt;
    }
    void findPairs(const DynamicList<eulerianFieldData>&,DynamicList<List<label>>&)const;
    void KkdTreeLikeSearch(const DynamicList<eulerianFieldData>&,label,label,std::vector<label>&,std::vector<label>&,std::vector<label>&)const;
    void SmixList();
'''
suffix=r'''
int main(){
    using M=Foam::mixParticleModel<int>;
    for(int dimensions : {1,4}) for(int n=0;n<=150;++n){
        M m;Foam::DynamicList<M::eulerianFieldData> data(n);Foam::DynamicList<Foam::List<int>> pairs;
        m.Xii_.resize(dimensions);
        for(int i=0;i<n;++i){data[i].r.resize(dimensions);data[i].r[0]=i%7;data[i].p[0]=i%11;}
        m.findPairs(data,pairs);std::vector<int> seen(n);
        for(const auto& g:pairs){assert(g.size()==2||g.size()==3);for(int i:g){assert(i>=0&&i<n);++seen[i];}}
        for(int count:seen)assert(count==(n<2?0:1));
    }
    M m;Foam::DynamicList<M::eulerianFieldData> data(4);Foam::DynamicList<Foam::List<int>> pairs;
    const double x[]={0,100,.01,100.01};
    for(int i=0;i<4;++i)data[i].p[0]=x[i];
    m.findPairs(data,pairs);
    for(const auto& g:pairs)assert(std::abs(x[g[0]]-x[g[1]])<.02);
    m.maxPairDistance_=.001;m.findPairs(data,pairs);assert(pairs.empty());
    m.maxPairDistance_=0;m.retainedPairFraction_=.5;m.findPairs(data,pairs);assert(pairs.size()==1);
    m.retainedPairFraction_=1;m.Xii_[0]=1e-6;
    for(int i=0;i<4;++i)data[i].r[0]=(i<2?0:1);
    m.findPairs(data,pairs);
    for(const auto& g:pairs)assert(data[g[0]].r[0]==data[g[1]].r[0]);
    // Exchange the same triplet on each rank, including its remote-remote edge.
    std::array<double,3> allValues{};
    for(int localRank=0;localRank<3;++localRank){
        M t;t.particleList_.resize(3);t.eulerianFieldDataList_.resize(3);
        t.particlePairs_.append(Foam::List<int>{0,1,2});
        double before=0;
        for(int i=0;i<3;++i){
            t.eulerianFieldDataList_[i].index=i;
            t.eulerianFieldDataList_[i].rank=(i==localRank?0:1);
            t.particleList_[i].value=i*.4;t.particleList_[i].weight=i+1;
            before+=t.particleList_[i].value*t.particleList_[i].weight;
        }
        t.SmixList();double after=0;
        for(int i=0;i<3;++i){
            assert(std::abs(t.particleList_[i].exposure-.2)<1e-14);
            after+=t.particleList_[i].value*t.particleList_[i].weight;
            if(localRank==0)allValues[i]=t.particleList_[i].value;
            else assert(t.particleList_[i].value==allValues[i]);
        }
        assert(std::abs(before-after)<1e-14);
    }
    std::cout<<"Pairing: N=0..150, no duplicates, physical/progress dominance, distance cap, rejection, triplet conservation/exposure and rank-copy consistency passed\n";
}
'''
with tempfile.TemporaryDirectory() as tmp:
    cpp=Path(tmp)/'pairing.cpp';exe=Path(tmp)/'pairing'
    cpp.write_text(prefix+compare+';\n};\n}\n'+methods+suffix)
    subprocess.run(['g++','-std=c++17','-O1','-g','-D_GLIBCXX_ASSERTIONS','-fsanitize=undefined',str(cpp),'-o',str(exe)],check=True)
    subprocess.run([str(exe)],check=True)
