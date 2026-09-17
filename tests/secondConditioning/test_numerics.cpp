#include "SecondConditioningNumerics.H"
#include <cassert>
#include <iostream>
#include <random>
using namespace Foam::secondConditioningNumerics;

double rhs(double y) { return 9060*(1-y)*std::exp(20*(y-1)); }
double referenceReaction(double y, double dt)
{
    const int n=30000; const double h=dt/n;
    for(int i=0;i<n;++i) {
        const double a=rhs(y),b=rhs(y+h*a/2),c=rhs(y+h*b/2),d=rhs(y+h*c);
        y += h*(a+2*b+2*c+d)/6;
    }
    return y;
}
int main()
{
    for(double y : {0., 0.2, 0.9, 0.999, 1.})
        for(double dt : {0., 1e-9, 1e-6, 1e-3, 0.01}) {
            const double next=progressReaction(y,dt,9060,20,.05);
            assert(next>=y && next<=1);
            assert(progressReaction(y,dt,0,20,.05)==y);
        }
    assert(std::abs(progressReaction(.3,.02,200,0,.05)-(1-.7*std::exp(-4)))<1e-13);
    const double ref=referenceReaction(.86,3e-4);
    const double e1=std::abs(progressReaction(.86,3e-4,9060,20,.1)-ref);
    const double e2=std::abs(progressReaction(.86,3e-4,9060,20,.05)-ref);
    const double e3=std::abs(progressReaction(.86,3e-4,9060,20,.025)-ref);
    assert(e1/e2>3 && e2/e3>3);
    bool caught=false;
    try { progressReaction(1.1,.001,1,20,.05); } catch(const std::domain_error&) { caught=true; }
    assert(caught);
    std::cout << "Reaction: bounded, zero-source, analytic Z=0 and second-order convergence passed\n";

    const std::vector<std::pair<int,int>> edges{{0,1},{0,3},{1,2},{1,4},{2,3},{3,4}};
    std::vector<bool> visited(edges.size(),false);
    for(unsigned phase=0;phase<edges.size();++phase) {
        const auto partner=neighbourPartners(6,edges,phase);
        assert(partner[5]==-1);
        for(int i=0;i<6;++i) if(partner[i]>=0) assert(partner[partner[i]]==i);
        for(unsigned e=0;e<edges.size();++e)
            visited[e]=visited[e] || partner[edges[e].first]==edges[e].second;
    }
    for(bool v:visited) assert(v);
    assert(neighbourPartners(1,{},0)[0]==-1);
    std::cout << "Processor groups: symmetric/disjoint matching, isolated rank and eventual edge coverage passed\n";

    std::mt19937_64 generator(451);
    std::normal_distribution<double> normal;
    std::uniform_real_distribution<double> uniform;
    const int n=250000;
    double mean=0,var=0,cov=0;
    for(int i=0;i<n;++i) {
        const double old=normal(generator), next=ouUpdate(old,.2,1,normal(generator));
        mean+=next; var+=next*next; cov+=old*next;
    }
    assert(std::abs(mean/n)<.012 && std::abs(var/n-1)<.015);
    assert(std::abs(cov/n-std::exp(-.2))<.015);
    assert(ouUpdate(2,0,1,3)==2);
    assert(ouUpdate(0,1e-18,1,1)>0); // cancellation regression
    assert(modifiedProgress(.4,0,12)==.4);
    assert(modifiedProgress(1,.05,2)>1);
    assert(modifiedProgress(1,0,0,.1)==1.1);
    std::cout << "OU: stationary mean/variance, correlation, tiny dt and modified-reference limits passed\n";

    for(double dt : {.001,.1,1.,10.}) {
        double second=0,fourth=0;
        for(int i=0;i<n;++i) {
            const double a=mixingExtent(dt,1,true,uniform(generator),uniform(generator));
            assert(a>=0 && a<=1);
            const double v=(1-a)*(1-a);
            second+=v; fourth+=v*v;
            const double p=.1+a*(.66-.1), q=.9+a*(.66-.9);
            assert(std::abs(.3*p+.7*q-.66)<1e-14);
        }
        second/=n; fourth/=n;
        const double se=std::sqrt(std::max(0.,fourth-second*second)/n);
        assert(std::abs(second-std::exp(-2*dt))<6*se+1e-5);
        const double a=mixingExtent(dt,1,false,0,0);
        assert(std::abs((1-a)*(1-a)-std::exp(-2*dt))<1e-14);
    }
    assert(mixingExtent(.2,1,true,.125,.5)==mixingExtent(.2,1,true,.125,.5));
    std::cout << "Mixing: random-event variance decay, deterministic limit, bounds, conservation and shared-pair draws passed\n";
}
