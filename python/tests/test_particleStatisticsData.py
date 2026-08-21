from mmcFoamReader import particleStatisticsData
import numpy as np

def test_particleStatisticsData():
    pData = particleStatisticsData('tests/particleStatisticsData','pData')
    x = pData.getfield('x')
    v = pData.getfield('value')
    assert (v == [1.0,2.0,3.0,4.0]).all()
    
    binPos,binVal, binStdDev = pData.conditionalStatistics("x","value",bins=3)
    
    assert (binPos == [1.5,2.5,3.5]).all()
    assert (binVal == [1.0,2.0,3.5]).all()

def test_particleStatisticsData_weighted():
    pData = particleStatisticsData('tests/particleStatisticsData','pData')
    v = pData.getfield('value')
    assert (v == [1.0,2.0,3.0,4.0]).all()
    # Default weight
    wt = np.ones(len(v))*0.5
    binPos,binVal, binStdDev = pData.conditionalStatistics("x","value",weight=wt,bins=3)
    
    assert (binPos == [1.5,2.5,3.5]).all()
    assert (binVal == [1.0,2.0,3.5]).all()

def test_particleStatisticsData_provideData():
    pData = particleStatisticsData('tests/particleStatisticsData','pData')
    v = pData.getfield('value')
    x = pData.getfield('x')
    assert (v == [1.0,2.0,3.0,4.0]).all()
    # Default weight
    wt = np.ones(len(v))*0.5
    binPos,binVal, binStdDev = pData.conditionalStatistics(x,v,weight=wt,bins=3)
    
    assert (binPos == [1.5,2.5,3.5]).all()
    assert (binVal == [1.0,2.0,3.5]).all()


def test_particleStatisticsData_reynoldsAverage():
    pData = particleStatisticsData('tests/particleStatisticsData','reynoldsAverageData')
    binPos,binVal, binStdDev = pData.reynoldsAveragedConditionalStatistics("x","value",bins=4)
    print(binPos)
    assert(binPos == [0.5,1.5,2.5,3.5]).all()
    assert(binVal == [2.5,2.5,2.5,2.5]).all()


