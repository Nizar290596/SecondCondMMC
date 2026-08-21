from mmcFoamReader import particleStatisticsData
import numpy as np

def test_particleStatistics():
    pData = particleStatisticsData('tests/particleStatisticsData','pData')
    x = pData.getfield('x')
    v = pData.getfield('value')
    assert (v == [1.0,2.0,3.0,4.0]).all()
    # Default weight
    weight = np.ones(len(x))
    binPos,binVal, binStdDev = pData.conditionalStatistics(v,x,weight,3)
    
    assert (binPos == [1.5,2.5,3.5]).all()
    assert (binVal == [1.0,2.0,3.5]).all()

