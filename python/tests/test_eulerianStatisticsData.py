from mmcFoamReader import eulerianStatisticsData
import numpy as np

"""Check that the data is read in correctly"""
def test_readEulerianStats():
    x,y,v = createFile()
    data = eulerianStatisticsData("tests/dataFile.dat")

    assert (data.values() == np.array(v)).all()

    # Create a position array from x and y
    pos = np.array([x,y])
    pos = pos.transpose()

    assert (data.pos() == pos).all()


def func(x,y):
    return x+y

def createFile():
    rng = np.random.default_rng()

    x = rng.random(10) - 0.5
    y = rng.random(10) - 0.5
    v = func(x, y)
    f = open("tests/dataFile.dat","w")
    f.write("(\n")
    for e in v:
        f.write(str(e)+"\n")
    f.write(")\n")
    f.write("(\n")
    for ind in range(len(x)):
        f.write("(")
        f.write(str(x[ind]))
        f.write(" ")
        f.write(str(y[ind]))
        f.write(")\n")
    f.write(")\n")
    f.close()
    return x,y,v


"""Check the plotting along line function"""
def test_plotAlongLine():
    createFile()
    
    data = eulerianStatisticsData("tests/dataFile.dat")
    maxX = 0.5*max(data.pos()[0])
    maxY = 0.5*max(data.pos()[1])
    point1 = np.array([0.0, 0.0])
    point2 = np.array([maxX, maxY])
    x,y = data.plotAlongLine(point1, point2, 100)
    assert x[0] == 0
    assert x[99] == np.sqrt(maxX**2+maxY**2)

    # Check the y value
    t = np.linspace(0,1,100)

    for i in range(len(t)):
        P = point1 + t[i]*(point2-point1)
        assert (y[i] - func(P[0],P[1])) < 1E-5
    
