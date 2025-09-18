import sys, string, os
import numpy as np

class QP:
    sourceID=np.empty(0)
    sinkID=np.empty(0)
    loopID=np.empty(0)
    position=np.empty([0,0])
    j=np.empty(0)
    tangent=np.empty([0,0])
    stress=np.empty([0,0])
    pkForce=np.empty([0,0])
    stackingFaultForce=np.empty([0,0])
    lineTensionForce=np.empty([0,0])
    velocity=np.empty([0,0])
    elasticEnergyPerLength=np.empty([0,0])
    coreEnergyPerLength=np.empty([0,0])
    cCD=np.empty([0,0])
    cDD=np.empty([0,0])

def readAUXtxt(filename):
    try:
        with open(filename+'.txt', "r") as auxFile:
            numMeshNodes=int(auxFile.readline().rstrip())
            numQP=int(auxFile.readline().rstrip())
            numBoundary=int(auxFile.readline().rstrip())
        
            qp=QP();
            
            qp.position=np.empty([numQP, 3])
            qp.j=np.empty([numQP, 3])
            qp.tangent=np.empty([numQP, 3])
            qp.velocity=np.empty([numQP, 3])

            # node data
            for k in range(numQP):
                data = np.fromstring(auxFile.readline().strip(), sep=' ')
                qp.position[k, :] = data[3:6]
                qp.j[k, :] = data[6]
                qp.tangent[k, :] = data[7:10]
                qp.velocity[k, :] = data[28:31]

        return qp

    except Exception as e:
        print(f"Error reading file {filename}: {e}")
    return None
