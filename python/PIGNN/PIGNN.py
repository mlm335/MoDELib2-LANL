import numpy as np
import torch as th
import dgl
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN

class PIGNN():
    def __init__(self):
        super(PIGNN, self).__init__()
        # create model
        # load model from memory

    def readdata(self, pos, stress, burgers_vector, CELLS, T):

### Read data from DDD, assuming number of nodes and GPs are n_node, n_GP, and n_seg number of segments
###  pos: double [n_node+n_GP, 3]
###  stress: double [n_node+n_GP, 3,3]
###  burgers_vector: double [n_seg, 3]
###  CELSS: list of lists of int, len(list) = n_seg
###  T: Double
        output = {'nPoints':None,
                  'nCells':None,
                  'CELLS':None,
                  'points':None,
                  'isPointANode':None,
                  'isPointAGP':None,
                  'T':None,
                  'stress':None}
        output['points'] = pos
        output['nPoints'] = pos.shape[0]
        output['CELLS'] = CELLS
        output['nCells'] = len(CELLS)
        output['stress'] =  stress
        output['T'] = T 
        output['burgers_vector_global'] = burgers_vector
        output['isPointANode'] = np.amax(np.abs(output['stress']),axis=(1,2))==0
        output['isPointAGP'] = np.invert(output['isPointANode'])
        return output


def VTKparser(filename):

    reader = vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()

    data = reader.GetOutput()
VTKparser(filename)
    output = {'nPoints':None,
              'nCells':None,
              'CELLS':None,
              'CELL_TYPES':None,
              'points':None,
              'isPointANode':None,
              'isPointAGP':None}

    pointFeatures = ['PK_force','velocity','stress']
    cellFeatures = ['dislocation_index','burgers_vector_global','burgers_vector_local']

    for feature in pointFeatures:
        output[feature] = None

    for feature in cellFeatures:
        output[feature] = None

    output['nPoints'] = data.GetNumberOfPoints()
    output['nCells'] = data.GetNumberOfCells()
    output['points'] = VN.vtk_to_numpy(data.GetPoints().GetData())
    output['CELL_TYPES']=VN.vtk_to_numpy(data.GetCellTypesArray())

    cells = []

    for cellID in range(output['nCells']):

        bufferNPoints = data.GetCell(cellID).GetNumberOfPoints()
        cells.append(np.array([ data.GetCell(cellID).GetPointId(pointID) for pointID in range(bufferNPoints)]))

    output['CELLS'] = cells

    for feature in pointFeatures:
        output[feature] = VN.vtk_to_numpy(data.GetPointData().GetArray(feature))

    for feature in cellFeatures:
        output[feature] = VN.vtk_to_numpy(data.GetCellData().GetArray(feature))

    output['isPointANode'] = np.amax(np.abs(output['PK_force']),axis=1)==0
    output['isPointAGP'] = np.invert(output['isPointANode'])
    return output



gnn = PIGNN()

pos = np.random.randn(100,3)
stress = np.random.randn(100,3,3)
burgers_vector_global = np.random.randn(20,3)
CELLS = [[1,2,3],[2,3,4]]
T = 300

out = gnn.readdata(pos, stress, burgers_vector_global, CELLS, T)

        
    

        
        
    
