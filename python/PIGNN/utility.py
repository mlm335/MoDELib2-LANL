import numpy as np
import torch as th
import dgl
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN

class GraphConstruct():
    def __init__(self):
        super(GraphConstruct, self).__init__()
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
        output['isPointANode'] = np.amax(np.abs(output['stress']),axis=(1))==0
        output['isPointAGP'] = np.invert(output['isPointANode'])
        return output
    def constructgraph(self, pos, stress, burgers_vector, CELLS, T):
        data = self.readdata(pos, stress, burgers_vector, CELLS, T)

    
        nodes = np.arange(data['nPoints'])
        cells = np.arange(data['nCells'])
        n_nodes = data['isPointANode'].sum()
        n_GP = data['isPointAGP'].sum()
        n_cells = len(cells)
        src = []
        dst = []
        isrc = []
        idst = []

        snormal = np.array([0.0,1.0,0.0])
        features_burger = np.zeros([n_GP,1])
        features_stress = np.zeros([n_GP,6])
        features_stresst = np.zeros([n_GP,12])

        features_T = np.ones([n_nodes,1]) *data['T']
        features_dis = np.zeros([n_nodes,6])
        features_pos = data['points'][data['isPointANode']]
        for j in range(n_cells):
            n_nodes_in_cell = len(data['CELLS'][j])
            dx1 = data['points'][data['CELLS'][j][n_nodes_in_cell-1]]- data['points'][data['CELLS'][j][0]]
            burger_ref = data['burgers_vector_global'][0]
            if burger_ref[0] < 0.0:
                burger_ref = -1.0*burger_ref
            if np.abs(((data['burgers_vector_global'][j] - burger_ref)**2).sum()) < 0.0001 :
                features_dis[data['CELLS'][j][n_nodes_in_cell-1], :3] = dx1
                features_dis[data['CELLS'][j][0], 3:6] = dx1
                is_align = True
            else:
                features_dis[data['CELLS'][j][n_nodes_in_cell-1], 3:6] = -dx1
                features_dis[data['CELLS'][j][0], :3] = -dx1
                is_align = False
            for i in range(1,n_nodes_in_cell-1):
                dst.append(data['CELLS'][j][0])
                dst.append(data['CELLS'][j][n_nodes_in_cell-1])
                src.append(data['CELLS'][j][i]-n_nodes)
                src.append(data['CELLS'][j][i]-n_nodes)


                
                isrc.append(data['CELLS'][j][0])
                isrc.append(data['CELLS'][j][n_nodes_in_cell-1])
                idst.append(data['CELLS'][j][i]-n_nodes)
                idst.append(data['CELLS'][j][i]-n_nodes)
                b = data['burgers_vector_global'][j]
                b = b/(b*b).sum()**0.5


                dx1 = dx1/(dx1*dx1).sum()**0.5
                if not is_align:
                    b=-b
                    dx1g = -dx1
                else:
                    b=b
                    dx1g = dx1
                features_stresst[data['CELLS'][j][i]-n_nodes, :3] = b
                features_stresst[data['CELLS'][j][i]-n_nodes, 3:6] = dx1g
                features_stresst[data['CELLS'][j][i]-n_nodes, 6:] = data['stress'][data['CELLS'][j][i]][[0,1,2,4,5,8]]

                sign_theta = np.sign(np.dot(np.cross(b,dx1g),snormal))
                features_burger[data['CELLS'][j][i]-n_nodes] = sign_theta*(b*dx1g).sum()
                stress = data['stress'][data['CELLS'][j][i]]

                eigva, eigvec = np.linalg.eigh(stress.reshape([3,3]))
                if ((dx1g*eigvec[:,0]).sum()) < 0:
                    eigvec[:,0] = -eigvec[:,0]
                if ((dx1g*eigvec[:,1]).sum()) < 0:
                    eigvec[:,1] = -eigvec[:,1]

                if np.dot(np.cross(eigvec[:,0], eigvec[:,1]),eigvec[:,2]) < 0:
                    eigvec[:,2] = -eigvec[:,2]
                features_stress[data['CELLS'][j][i]-n_nodes, :3] = eigva
                features_stress[data['CELLS'][j][i]-n_nodes, 3] = ((dx1g*eigvec[:,0]).sum())
                features_stress[data['CELLS'][j][i]-n_nodes, 4] = ((dx1g*eigvec[:,1]).sum())
                features_stress[data['CELLS'][j][i]-n_nodes, 5] = ((dx1g*eigvec[:,2]).sum())


        features_mag = []
        features_tangent = []
        center  = data['points'][data['isPointANode']].mean(axis=0)
        for i in range(n_nodes):
            dx1 = (features_dis[i][:3]**2).sum()**0.5
            dx2 = (features_dis[i][3:6]**2).sum()**0.5
            if dx1 == 0:
                dx1 = dx2
                features_dis[i][:3] = features_dis[i][3:6]
            if dx2 == 0:
                dx2 = dx1
                features_dis[i][3:6] = features_dis[i][:3]
            tan1 = features_dis[i][:3]/dx1
            tan2 = features_dis[i][3:6]/dx2
            dx12 = (tan1*tan2).sum()
            orth1 = tan1  + tan2
            orth1[1] = 0.0
            orth2 = np.array([1.0,0.0,1.0])#data['points'][i] - center                                                                                                        
            orth1 = orth1/(orth1*orth1).sum()**0.5
            orth2 = orth2/(orth2*orth2).sum()**0.5
            orth2 = orth2 - (orth1*orth2).sum()*orth1
            orth2 = orth2/(orth2*orth2).sum()**0.5
            if np.cross(orth1, orth2).sum() < 0:
                orth2 = -orth2
            orth3 = np.cross(orth1, orth2)
            features_mag.append(np.array([dx1/1000,dx2/1000,dx12]))
            features_tangent.append(np.array([orth1, orth2, orth3]))



        features_mag = np.stack(features_mag)
        features_tangent = np.stack(features_tangent)

        dst = np.array(dst)
        src = np.array(src)
        
        idst = np.array(idst)
        isrc = np.array(isrc)

 
        features_GP = np.concatenate([ features_burger,features_stress], axis=1 )
        features_Node = features_mag#np.concatenate([ features_mag, features_T], axis=1)
        graph = dgl.heterograph({('GP','edge','Node'):(src,dst), ('Node','iedge','GP'):(isrc,idst)})
        graph=graph
        graph.nodes['GP'].data['GPfeat'] = th.from_numpy(features_GP)
        graph.nodes['GP'].data['bdott'] = th.from_numpy(features_burger)
        
        graph.nodes['GP'].data['GPfeatt'] = th.from_numpy(features_stresst)

        graph.nodes['Node'].data['Nodefeat'] = th.from_numpy(features_Node)
        graph.nodes['Node'].data['Nodetan'] = th.from_numpy(features_tangent)
        graph.nodes['Node'].data['pos'] = th.from_numpy(features_pos)
        graph.nodes['Node'].data['T'] = th.from_numpy(features_T)

        graph.nodes['GP'].data['GPfeat'][:,:1] = -graph.nodes['GP'].data['GPfeat'][:,:1]
        graph.nodes['GP'].data['bdott'][:,:1] = -graph.nodes['GP'].data['bdott'][:,:1]


        return graph


def VTKparser(filename):

    reader = vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()

    data = reader.GetOutput()

    output = {'nPoints':None,
              'nCells':None,
              'CELLS':None,
              'CELL_TYPES':None,
              'points':None,
              'isPointANode':None,
              'isPointAGP':None}

    pointFeatures = ['PK_force','stress']
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




        
    

        
        
    
