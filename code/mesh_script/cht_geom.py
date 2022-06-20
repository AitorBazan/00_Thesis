from ofblockmeshdicthelper import BlockMeshDict, Vertex, Point, SimpleGrading
import numpy as np
import os

from tkinter import Tk
from tkinter.filedialog import askdirectory
path = askdirectory(title='Select Folder') 

bmd = BlockMeshDict()
r = 41.4
Hp = 1
Nx_pipe = 300
Ny_pipe = 300
Nz_pipe = 1
r2 = r + 10
rati = r/r2


bmd.add_vertex(-r/2,-r/2,-Hp,'v0')
bmd.add_vertex(0,-r/2,-Hp,'v1')
bmd.add_vertex(0,r/2,-Hp,'v2')
bmd.add_vertex(-r/2,r/2,-Hp,'v3')
bmd.add_vertex(-r/2,-r/2,0,'v4')
bmd.add_vertex(0,-r/2,0,'v5')
bmd.add_vertex(0,r/2,0,'v6')
bmd.add_vertex(-r/2,r/2,0,'v7')
bmd.add_vertex(-r*np.sin(np.pi/4),-r*np.cos(np.pi/4),-Hp,'v8')
bmd.add_vertex(0,-r,-Hp,'v9')
bmd.add_vertex(-r*np.sin(np.pi/4),r*np.cos(np.pi/4),-Hp,'v10')
bmd.add_vertex(0,r,-Hp,'v11')
bmd.add_vertex(-r*np.sin(np.pi/4),-r*np.cos(np.pi/4),0,'v12')
bmd.add_vertex(0,-r,0,'v13')
bmd.add_vertex(-r*np.sin(np.pi/4),r*np.cos(np.pi/4),0,'v14')
bmd.add_vertex(0,r,0,'v15')



bmd.add_vertex(-r2*np.sin(np.pi/4),-r2*np.cos(np.pi/4),-Hp,'v16')
bmd.add_vertex(0,-r2,-Hp,'v17')
bmd.add_vertex(-r2*np.sin(np.pi/4),r2*np.cos(np.pi/4),-Hp,'v18')
bmd.add_vertex(0,r2,-Hp,'v19')
bmd.add_vertex(-r2*np.sin(np.pi/4),-r2*np.cos(np.pi/4),0,'v20')
bmd.add_vertex(0,-r2,0,'v21')
bmd.add_vertex(-r2*np.sin(np.pi/4),r2*np.cos(np.pi/4),0,'v22')
bmd.add_vertex(0,r2,0,'v23')
bmd.add_arcedge(('v8','v9'),'arc1',Vertex(-r*np.sin(np.pi/8),-r*np.cos(np.pi/8),-Hp,'v_arc1'))
bmd.add_arcedge(('v12','v13'),'arc2',Vertex(-r*np.sin(np.pi/8),-r*np.cos(np.pi/8),0,'v_arc2'))
bmd.add_arcedge(('v8','v10'),'arc3',Vertex(-r,0,-Hp,'v_arc3'))
bmd.add_arcedge(('v12','v14'),'arc4',Vertex(-r,0,0,'v_arc4'))
bmd.add_arcedge(('v10','v11'),'arc5',Vertex(-r*np.sin(np.pi/8),r*np.cos(np.pi/8),-Hp,'v_arc5'))
bmd.add_arcedge(('v14','v15'),'arc6',Vertex(-r*np.sin(np.pi/8),r*np.cos(np.pi/8),0,'v_arc6'))
bmd.add_arcedge(('v16','v17'),'arc7',Vertex(-r2*np.sin(np.pi/8),-r2*np.cos(np.pi/8),-Hp,'v_arc7'))
bmd.add_arcedge(('v20','v21'),'arc8',Vertex(-r2*np.sin(np.pi/8),-r2*np.cos(np.pi/8),0,'v_arc8'))
bmd.add_arcedge(('v16','v18'),'arc9',Vertex(-r2,0,-Hp,'v_arc9'))
bmd.add_arcedge(('v20','v22'),'arc10',Vertex(-r2,0,0,'v_arc10'))
bmd.add_arcedge(('v18','v19'),'arc11',Vertex(-r2*np.sin(np.pi/8),r2*np.cos(np.pi/8),-Hp,'v_arc11'))
bmd.add_arcedge(('v22','v23'),'arc12',Vertex(-r2*np.sin(np.pi/8),r2*np.cos(np.pi/8),0,'v_arc12'))


prism_pipe = bmd.add_hexblock(('v0','v1','v2','v3','v4','v5','v6','v7'), (int(Nx_pipe/2),int(Ny_pipe/2),Nz_pipe),'prims_pipe',
                 grading=SimpleGrading(1,1,1))
south_pipe = bmd.add_hexblock(('v8','v9','v1','v0','v12','v13','v5','v4'),(int(Nx_pipe/2),int(Ny_pipe/4),Nz_pipe),'south_pipe',
                 grading=SimpleGrading(1,1,1))
east_pipe = bmd.add_hexblock(('v8','v0','v3','v10','v12','v4','v7','v14'),(int(Ny_pipe/4),int(Ny_pipe/2),Nz_pipe),'east_pipe',
                 grading=SimpleGrading(1,1,1))
north_pipe = bmd.add_hexblock(('v3','v2','v11','v10','v7','v6','v15','v14'),(int(Nx_pipe/2),int(Ny_pipe/4),Nz_pipe),'north_pipe',
                 grading=SimpleGrading(1,1,1))

#Region #2
south_pipe_region2 = bmd.add_hexblock(('v16','v17','v9','v8','v20','v21','v13','v12'),(int(Nx_pipe/2),int(Ny_pipe/4),Nz_pipe),'south_pipe_region2',
                 grading=SimpleGrading(1,1,1))
east_pipe_region2 = bmd.add_hexblock(('v16','v8','v10','v18','v20','v12','v14','v22'),(int(Ny_pipe/4),int(Ny_pipe/2),Nz_pipe),'east_pipe_region2',
                 grading=SimpleGrading(1,1,1))
north_pipe_region2 = bmd.add_hexblock(('v10','v11','v19','v18','v14','v15','v23','v22'),(int(Nx_pipe/2),int(Ny_pipe/4),Nz_pipe),'north_pipe_region2',
                 grading=SimpleGrading(1,1,1))

bmd.add_boundary('wall','solidWall',[south_pipe_region2.face('s'),east_pipe_region2.face('w'),north_pipe_region2.face('n')])
bmd.add_boundary('empty','fluidFrontAndBack',[prism_pipe.face('b'),south_pipe.face('b'),east_pipe.face('b'),north_pipe.face('b'),prism_pipe.face('t'),south_pipe.face('t'),east_pipe.face('t'),north_pipe.face('t')])
bmd.add_boundary('empty','solidFrontAndBack',[south_pipe_region2.face('b'),east_pipe_region2.face('b'),north_pipe_region2.face('b'),south_pipe_region2.face('t'),east_pipe_region2.face('t'),north_pipe_region2.face('t')])
bmd.add_boundary('symmetryPlane','fluidSymmetryBC',[prism_pipe.face('e'),south_pipe.face('e'),north_pipe.face('e')])
bmd.add_boundary('symmetryPlane','solidSymmetryBC',[south_pipe_region2.face('e'),north_pipe_region2.face('e')])


bmd.set_metric('mm')
bmd.assign_vertexid()
print(bmd.format())

print(bmd.format())

filetxt = bmd.format()
nameOfFile = 'blockMeshDict'
completeName = os.path.join(path, nameOfFile)
with open(completeName,'w') as f:
    f.writelines(filetxt)

f.close()