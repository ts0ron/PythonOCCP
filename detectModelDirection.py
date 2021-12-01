"""

Find the longest edge of the part - returns the direction vector and the length

@author: Ron Tsabary

"""

import random
import os
import os.path
import sys
import time
import numpy as np

from OCC.Core.BRepLib import BRepLib_MakeFace, BRepLib_MakeShell

from OCC.Extend.ShapeFactory import point_list_to_TColgp_Array1OfPnt, make_face
from OCC.Core.TopoDS import TopoDS_Wire, TopoDS_Face, TopoDS_Builder
from OCC.Extend.TopologyUtils import list_of_shapes_to_compound
from OCCUtils.Topology import (WireExplorer, Topo)
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity
# from OCC.Display.OCCViewer import
from OCC.Display.SimpleGui import init_display
from OCC.Core.BRep import BRep_Tool_Surface
from OCC.Core.GeomAbs import GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_BSplineSurface, GeomAbs_Cone, GeomAbs_Sphere, \
    GeomAbs_Torus, GeomAbs_BezierSurface, GeomAbs_SurfaceOfRevolution, GeomAbs_SurfaceOfExtrusion, \
    GeomAbs_OffsetSurface, GeomAbs_OtherSurface, GeomAbs_Circle
# from OCC.Core.TopoDS import topods_Face
from OCC.Core.TopoDS import topods_Vertex
from OCC.Core.TopoDS import topods_Face
from OCC.Core.TopoDS import topods_Edge
from OCC.Core.TopoDS import topods_Wire
from OCC.Core.TopoDS import topods_Shell
from OCC.Core.TopoDS import topods_Solid
from OCC.Core.TopoDS import topods_Compound
from OCC.Core.TopoDS import topods_CompSolid
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.TopExp import (TopExp_Explorer, topexp_MapShapesAndAncestors)
from OCC.Core.TopAbs import (TopAbs_VERTEX, TopAbs_EDGE, TopAbs_FACE, TopAbs_WIRE,
                             TopAbs_SHELL, TopAbs_SOLID, TopAbs_COMPOUND,
                             TopAbs_COMPSOLID)
# from OCC.TopExp import TopExp_Explorer, topexp_MapShapesAndAncestors
from OCC.Core.TopTools import (TopTools_ListOfShape,
                               TopTools_IndexedDataMapOfShapeListOfShape)
from OCC.Core.ShapeAnalysis import ShapeAnalysis_CheckSmallFace, ShapeAnalysis_Curve
from OCC.Core.ShapeAnalysis import shapeanalysis
from OCC.Core.BRepTools import breptools, breptools_UVBounds
from OCC.Core.TopOpeBRepTool import TopOpeBRepTool_ShapeTool
from OCC.Core.BOPTools import BOPTools_AlgoTools
from OCC.Core.BRep import BRep_Tool


def read_step_file(filename, as_compound=True):
    """ read the STEP file and returns a compound
    """
    step_reader = STEPControl_Reader()
    status = step_reader.ReadFile(filename)

    if status == IFSelect_RetDone:  # check status
        failsonly = False
        step_reader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity)
        step_reader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity)

        ok = step_reader.TransferRoot(1)
        _nbs = step_reader.NbShapes()
        if _nbs == 0:
            raise AssertionError("No shape to transfer.")
        elif _nbs == 1:  # most cases
            return step_reader.Shape(1)
        elif _nbs > 1:
            print("Number of shapes:", _nbs)
            shps = []
            # loop over root shapes
            for k in range(1, _nbs + 1):
                new_shp = step_reader.Shape(k)
                if not new_shp.IsNull():
                    shps.append(new_shp)
            if as_compound:
                compound, result = list_of_shapes_to_compound(shps)
                if not result:
                    print("Warning: all shapes were not added to the compound")
                return compound
            else:
                print("Warning, returns a list of shapes.")
                return shps
    else:
        raise AssertionError("Error: can't read file.")
    return None


def isEmptyList(list):
    """
    Boolean Function to check if a list is empty
    """
    if len(list) == 0:
        return True
    return False


def distance(p1, p2):
    """
    Calculates the distance between two gp_Pnts
    """
    x1, y1, z1 = p1.X(), p1.Y(), p1.Z()
    x2, y2, z2 = p2.X(), p2.Y(), p2.Z()
    x = x1 - x2
    y = y1 - y2
    z = z1 - z2
    return np.sqrt(x ** 2 + y ** 2 + z ** 2), [x, y, z]


if __name__ == '__main__':
    display, start_display, add_menu, add_function_to_menu = init_display()
    display.SetSelectionModeFace()
    # BJB_Terminal_block_3_poles_46.413.1115.50.STEP
    # threads_cubic.STEP
    # 543-03-A001-09_assembled_panel.STEP
    # reading the shape
    source_dir = "C:\\Users\\Administrator\\Desktop\\Python-Occ-Algo-Files\\files-dir\\"
    shp = read_step_file(source_dir + "BJB_Terminal_block_3_poles_46.413.1115.50.STEP")

    p = BOPTools_AlgoTools()
    g = BRep_Tool()
    t = Topo(shp)

    faceExplorer = TopExp_Explorer(shp, TopAbs_FACE)

    faces = []

    while faceExplorer.More():
        curr_face = faceExplorer.Current()
        faces.append(curr_face)
        faceExplorer.Next()

    st = TopOpeBRepTool_ShapeTool()
    max_length = 0
    main_direction = [0, 0, 0]
    max_p1 = None
    max_p2 = None
    for f in faces:
        brep_t = breptools()
        outer_wire = brep_t.OuterWire(f)
        w_exp = WireExplorer(outer_wire)
        edge_exp = w_exp.ordered_edges()

        while True:
            """
            Entering the loop to find the longest edge
            and calculate her length and direction
            """
            try:
                # next : returns an edge if the explorer is not empty, else raises an exception
                edge = edge_exp.__next__()
                e_topo = Topo(edge)
                ver_exp = e_topo.vertices()
                vertices = []
                while True:
                    # Extracting the vertices of an edge for calculation
                    try:
                        vertex = ver_exp.__next__()
                        vertices.append(vertex)
                    except Exception as e:
                        break

                # Converting the verices to a geometric point
                p1 = st.Pnt(vertices[0])
                p2 = st.Pnt(vertices[1])
                # Caclulating the length of the edge and setting a new maximum length if needed
                dist, direction = distance(p1, p2)
                if dist > max_length:
                    max_length = dist
                    max_direction = direction
                    max_p1 = p1
                    max_p2 = p2
            except Exception as e:
                break

    print("########################################\n##############################\n############\n######\n##")
    print("Finished running on all edges")
    print("First point: X-", max_p1.X(), " Y -", max_p1.Y(), "Z -", max_p1.Z())
    print("Second point: X-", max_p2.X(), " Y -", max_p2.Y(), "Z -", max_p2.Z())
    print("Direction of the edge - ", max_direction / max_length)    # printing a normalized direction
    print("Distance of the edge - ", max_length)

    """
    Vertex to point - The ShapeTool is the way to convert from an topo-entity to a geometric object
    """
    # point = st.Pnt(vertices.__next__())
    # print(point.X())
    # print(point.Y())
    # print(point.Z())

    display.DisplayShape(shp, update=True)
    start_display()
