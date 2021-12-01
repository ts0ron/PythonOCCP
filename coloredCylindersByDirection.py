"""

Extract cylinders from the assembly - in colors for viisualization
Checking the type of each face - for a cylinderical type
The name of the file s- in main
need to enter direction of cylinder
Shows just cylinder

@author: Ron Tsabary

"""

import random
import os
import os.path
import sys
import time
import numpy as np
from OCCUtils.Topology import Topo
from OCC.Core.BRepLib import BRepLib_MakeFace, BRepLib_MakeShell
from OCC.Core.Quantity import Quantity_Color, Quantity_TOC_RGB
from OCC.Extend.ShapeFactory import point_list_to_TColgp_Array1OfPnt, make_face
from OCC.Core.TopoDS import TopoDS_Wire, TopoDS_Face, TopoDS_Builder, TopoDS_Shape
from OCC.Extend.TopologyUtils import list_of_shapes_to_compound
from OCC.Extend.TopologyUtils import WireExplorer
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
from OCC.Core.BOPTools import BOPTools_AlgoTools
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
# from OCC.Geom import Handle_Geom_Curve
from OCC.Core.GeomAdaptor import GeomAdaptor_Curve
from OCC.Core.GCPnts import GCPnts_QuasiUniformDeflection, GCPnts_UniformDeflection, GCPnts_UniformAbscissa
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop_LinearProperties
from OCC.Core.gp import gp_Pnt, gp_Cylinder
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse, BRepAlgoAPI_Section
from OCC.Core.Geom import Geom_Line
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace
from OCC.Core.IntCurvesFace import IntCurvesFace_ShapeIntersector
from OCC.Core.ShapeExtend import ShapeExtend_WireData

# Precision module allows us to access const variables such as Confusion, Angular etc. (epsilon of distance or angle respectively)
from OCC.Core.Precision import precision_Confusion


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


def recognize_face(a_face):
    """ Takes a TopoDS shape and tries to identify its nature
    whether it is a plane a cylinder a torus etc.
    if a plane, returns the normal
    if a cylinder, returns the radius
    """
    # turn the face (inifinte area) into surface (finite area)
    surf = BRepAdaptor_Surface(a_face, True)
    surf_type = surf.GetType()
    if surf_type == GeomAbs_Plane:
        # print("--> plane")
        # look for the properties of the plane
        # first get the related gp_Pln
        gp_pln = surf.Plane()
        location = gp_pln.Location()  # a point of the plane
        normal = gp_pln.Axis().Direction()  # the plane normal
        # then export location and normal to the console output
        # print("--> Location (global coordinates)", location.X(), location.Y(), location.Z())
        # print("--> Normal (global coordinates)", normal.X(), normal.Y(), normal.Z())

    elif surf_type == GeomAbs_Cylinder:
        # print("--> cylinder")
        # look for the properties of the cylinder
        # first get the related gp_Cyl
        gp_cyl = surf.Cylinder()
        location = gp_cyl.Location()  # a point of the axis
        axis = gp_cyl.Axis().Direction()  # the cylinder axis
        # then export location and normal to the console output
        # print("--> Location (global coordinates)", location.X(), location.Y(), location.Z())
        # print("--> Axis (global coordinates)", axis.X(), axis.Y(), axis.Z())
        return True, location, axis

    elif surf_type == GeomAbs_BSplineSurface:
        # print("--> BSplineSurface")
        gp_bsrf = surf.Surface()
        # degree = gp_bsrf.NbUKnots()

    elif surf_type == GeomAbs_Cone:
        kind = "Cone"
        # print(kind)
        gp_bsrf = surf.Surface()
        # degree = gp_bsrf.NbUKnots()

    elif surf_type == GeomAbs_Sphere:
        kind = "Sphere"
        # print(kind)
        gp_bsrf = surf.Surface()
        # degree = gp_bsrf.NbUKnots()
    elif surf_type == GeomAbs_Torus:
        # print(a_face)
        kind = "Torus"
        return kind, None, None
    elif surf_type == GeomAbs_BezierSurface:
        # print(a_face)
        kind = "Bezier"
        return kind, None, None
    elif surf_type == GeomAbs_SurfaceOfRevolution:
        # print(a_face)
        kind = "Revolution"
        return kind, None, None
    elif surf_type == GeomAbs_SurfaceOfExtrusion:
        # print(a_face)
        kind = "Extrusion"
        return kind, None, None
    elif surf_type == GeomAbs_OffsetSurface:
        # print(a_face)
        kind = "Offset"
        return kind, None, None
    elif surf_type == GeomAbs_OtherSurface:
        # print(a_face)
        kind = "Other"
    else:
        print()
        # TODO - there are plenty other type that can be checked
        # see documentation for the BRepAdaptor class
        # https://www.opencascade.com/doc/occt-6.9.1/refman/html/class_b_rep_adaptor___surface.html
        # print(surf_type)
        # print("not implemented")
    return False, None, None


def isPlane(aFace):
    surf = BRepAdaptor_Surface(aFace, True)
    surf_type = surf.GetType()

    if surf_type == GeomAbs_Plane:
        return True
    return False


def getPlaneInfo(aFace):
    surf = BRepAdaptor_Surface(aFace, True)
    surf_type = surf.GetType()
    if surf_type == GeomAbs_Plane:
        print("--> plane")
        # look for the properties of the plane
        # first get the related gp_Pln
        gp_pln = surf.Plane()
        location = gp_pln.Location()  # a point of the plane
        normal = gp_pln.Axis().Direction()  # the plane normal
        # then export location and normal to the console output
        print("--> Location (global coordinates)", location.X(), location.Y(), location.Z())
        print("--> Normal (global coordinates)", normal.X(), normal.Y(), normal.Z())
        return location, normal
    return False


def isEmptyList(lst):
    if len(lst) == 0:
        return True
    return False


def equalNormals(normal_a, normal_b):
    if normal_a.X() == normal_b.X() and normal_a.Y() == normal_b.Y() and normal_a.Z() == normal_b.Z():
        return True
    return False


def is_face_explored(face, explored):
    for base in explored:
        for f in base:
            if TopoDS_Face.IsEqual(face, f):
                return True
    return False


def fuse_walls(walls):
    """
    This function recieves a list of walls and combine it all to 1 Shell.
    @input: walls - list of faces
    @returns: TopoDS_Shell
    """
    if isEmptyList(walls):
        return None
    wall_shell = BRepLib_MakeShell(walls[0])
    shell = wall_shell.Shell()
    walls.remove(walls[0])
    for x in walls:  # Adding all the suspected walls to be one shape before making an intersection
        TopoDS_Builder.Add(shell, x)
    return shp


def is_base_of_channel(base, walls):
    shp = fuse_walls(walls)


def base_builder(base, base_to_explore, walls):
    """
    Builds the base from the given base_to_explore

    @returns: Dictionary of the kind :  { is_channel: Boolean, Base: TopoDS_Face[] }
    """
    if isEmptyList(base_to_explore) == 0:
        return is_base_of_channel(base, walls)  # TODO - the intersection thingy and closed wire validation
    else:
        curr_face = base_to_explore[0]
        base.append(curr_face)
        base_to_explore.remove(curr_face)

        base_location, base_normal = getPlaneInfo(curr_face)

        face_edges = Topo.edges_from_face(curr_face)
        for e in face_edges:
            faces_connected_to_edge = Topo.faces_from_edge(e)
            for f in faces_connected_to_edge:
                if isPlane(f):
                    if TopoDS_Face.IsEqual(f, curr_face) or (f in base):  # If the current checked face is
                        continue  # curr_face or it is in the base list, skip.
                    if isEmptyList(base):
                        base.append(f)
                        continue
                    curr_location, curr_normal = getPlaneInfo(f)
                    if equalNormals(base_normal, curr_normal):
                        base_to_explore.append(f)
                    else:
                        walls.append(f)
                    base_builder(base, base_to_explore, walls)


if __name__ == '__main__':
    # arrays for internal wires information.
    # arrays for cylindrical faces information
    source_dir = "C:\\Users\\Administrator\\Desktop\\Python-Occ-Algo-Files\\files-dir\\"
    channel_bases = []  # saves each channle's base plan
    explored_bases = []
    display, start_display, add_menu, add_function_to_menu = init_display()
    display.SetSelectionModeFace()

    # reading the shape
    #threads_cubic.STEP
    shp = read_step_file(
       source_dir + "BJB_Terminal_block_3_poles_46.413.1115.50.STEP")

    p = BOPTools_AlgoTools()
    g = BRep_Tool()
    t = Topo(shp)
    t_shape = TopoDS_Shape()
    w_t = TopoDS_Wire
    f_t = TopoDS_Face

    """
    Identifying any Cylindrical Faces and displaying them with random color each time.
    """
    input_direction = input("Please state one of the option: x, y, z for the direction filter in lower case!!!!")

    faceExplorer = TopExp_Explorer(shp, TopAbs_FACE)
    while faceExplorer.More():
        curr_face = faceExplorer.Current()
        ans1, location1, axis1 = recognize_face(curr_face)
        if ans1:
            if not location1 is None and not axis1 is None:
                # print("Location: ", location1.X(), location1.Y(), location1.Z())
                # print("Axis: ", axis1.X(), axis1.Y(), axis1.Z())
                if input_direction == 'x' and np.absolute(axis1.X()) == 1 and (np.absolute(axis1.Y()) < 1e-13) and (np.absolute(axis1.Z()) < 1e-13):
                    color = Quantity_Color(random.random(),
                                           random.random(),
                                           random.random(),
                                           Quantity_TOC_RGB)
                    display.DisplayColoredShape(curr_face, color)
                    print("################################################################## Found a Cylinder:::::::")
                elif input_direction == 'y' and np.absolute(axis1.Y()) == 1 and (np.absolute(axis1.X()) < 1e-13) and (np.absolute(axis1.Z()) < 1e-13):
                    color = Quantity_Color(random.random(),
                                           random.random(),
                                           random.random(),
                                           Quantity_TOC_RGB)
                    display.DisplayColoredShape(curr_face, color)
                    print("################################################################## Found a Cylinder:::::::")
                elif input_direction == 'z' and np.absolute(axis1.Z()) == 1 and (np.absolute(axis1.Y()) < 1e-13) and (np.absolute(axis1.X()) < 1e-13):
                    color = Quantity_Color(random.random(),
                                           random.random(),
                                           random.random(),
                                           Quantity_TOC_RGB)
                    display.DisplayColoredShape(curr_face, color)
                    print("################################################################## Found a Cylinder:::::::")
        faceExplorer.Next()
    start_display()
