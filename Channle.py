"""
The algorithm for finding a channel

Collect for a plannar base the adjacent strutures with the same normal and suspected walls (adjacent surfaces with a different normal)-

to finish
Create a compound from a list of shapes (the walls)

Compute intersection



@author: Ron Tsabary
"""


import random
import os
import os.path
import sys
import time

from OCC.Core.BRepLib import BRepLib_MakeFace, BRepLib_MakeShell
from OCC.Core.TopTools import TopTools_IndexedDataMapOfShapeListOfShape, TopTools_ListOfShape
from OCC.Extend.ShapeFactory import point_list_to_TColgp_Array1OfPnt, make_face
from OCC.Core.TopoDS import TopoDS_Wire, TopoDS_Face, TopoDS_Builder
from OCC.Extend.TopologyUtils import list_of_shapes_to_compound
from OCC.Core.TopoDS import topods_Vertex, topods_Face, topods_Edge, topods_Wire, topods_Shell, topods_Compound, \
    topods_Solid, topods_CompSolid, TopoDS_Iterator
from OCC.Core.TopExp import topexp_MapShapesAndAncestors
from OCCUtils.Topology import Topo, WireExplorer
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity
from OCC.Display.SimpleGui import init_display
from OCC.Core.GeomAbs import GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_BSplineSurface, GeomAbs_Cone, GeomAbs_Sphere, GeomAbs_Torus, GeomAbs_BezierSurface, GeomAbs_SurfaceOfRevolution, GeomAbs_SurfaceOfExtrusion, GeomAbs_OffsetSurface, GeomAbs_OtherSurface, GeomAbs_Circle
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.TopExp import (TopExp_Explorer, topexp_MapShapesAndAncestors)
from OCC.Core.TopAbs import (TopAbs_VERTEX, TopAbs_EDGE, TopAbs_FACE, TopAbs_WIRE,
                        TopAbs_SHELL, TopAbs_SOLID, TopAbs_COMPOUND,
                        TopAbs_COMPSOLID)
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
        print("--> plane")
        # look for the properties of the plane
        # first get the related gp_Pln
        gp_pln = surf.Plane()
        location = gp_pln.Location()  # a point of the plane
        normal = gp_pln.Axis().Direction()  # the plane normal
        # then export location and normal to the console output
        print("--> Location (global coordinates)", location.X(), location.Y(), location.Z())
        print("--> Normal (global coordinates)", normal.X(), normal.Y(), normal.Z())

    elif surf_type == GeomAbs_Cylinder:
        print("--> cylinder")
        # look for the properties of the cylinder
        # first get the related gp_Cyl
        gp_cyl = surf.Cylinder()
        location = gp_cyl.Location()  # a point of the axis
        axis = gp_cyl.Axis().Direction()  # the cylinder axis
        # then export location and normal to the console output
        print("--> Location (global coordinates)", location.X(), location.Y(), location.Z())
        print("--> Axis (global coordinates)", axis.X(), axis.Y(), axis.Z())
        return True, location, axis

    elif surf_type == GeomAbs_BSplineSurface:
        print("--> BSplineSurface")
        gp_bsrf = surf.Surface()
        # degree = gp_bsrf.NbUKnots()

    elif surf_type == GeomAbs_Cone:
        kind = "Cone"
        print(kind)
        gp_bsrf = surf.Surface()
        # degree = gp_bsrf.NbUKnots()

    elif surf_type == GeomAbs_Sphere:
        kind = "Sphere"
        print(kind)
        gp_bsrf = surf.Surface()
        # degree = gp_bsrf.NbUKnots()
    elif surf_type == GeomAbs_Torus:
        print(a_face)
        kind = "Torus"
        return kind, None, None
    elif surf_type == GeomAbs_BezierSurface:
        print(a_face)
        kind = "Bezier"
        return kind, None, None
    elif surf_type == GeomAbs_SurfaceOfRevolution:
        print(a_face)
        kind = "Revolution"
        return kind, None, None
    elif surf_type == GeomAbs_SurfaceOfExtrusion:
        print(a_face)
        kind = "Extrusion"
        return kind, None, None
    elif surf_type == GeomAbs_OffsetSurface:
        print(a_face)
        kind = "Offset"
        return kind, None, None
    elif surf_type == GeomAbs_OtherSurface:
        print(a_face)
        kind = "Other"
    else:
        print(a_face)
        # TODO - there are plenty other type that can be checked
        # see documentation for the BRepAdaptor class
        # https://www.opencascade.com/doc/occt-6.9.1/refman/html/class_b_rep_adaptor___surface.html
        print(surf_type)
        print("not implemented")
    return False, None, None


def isPlane(aFace):
    """
    Check if a face is a plane
    """
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
    Attempt to fuse the walls - was unsuccessful - TODO
    Recommendation is to check for the list_of_shapes_to_compound function for connecting the walls into
    one entity.
    """
    if isEmptyList(walls):
        return None
    wall_shell = BRepLib_MakeShell(walls[0])
    shell = wall_shell.Shell()
    walls.remove(walls[0])
    for x in walls:     # Adding all the suspected walls to be one shape before making an intersection
        TopoDS_Builder.Add(shell, x)
    return shp


def make_plane_for_intersection(base):
    """
    This function should take a face in the base (randomly - or just the first one)
    and returns an infinite plane with normal indentical to the normal of the base_face
    and it's location is ( base_face.location + 1e-5 * base_face.normal )
    meaning an infinite plane move by 'epsilon' from the base.
    """
    raise NotImplementedError


def intersect(wall, intersection_plane):
    """
    This function will take the wall and the intersection plane and return the intersection
    Notes - the intersection is a wire\edges\curves - need to check the type returned and
    build a wire from it.
    return: that wire.
    """
    raise NotImplementedError

def validate_is_channel(intersected_wire):
    """
    This function will take the intersected wire and check if it is a closed wire
    if it is, will return True and the base related to the intersected wire is a base of a channel
    if it is not, will return False and the base related to the intersected wire is an open wall.
    """
    raise NotImplementedError


def base_builder(base, base_to_explore, walls):
    """
    Builds the base from the given base_to_explore

    @base_to_explore: is a list of faces that are suspected to be a base of a channel that we haven't explored yet
    @base: is a list of faces that are already a base of a channel that we have explored already
    @walls: is a list of faces that are connected to a base face and has a different normal

    @returns: Dictionary of the kind :  {Base: TopoDS_Face[], walls:TopoDS_Face[] }
    """
    if isEmptyList(base_to_explore):
        return {"base": base,
                "walls": walls}
    else:
        # Extracting a face from the base_to_explore list and adding it to the base list
        curr_face = base_to_explore[0]
        base.append(curr_face)
        base_to_explore.remove(curr_face)

        # Getting the data of the current base face we have extracted
        base_location, base_normal = getPlaneInfo(curr_face)

        # Getting the edges of the curr_face we're exploring
        face_edges = Topo.edges_from_face(curr_face)
        for e in face_edges:
            # Getting the faces that are connected to the current edge
            faces_connected_to_edge = Topo.faces_from_edge(e)
            for f in faces_connected_to_edge:
                """
                for each face connected to the edge, check if it is a plane and if so, catagorize it by its normal.
                Catogarize types - base_to_explore, walls, or not relavent. 
                """
                if isPlane(f):
                    if TopoDS_Face.IsEqual(f, curr_face) or (f in base):   # If the current checked face is
                        continue                                            # curr_face or it is in the base list, skip.
                    if isEmptyList(base):
                        base.append(f)
                        continue
                    curr_location, curr_normal = getPlaneInfo(f)
                    if equalNormals(base_normal, curr_normal):
                        base_to_explore.append(f)
                    else:
                        walls.append(f)
                    return base_builder(base, base_to_explore, walls) # the recursive call



if __name__ == '__main__':
    # arrays for internal wires information.
    # arrays for cylindrical faces information
    channel_bases = []          # saves each channle's base plan
    explored_bases = []
    display, start_display, add_menu, add_function_to_menu = init_display()
    display.SetSelectionModeFace()

    # reading the shape
    shp = read_step_file("/step_files_dir/threads_cubic.STEP")

    p = BOPTools_AlgoTools()
    g = BRep_Tool()
    t = Topo(shp)
    w_t = TopoDS_Wire
    f_t = TopoDS_Face

    faceExplorer = TopExp_Explorer(shp, TopAbs_FACE)

    while faceExplorer.More():
        # TODO
        """
        While there are face to explore - check if we have encountered that face before in the plane_walk
        If so, continue, else build base and walls
        """
        curr_face = faceExplorer.Current()
        ans, location, axis = recognize_face(curr_face)
        if is_face_explored(curr_face, explored_bases):
            faceExplorer.Next()
            continue
        base_builder([], [faceExplorer.Current()], [])
        break

    display.DisplayShape(shp, update=True)
    start_display()




