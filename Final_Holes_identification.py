# -*- coding: utf-8 -*-
"""
Find the holes based on two methods (evaluate holes (checks adjacent cylinders) and internal wire (checks center of internal wire))
Returns Json file with a dictionary with two fields - cylinder data center of wire - for the holes found.
Prompt a question - how to visualize results - wireframe of solid (full part)

Created on Wed Jan 27 18:44:02 2021

@author: Zofit Allouche
"""

import json
import random
import os
import os.path
import sys
import numpy as np
from OCCUtils.Topology import Topo, WireExplorer
from OCC.Core.Quantity import Quantity_Color
from OCC.Core.Quantity import Quantity_TOC_RGB
from OCC.Extend.TopologyUtils import list_of_shapes_to_compound
from OCC.Extend.TopologyUtils import WireExplorer
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity
from OCC.Display.SimpleGui import init_display
from OCC.Core.BRep import BRep_Tool_Surface
from OCC.Core.GeomAbs import GeomAbs_Plane, GeomAbs_Cylinder, GeomAbs_BSplineSurface, GeomAbs_Cone, GeomAbs_Sphere, \
    GeomAbs_Torus, GeomAbs_BezierSurface, GeomAbs_SurfaceOfRevolution, GeomAbs_SurfaceOfExtrusion, \
    GeomAbs_OffsetSurface, GeomAbs_OtherSurface, GeomAbs_Circle
from OCC.Core.TopoDS import topods_Vertex, topods_Face, topods_Edge, topods_Wire, topods_Shell, topods_Compound, \
    topods_Solid, topods_CompSolid, TopoDS_Iterator
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import (TopAbs_VERTEX, TopAbs_EDGE, TopAbs_FACE, TopAbs_WIRE,
                             TopAbs_SHELL, TopAbs_SOLID, TopAbs_COMPOUND,
                             TopAbs_COMPSOLID)
from OCC.Core.TopTools import (TopTools_ListOfShape,
                               TopTools_IndexedDataMapOfShapeListOfShape)
from OCC.Core.ShapeAnalysis import ShapeAnalysis_CheckSmallFace, ShapeAnalysis_Curve, shapeanalysis
from OCC.Core.BRepTools import breptools, breptools_UVBounds
from OCC.Core.BOPTools import BOPTools_AlgoTools
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
from OCC.Core.GeomAdaptor import GeomAdaptor_Curve
from OCC.Core.GCPnts import GCPnts_QuasiUniformDeflection, GCPnts_UniformDeflection, GCPnts_UniformAbscissa
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRepGProp import brepgprop_LinearProperties
from OCC.Core.gp import gp_Pnt, gp_Cylinder
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse, BRepAlgoAPI_Section
from OCC.Core.Geom import Geom_Line
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace
from OCC.Core.IntCurvesFace import IntCurvesFace_ShapeIntersector
from OCC.Core.ShapeExtend import ShapeExtend_WireData, ShapeExtend_WireData

output_dict = {}
output_dict_index = 0

data_dic = {"wires_center": [], "cylinders_data": []}

def dumpTopology(shape, level=0):
    """
     Print the details of an object from the top down
    """
    brt = BRep_Tool()
    s = shape.ShapeType()
    if s == TopAbs_VERTEX:
        pnt = brt.Pnt(topods_Vertex(shape))
        # print(".." * level + "<Vertex %i: %s %s %s>" % (hash(shape), pnt.X(), pnt.Y(), pnt.Z()))
    else:
        print(".." * level, end="")
        print(shapeTypeString(shape))
    it = TopoDS_Iterator(shape)
    while it.More():
        shp = it.Value()
        it.Next()
        dumpTopology(shp, level + 1)


def shapeTypeString(shape):
    st = shape.ShapeType()
    s = "?"
    if st == TopAbs_VERTEX:
        s = "Vertex"
    if st == TopAbs_SOLID:
        s = "Solid"
    if st == TopAbs_EDGE:
        s = "Edge"
    if st == TopAbs_FACE:
        s = "Face"
    if st == TopAbs_SHELL:
        s = "Shell"
    if st == TopAbs_WIRE:
        s = "Wire"
    if st == TopAbs_COMPOUND:
        s = "Compound."
    if st == TopAbs_COMPSOLID:
        s = "Compsolid."
    return "%s: %i" % (s, hash(shape))


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

    # aResShape = step_reader.Shape(1)
    # else:
    # print("Error: can't read file.")
    # sys.exit(0)
    # return aResShape


def recognize_face(a_face, output_cylinder_axis_location, toPrint=False):
    """ Takes a TopoDS shape and tries to identify its nature
    whether it is a plane a cylinder a torus etc.
    if a plane, returns the normal
    if a cylinder, returns the radius
    """
    # turn the face (inifinte area) into surface (finite area)
    if toPrint:
        print("##################################")
        print("Face recognition function says: ")
    surf = BRepAdaptor_Surface(a_face, True)
    surf_type = surf.GetType()
    if surf_type == GeomAbs_Plane:

        # look for the properties of the plane
        # first get the related gp_Pln
        gp_pln = surf.Plane()
        location = gp_pln.Location()  # a point of the plane
        normal = gp_pln.Axis().Direction()  # the plane normal
        # then export location and normal to the console output
        if toPrint:
            print("--> plane")
            print("--> Location (global coordinates)", location.X(), location.Y(), location.Z())
            print("--> Normal (global coordinates)", normal.X(), normal.Y(), normal.Z())
    elif surf_type == GeomAbs_Cylinder:

        # look for the properties of the cylinder
        # first get the related gp_Cyl
        gp_cyl = surf.Cylinder()
        location = gp_cyl.Location()  # a point of the axis
        axis = gp_cyl.Axis().Direction()  # the cylinder axis
        # then export location and normal to the console output
        if toPrint:
            print("--> cylinder")
            print("--> Location (global coordinates)", location.X(), location.Y(), location.Z())
            print("--> Axis (global coordinates)", axis.X(), axis.Y(), axis.Z())
            output_cylinder_axis_location.append(
                {
                    "location": [location.X(), location.Y(), location.Z()],
                    "axis": [axis.X(), axis.Y(), axis.Z()]
                })
            return gp_cyl

    elif surf_type == GeomAbs_BSplineSurface:
        if toPrint:
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
        if toPrint:
            print(kind)
        gp_bsrf = surf.Surface()
        # degree = gp_bsrf.NbUKnots()

    elif surf_type == GeomAbs_Torus:
        if toPrint:
            print(a_face)
        kind = "Torus"
        return kind, None, None
    elif surf_type == GeomAbs_BezierSurface:
        if toPrint:
            print(a_face)
        kind = "Bezier"
        return kind, None, None
    elif surf_type == GeomAbs_SurfaceOfRevolution:
        if toPrint:
            print(a_face)
        kind = "Revolution"
        return kind, None, None
    elif surf_type == GeomAbs_SurfaceOfExtrusion:
        if toPrint:
            print(a_face)
        kind = "Extrusion"
        return kind, None, None
    elif surf_type == GeomAbs_OffsetSurface:
        if toPrint:
            print(a_face)
        kind = "Offset"
        return kind, None, None
    elif surf_type == GeomAbs_OtherSurface:
        if toPrint:
            print(a_face)
        kind = "Other"
    else:
        if toPrint:
            print(a_face)
            print(surf_type)
            print("not implemented")
        # TODO there are plenty other type that can be checked
        # see documentation for the BRepAdaptor class
        # https://www.opencascade.com/doc/occt-6.9.1/refman/html/class_b_rep_adaptor___surface.html
    return None


def HoleByWire(aFace):
    # produce wires from face (2)
    d = t.wires_from_face(aFace)
    # enter wires into array (3)
    wires = []
    for dwire in d:
        wires.append(dwire)
    u = k.OuterWire(aFace)
    # find no. of wires(4)
    g = t.number_of_wires_from_face(aFace)
    print(g)
    # if there is more than one wire then there is a hole and the following actions to be taken:
    if g > 1:
        for ijk in range(1, g):
            # add each internal wire in an array of wires (5)
            wires_of_holes.append(wires[ijk])
            # calculate the area size of each internal wire(6)
            recognize_surface(wires[ijk])
            # calculate the center point of each internal wire (7)
            print("Calculating the center of the wire")
            CenterOfWire(wires[ijk])

# This function identify internal wires and their properties
def recognize_surface(new_wire):
    # print(new_wire)
    # new_face=BRepBuilderAPI_MakeFace(wires[i])
    face1 = BRepBuilderAPI_MakeFace(new_wire).Face()
    # print("face1", face1)
    wire_data = ShapeExtend_WireData(new_wire, True, True)
    # wire_data_handle = ShapeExtend_WireData(wire_data)
    # The surface area of the face is calculated and appended to the list
    surface_area = abs(k.TotCross2D(wire_data, face1))
    # print("surface area", surface_area)


# This function calculate the wire surface - -relevant for hole volume calculation
def CenterOfWire(new_wire):
    """
    This function calculate the center of the wire - relevant for hole depth calculation
    """
    global output_dict
    global output_dict_index
    color = Quantity_Color(random.random(),
                           random.random(),
                           random.random(),
                           Quantity_TOC_RGB)
    display.DisplayColoredShape(new_wire, color)
    center = gp_Pnt()
    ordered_edges = t.edges_from_wire(new_wire)
    btool = BRep_Tool()
    for edge in ordered_edges:
        g1 = GProp_GProps()
        # Create a Geom_curve from Edge
        brepgprop_LinearProperties(edge, g1)
        length = g1.Mass()
        # print(length)
        hc = btool.Curve(edge)
        ad = GeomAdaptor_Curve(hc[0])
        discretizer = GCPnts_UniformAbscissa()
        discretizer.Initialize(ad, 10, hc[1], hc[2])
        num_of_points = discretizer.NbPoints()
        # print("number of point", num_of_points)
        for c in range(1, num_of_points):
            p = ad.Value(discretizer.Parameter(c))
            # print(p.X(),p.Y(),p.Z())
            center = gp_Pnt(center.X() + p.X(), center.Y() + p.Y(), center.Z() + p.Z())
    center = gp_Pnt(center.X() / num_of_points, center.Y() / num_of_points, center.Z() / num_of_points)
    # with open('output.json', 'r') as fp_in:
    #     data_dic = json.load(fp_in)
    #     data_dic['wires_center'].append([center.X(), center.Y(), center.Z()])
    #     with open('output.json', 'w') as fp_out:
    #         json.dump(data_dic, fp_out, indent=2)
    #         print("####### JUST WRITTEN a")
    data_dic['wires_center'].append([center.X(), center.Y(), center.Z()])
    print("Center of internal wire", center.X(), center.Y(), center.Z())
    output_dict["center_" + str(output_dict_index)] = [center.X(), center.Y(), center.Z()]
    output_dict_index = output_dict_index + 1
    wire_of_holes_x.append(center.X())
    wire_of_holes_y.append(center.Y())
    wire_of_holes_z.append(center.Z())


def distance(p1, p2):
    x1, y1, z1 = p1.X(), p1.Y(), p1.Z()
    x2, y2, z2 = p2.X(), p2.Y(), p2.Z()
    x = x1 - x2
    y = y1 - y2
    z = z1 - z2
    return np.sqrt(x ** 2 + y ** 2 + z ** 2), [x, y, z]


def EvaluateHole(aFace, toPrint=False):
    """
    This function identify cylindrical holes
    """
    isHole = False
    edges = t.edges_from_face(aFace)
    axisList = []
    boundaryEdges = []

    # 1 check
    # if the boundary edges are closed, the face is a "possible" cylindrical hole (10)
    bIsClosedCurve = False
    isHole = False
    for edge in edges:
        m = g.IsReallyClosed(edge, aFace)
        curve = BRepAdaptor_Curve(edge)
        answer = curve.IsClosed()
        if answer:
            isHole = True
            bIsClosedCurve = True
            break

            # 2 check
    # could be a sectioned cylinder. Collect group of faces contributing to hole  (11)
    if bIsClosedCurve == False and isHole == False:
        surf = BRepAdaptor_Surface(aFace, True)
        surf_type = surf.GetType()
        if surf_type == GeomAbs_Cylinder:
            gp_cyl1 = surf.Cylinder()
            location1 = gp_cyl1.Location()
            # produce for the cylinder face all edges
            edges = t.edges_from_face(aFace)
            for edge in edges:
                # identify all faces which connect to the cylindrical face
                commonfaces = AskFacesFromEdge(edge, aFace)
                for cface in commonfaces:
                    if cface != aFace:  # TODO - implement isEqual for faces
                        surf = BRepAdaptor_Surface(cface, True)
                        surf_type = surf.GetType()
                        # if the attached face to checked face is also a cylinder, check if have the same center
                        if surf_type == GeomAbs_Cylinder:
                            gp_cyl2 = surf.Cylinder()
                            location2 = gp_cyl2.Location()
                            direction2 = gp_cyl2.Axis().Direction()
                            data_dic['cylinders_data'].append({"location": [location1.X(), location1.Y(), location1.Z()],
                                                               "axis": [direction2.X(), direction2.Y(), direction2.Z()]})
                            if location1.Distance(location2) < 1E-3:
                                # with open("output.json", 'r') as fp_in:
                                #     data_dic = json.load(fp_in)
                                #     data_dic['cylinders_locations'].append([location1.X(),location1.Y(), location1.Z() ])
                                #     with open("output.json", 'w') as fp_output:
                                #         json.dump(data_dic, fp_output, indent=2)

                                hole1 = [aFace, cface]
                                print("ADDING A FACE TO A COLORED DISPLAY - IT IS THE HOLE")
                                color = Quantity_Color(random.random(),
                                                       random.random(),
                                                       random.random(),
                                                       Quantity_TOC_RGB)
                                display.DisplayColoredShape(aFace, color)
                                display.DisplayColoredShape(cface, color)
                                HoleFacePair.append(hole1)
                                isHole = True
                                break
                if isHole:
                    break

    # 3 check
    # Identify if the cylinderical shape is a full cylinder or a hole with form of cylinder (12)
    if isHole:
        surf = BRepAdaptor_Surface(aFace, True)
        surf_type = surf.GetType()
        if surf_type == GeomAbs_Cylinder:
            # creation of line between the center of face and the face. line passes through hole center.length of line is holeFaces diagonal
            gp_cyl2 = surf.Cylinder()
            pnt1 = gp_cyl2.Location()
            u_min, u_max, v_min, v_max = breptools_UVBounds(aFace)

            v_center = (v_min + v_max) * 0.5
            u_center = (u_min + u_max) * 0.5

            pnt2 = surf.Value(u_center, v_center)

            dis, direction = distance(pnt1, pnt2)

            my_line = BRepBuilderAPI_MakeEdge(pnt1, pnt2).Edge()
            # Check if line intersection with the shape
            Intersection = Perform_Intersection(shp, my_line)
            if toPrint:
                print("Params ::", u_min, u_max, v_min, v_max)
                print("center params:", v_center, u_center)
                print(pnt1.X(), pnt1.Y(), pnt1.Z(), "\n", "Point on surface", pnt2.X(), pnt2.Y(), pnt2.Z())
                print("The radius is : ", dis)
                print("Intersections inside the suspected hole has ###", Intersection, " intersection points")
            display.DisplayShape(my_line)
            # if line does not intersect, the face is a hole
            if Intersection == 1:
                holeFaces.append(aFace)


def AskFacesFromEdge(iEdge, iFace):
    """
    This function returns all attached faces to an edge of face
    """
    iShape = t
    allFaces = t.faces()
    selFaces = []
    for face in allFaces:
        edges = t.edges_from_face(face)
        for edgex in edges:
            if face != iFace:
                # print(iFace, face,edgex, iEdge)
                if edgex.IsSame(iEdge) == True:
                    # print("edgex", edgex)
                    selFaces.append(face)
    # print(selFaces)
    return selFaces


def Perform_Intersection(shape1, shape2):
    # perform intersection between the two shapes
    section_shp = BRepAlgoAPI_Section(shape1, shape2, False)
    # compute a point on the curve
    section_shp.ComputePCurveOn1(True)
    section_shp.Approximation(True)
    # building the intersection
    section_shp.Build()
    # output of boolean operation
    inter_shp = section_shp.Shape()
    # this explore all the vertex inside the inter_shp
    anVertExplorer = TopExp_Explorer(inter_shp, TopAbs_VERTEX)
    nIntersection = 0
    while anVertExplorer.More():
        nIntersection = nIntersection + 1
        # check for the next vertex
        anVertExplorer.Next()
    display.DisplayShape(inter_shp)
    return nIntersection


# This function checks intersection between two shapes and return the number of intersection' points
if __name__ == '__main__':
    # arrays for internal wires information.

    source_dir = "C:\\Users\\Administrator\\Desktop\\Python-Occ-Algo-Files\\files-dir\\"
    output_cylinder_axis_location = []
    wires_of_holes = []
    wire_of_holes_x = []
    wire_of_holes_y = []
    wire_of_holes_z = []
    # arrays for cylindrical faces information
    holeFaces = []
    HoleFacePair = []
    holeListCopy = []
    finalHoleList = []
    holeShapes = []
    display, start_display, add_menu, add_function_to_menu = init_display()
    display.SetSelectionModeFace()  # switch to Face selection mode
    # first loads the STEP file and display
    # BJB_Terminal_block_3_poles_46.413.1115.50.STEP
    # threads_cubic.STEP
    # 543-03-A001-09_assembled_panel.STEP
    # Philips_driver_9290-008-52100_3Dsolid _2013-04-11.stp
    shp = read_step_file(source_dir + "threads_cubic.STEP")

    # parameters which enable access to relevant function in imported classes
    l = ShapeAnalysis_CheckSmallFace()
    k = shapeanalysis()
    p = BOPTools_AlgoTools()
    g = breptools()
    a = BRep_Tool()
    r = gp_Cylinder()
    curveanalysis = ShapeAnalysis_Curve()
    # turning the download shape into Topo object
    t = Topo(shp)
    # loop over faces only
    for f in t.faces():
        # call the recognition function for identifying surface type
        recognize_face(f, output_cylinder_axis_location, False)
        # call the HoleByWires function for identifying the internal wires which compose the holes(1)
        HoleByWire(f)
        # call the EvaluateHole function for identifying holes (9)
        EvaluateHole(f)

    # Identify number of complete cylinders by comparing list of holeface to thel list of pairs of half cylinders (13)
    holeListCopy = holeFaces.copy()
    nH = len(holeListCopy)
    while nH > 0:
        hi = nH - 1
        h = holeListCopy[hi]
        isPair = False
        for p in HoleFacePair:
            if h.IsEqual(p[0]) or h.IsEqual(p[1]):
                for ih in range(0, len(holeListCopy)):
                    h1 = holeListCopy[ih]
                    if p[0].IsEqual(h1):
                        del holeListCopy[ih]
                        break

                for ih in range(0, len(holeListCopy)):
                    h1 = holeListCopy[ih]
                    if p[1].IsEqual(h1):
                        del holeListCopy[ih]
                        break

                isPair = True
                finalHoleList.append(h)
                break
        if isPair == False:
            for ih in range(0, len(holeListCopy)):
                h1 = holeListCopy[ih]
                if h.IsEqual(h1):
                    finalHoleList.append(h)
                    del holeListCopy[ih]
                    break

        nH = len(holeListCopy)

    wires_centers = []
    for i in range(len(wire_of_holes_x)):
        wires_centers.append([wire_of_holes_x[i], wire_of_holes_y[i], wire_of_holes_z[i]])

    # create full hole from two cylindrical faces
    for i in HoleFacePair:
        holeShape = BRepAlgoAPI_Fuse(i[0], i[1])
        holeShapes.append(holeShape.Shape())
        holeShape1 = holeShape.Shape()
        anVertExplorer = TopExp_Explorer(holeShape1, TopAbs_EDGE)
        EdgeList = []
        while anVertExplorer.More():
            EdgeList.append(anVertExplorer.Current())
            # check for the next edge
            anVertExplorer.Next()
        # print("---->", len(EdgeList))
        BounderyEdge = []
        for e in EdgeList:
            nFaces = t.number_of_edges_from_face(e)
            if nFaces == 1:
                BounderyEdge.append(e)
        # Random color generation
        color = Quantity_Color(random.random(),
                               random.random(),
                               random.random(),
                               Quantity_TOC_RGB)
        # Adding to display the faces of the hole and the boundry edges
        display.DisplayColoredShape(BounderyEdge, color)
        display.DisplayColoredShape(i[0], color)
        display.DisplayColoredShape(i[1], color)

    # calculate number of holes by wire (8)
    m = len(wires_of_holes)
    n = round(m / 2)
    # calculate the depth of the hole
    hezim = m / 2 + 1

    # calculating the depth of the hole - relevant for holes in y directon - possible to develope for other hole directions as well.
    for i in range(0, m):
        copy_wire_of_holes_x = wire_of_holes_x
        copy_wire_of_holes_y = wire_of_holes_y
        copy_wire_of_holes_z = wire_of_holes_z
    count = 0
    for i in range(0, m):
        n = int(i + 1)
        for j in range(n, m):
            if (copy_wire_of_holes_x[i] - copy_wire_of_holes_x[j]) < 0.001 and (
                    copy_wire_of_holes_x[i] - copy_wire_of_holes_x[j]) > (-1):
                if (copy_wire_of_holes_z[i] - copy_wire_of_holes_z[j]) < 0.001 and (
                        copy_wire_of_holes_z[i] - copy_wire_of_holes_z[j]) > (-1):
                    count = count + 1
                    length = copy_wire_of_holes_y[i] - copy_wire_of_holes_y[j]
                    print("hole number ", count)
                    print("depth = ", length)

    print(data_dic)
    with open('output.json', 'w') as fp_out:
        json.dump(data_dic, fp_out, indent=2)
    print(wires_of_holes)
    if input("Display Shape ? (Y - yes , N - shows only colored wires)") == "Y":
        display.DisplayShape(shp, update=True)
    start_display()