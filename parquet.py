"""     Renders a Voronoi diagram from a list of points - https://github.com/timhutton/voronoi-honeycombs
        Copyright (C) 2022 Tim Hutton

        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
import itertools
import random
import scipy
import vtk


def makeVTKScene(width, height, red, green, blue):
    """Creates a VTK window to render into."""
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    track = vtk.vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(track)
    ren.SetBackground(red, green, blue)
    renWin.SetSize(width, height)
    return ren, iren


def makePolyData(verts, faces):
    """Returns a vtkPolyData constructed from the supplied verts and faces."""
    pd = vtk.vtkPolyData()
    pts = vtk.vtkPoints()
    for pt in verts:
        pts.InsertNextPoint( pt[0], pt[1], pt[2] )
    cells = vtk.vtkCellArray()
    for f in faces:
        cells.InsertNextCell( len(f) )
        for v in f:
            cells.InsertCellPoint( v )
    pd.SetPoints(pts)
    pd.SetPolys(cells)
    return pd


def addSurface(renderer, verts, faces, red, green, blue, opacity=1):
    """Add the specified surface to the renderer scene with the specified color."""
    surface = makePolyData(verts, faces)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(surface)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(red, green, blue)
    actor.GetProperty().SetOpacity(opacity)
    #actor.GetProperty().SetRepresentationToWireframe()
    renderer.AddActor(actor)


def addSphere(renderer, pos, radius, red, green, blue):
    """Add a sphere to the renderer scene with the specified color."""
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(radius)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(sphere.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.SetPosition(pos[0], pos[1], pos[2])
    actor.GetProperty().SetColor(red, green, blue)
    renderer.AddActor(actor)


def lerp(a, b, u):
    """Linear interpolation between a (at u=0) and b (at u=1)."""
    return a + (b - a) * u


def getInterpolatedBCCUnitCell(u):
    """Returns a list of points for the specified unit cell, plus its xyz size."""
    # u=0: cubic, u=1: BCC
    p0 = (0, 0, 0) # corner of the cube (doesn't change)
    p1 = (lerp(0, 0.5, u), lerp(0, 0.5, u), lerp(1, 0.5, u)) # 0,0,1 corner or center of cube
    x_size = 1
    y_size = 1
    z_size = lerp(2, 1, u)
    return [p0, p1], (x_size, y_size, z_size)


def main():
    ren, iren = makeVTKScene(800, 600, 0.95, 0.9, 0.85)

    # make a list of 3D points
    #unit_cell = [(0, 0, 0)] # cubic lattice => cubes
    #unit_cell = [(0, 0, 0), (0.5, 0.5, 0.5)] # BCC => truncated octahedra
    #unit_cell = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)] # FCC => rhombic dodecahedra
    #unit_cell = [(0,0,0), (0,0.5,0.5), (0.5,0,0.5), (0.5,0.5,0),
    #             (0.75,0.75,0.75), (0.75,0.25,0.25), (0.25,0.75,0.25), (0.25,0.25,0.75)] # diamond cubic => triakis truncated tetrahedra
    unit_cell = [(0,0,0), (0.5,0.5,0.5), (0,0.25,0.5), (0.25,0.5,0), (0.5,0,0.25),
                 (0,0.75,0.5), (0.75,0.5,0), (0.5,0, 0.75)] # approximation of Weaire–Phelan
    size = (1, 1, 1)
    #unit_cell, size = getInterpolatedBCCUnitCell(1)
    nx, ny, nz = (4, 3, 2)
    genpt = lambda p, po: (p[0] + po[0]*size[0], p[1] + po[1]*size[1], p[2] + po[2]*size[2])
    internal_offsets = list(itertools.product(range(nx), range(ny), range(nz)))
    internal_pts = [genpt(p, po) for p in unit_cell for po in internal_offsets]
    all_offsets = list(itertools.product(range(-1, nx + 1), range(-1, ny + 1), range(-1, nz + 1)))
    external_pts = [genpt(p, po) for p in unit_cell for po in all_offsets if not genpt(p, po) in internal_pts]
    pts = internal_pts + external_pts

    # make a Voronoi structure from them
    print("Finding Voronoi...")
    v = scipy.spatial.Voronoi(pts)

    # add the finite cells to the scene
    print("Rendering...")
    for iVert, reg_num in enumerate(v.point_region[:len(internal_pts)]): # only interested in the internal cells, to avoid boundary effects
        indices = v.regions[reg_num]
        if -1 in indices: # external region including point at infinity
            continue
        verts = v.vertices[indices]
        cell = scipy.spatial.ConvexHull(verts)
        faces = cell.simplices
        addSurface(ren, verts, faces, random.random(), random.random(), random.random(), opacity=1)
        addSphere(ren, internal_pts[iVert], 0.05, random.random(), random.random(), random.random())

    if False:
        # add the unit cube
        cube = vtk.vtkCubeSource()
        cube.SetCenter(0.5, 0.5, 0.5)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(cube.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(0, 0, 0)
        actor.GetProperty().SetRepresentationToWireframe()
        ren.AddActor(actor)

    # render the scene and start the interaction loop
    ren.GetActiveCamera().SetPosition(-7.5, 10.5, 20.5)
    ren.ResetCamera()
    iren.Start()


if __name__ == "__main__":
    main()