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

print("Loading libraries...")

import functools
import itertools
import random

import scipy
import vtk


def makeVTKWindow(window_size, background_color, title):
    """Creates a VTK window to render into."""
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    renWin.SetWindowName(title)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    track = vtk.vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(track)
    ren.SetBackground(*background_color)
    renWin.SetSize(*window_size)
    return ren, renWin, iren


def makePolyData(verts, faces):
    """Returns a vtkPolyData constructed from the supplied verts and faces."""
    pd = vtk.vtkPolyData()
    pts = vtk.vtkPoints()
    for pt in verts:
        pts.InsertNextPoint(*pt)
    cells = vtk.vtkCellArray()
    for f in faces:
        cells.InsertNextCell(len(f))
        for v in f:
            cells.InsertCellPoint(v)
    pd.SetPoints(pts)
    pd.SetPolys(cells)
    return pd


def addSurface(renderer, verts, faces, color, opacity=1, wireframe=False):
    """Add the specified surface to the renderer scene with the specified color."""
    surface = makePolyData(verts, faces)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(surface)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*color)
    actor.GetProperty().SetOpacity(opacity)
    if wireframe:
        actor.GetProperty().SetRepresentationToWireframe()
    renderer.AddActor(actor)


def addSphere(renderer, pos, radius, color):
    """Add a sphere to the renderer scene with the specified color."""
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(radius)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(sphere.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.SetPosition(*pos)
    actor.GetProperty().SetColor(*color)
    renderer.AddActor(actor)


def lerp(a_vals, b_vals, u):
    """Linear interpolation between a_vals at u=0 and b_vals at u=1."""
    return [a + (b - a) * u for a,b in zip(a_vals, b_vals)]


def lerpUnitCell(a, b, u):
    """Interpolates two unit cells, returning 3D points and dimensions of the resulting unit cell."""
    a_pts = [[x / a["scale"] for x in pt] for pt in a["unit_cell"]]
    b_pts = [[x / b["scale"] for x in pt] for pt in b["unit_cell"]]
    return [lerp(pa, pb, u) for pa,pb in zip(a_pts, b_pts)], lerp(a["size"], b["size"], u)


@functools.lru_cache
def temporallyConsistentRandomColor(i):
    """Return a random RGB [0,1] color, and remember previously answers."""
    return random.random(), random.random(), random.random()


#TODO: code for parquet deformation

def main():

    window_size = 800, 600
    background_color = 0.95, 0.9, 0.85
    ren, renWin, iren = makeVTKWindow(window_size, background_color, "Voronoi Honeycombs")

    # Define the animation sequence
    num_frames = 30
    u_values = [iFrame / num_frames for iFrame in range(num_frames + 1)]
    u_values.extend(u_values[::-1])
    pauses = [num_frames]

    # A family of unit cells each with 8 cell centers that we can interpolate between
    #   unit_cell : List of 3D points expressed in some convenient scale.
    #   scale     : Factor to divide by to obtain the points we want.
    #   size      : Dimensions of the resulting unit cell after dividing points by scale.
    cubic2x2x2   = { "unit_cell": [(0,0,0),(0,1,0),(1,1,0),(1,0,0),(0,0,1),(0,1,1),(1,1,1),(1,0,1)], "scale": 1, "size": [2,2,2] }
    bcc2x2x1     = { "unit_cell": [(0,0,0),(0,2,0),(2,2,0),(2,0,0),(1,1,1),(1,3,1),(3,3,1),(3,1,1)], "scale": 2, "size": [2,2,1] }
    fcc1x1x2     = { "unit_cell": [(0,0,0),(0,1,1),(1,1,0),(1,0,1),(0,0,2),(0,1,3),(1,1,2),(1,0,3)], "scale": 1, "size": [2,2,4] }
    diamondCubic = { "unit_cell": [(0,0,0),(1,3,1),(2,2,0),(3,1,1),(1,1,3),(0,2,2),(3,3,3),(2,0,2)], "scale": 2, "size": [2,2,2] }
    weairePhelan = { "unit_cell": [(0,0,0),(1,2,0),(3,2,0),(2,0,1),(0,1,2),(0,3,2),(2,2,2),(2,0,3)], "scale": 2, "size": [2,2,2] }

    for iFrame, u in enumerate(u_values):

        # Make a list of 3D points
        unit_cell, size = lerpUnitCell(cubic2x2x2, diamondCubic, u)
        nx, ny, nz = (2, 2, 2)
        genpt = lambda p, offset: [p[i] + offset[i]*size[i] for i in range(3)]
        internal_offsets = list(itertools.product(range(nx), range(ny), range(nz)))
        internal_pts = [genpt(p, offset) for p in unit_cell for offset in internal_offsets]
        all_offsets = list(itertools.product(range(-1, nx + 1), range(-1, ny + 1), range(-1, nz + 1)))
        external_pts = [genpt(p, offset) for p in unit_cell for offset in all_offsets if not genpt(p, offset) in internal_pts]
        pts = internal_pts + external_pts

        # Compute a Voronoi structure from the list of 3D points
        v = scipy.spatial.Voronoi(pts)

        # Add the Voronoi cells to the scene
        for iVert, reg_num in enumerate(v.point_region[:len(internal_pts)]): # only interested in the internal cells, to avoid boundary effects
            indices = v.regions[reg_num]
            if -1 in indices: # external region including point at infinity
                continue
            verts = v.vertices[indices]
            faces = scipy.spatial.ConvexHull(verts).simplices
            addSurface(ren, verts, faces, temporallyConsistentRandomColor(iVert), opacity=1, wireframe=False)
            #addSphere(ren, internal_pts[iVert], 0.05, temporallyConsistentRandomColor(iVert))

        if iFrame == 0:
            print("Controls:")
            print("  mouse left drag: rotate the scene")
            print("  mouse right drag up/down, or mouse wheel: zoom in/out")
            print("  shift + mouse left drag: pan")
            print("\nPress 'q' to start the animation")
            ren.GetActiveCamera().SetPosition(-7.5, 5.5, -20.5)
            ren.ResetCamera()
            iren.Start()
            ren.RemoveAllViewProps()
            print("Animating...")
        elif iFrame in pauses:
            print("Press 'q' to continue")
            renWin.Render()
            iren.Start()
            ren.RemoveAllViewProps()
            print("Animating...")
        elif iFrame == len(u_values) - 1:
            renWin.Render()
            print("Press 'q' to exit")
            iren.Start()
        else:
            renWin.Render()
            ren.RemoveAllViewProps()


if __name__ == "__main__":
    main()
