"""     Making space-filling polyhedra using Voronoi cells - https://github.com/timhutton/voronoi-honeycombs
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

import itertools
import math
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
    """Add the specified surface to the renderer scene."""
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
    return actor


def addSphere(renderer, pos, radius, color):
    """Add a sphere to the renderer scene."""
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(radius)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(sphere.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.SetPosition(*pos)
    actor.GetProperty().SetColor(*color)
    renderer.AddActor(actor)


def addUnitCube(renderer, size):
    """Add the unit cube to the renderer scene."""
    cube = vtk.vtkCubeSource()
    cube.SetBounds(0, size[0], 0, size[1], 0, size[2])
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(cube.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0, 0, 0)
    actor.GetProperty().SetRepresentationToWireframe()
    renderer.AddActor(actor)


def addLabel(pos, text, renderer, color, scale):
    """Add a text label to the renderer scene."""
    vector_text = vtk.vtkVectorText()
    vector_text.SetText(text);
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(vector_text.GetOutputPort())
    actor = vtk.vtkFollower()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*color)
    actor.SetScale(scale)
    actor.SetPosition(*pos)
    renderer.AddActor(actor)


def lerp(a_vals, b_vals, u):
    """Linear interpolation between a_vals at u=0 and b_vals at u=1."""
    return [a + (b - a) * u for a,b in zip(a_vals, b_vals)]


def lerpUnitCell(a, b, u):
    """Interpolates two unit cells, returning 3D points and dimensions of the resulting unit cell."""
    a_pts = [[x / a["scale"] for x in pt] for pt in a["unit_cell"]]
    b_pts = [[x / b["scale"] for x in pt] for pt in b["unit_cell"]]
    return [lerp(pa, pb, u) for pa,pb in zip(a_pts, b_pts)], lerp(a["size"], b["size"], u)


def temporallyConsistentRandomColor(i, colors):
    """Return a random RGB [0,1] color, and remember previous answers."""
    if not i in colors:
        colors[i] = random.random(), random.random(), random.random()
    return colors[i]


def getUnitCell(label):
    """ Returns one of a family of unit cells each with 8 cell centers that we can interpolate between.
           unit_cell : List of 3D points expressed in some convenient scale
           scale     : Factor to divide by to obtain the points we want
           size      : Dimensions of the resulting unit cell after dividing points by scale
    """
    unit_cells = {
        "cubic":   { "unit_cell": [(0,0,0),(0,1,0),(1,1,0),(1,0,0),(0,0,1),(0,1,1),(1,1,1),(1,0,1)], "scale": 1, "size": [2,2,2], "name": "Cubic" },
        "bcc":     { "unit_cell": [(0,0,0),(0,2,0),(2,2,0),(2,0,0),(1,1,1),(1,3,1),(3,3,1),(3,1,1)], "scale": 2, "size": [2,2,1], "name": "BCC" },
        "fcc":     { "unit_cell": [(0,0,0),(0,1,1),(1,1,0),(1,0,1),(0,0,2),(0,1,3),(1,1,2),(1,0,3)], "scale": 1, "size": [2,2,4], "name": "FCC" },
        "diamond": { "unit_cell": [(0,0,0),(1,3,1),(2,2,0),(3,1,1),(1,1,3),(0,2,2),(3,3,3),(2,0,2)], "scale": 2, "size": [2,2,2], "name": "Diamond cubic" },
        "a15":     { "unit_cell": [(0,0,0),(1,2,0),(3,2,0),(2,0,1),(0,1,2),(0,3,2),(2,2,2),(2,0,3)], "scale": 2, "size": [2,2,2], "name": "A15 crystal (Weaire-Phelan)" },
        "laves":   { "unit_cell": [(0,0,0),(1,3,0),(2,3,1),(3,0,1),(0,1,3),(1,2,3),(2,2,2),(3,1,2)], "scale": 2, "size": [2,2,2], "name": "Laves graph" } }
    return unit_cells[label]


def animateTransitions(a="cubic", b="laves", showUnitCell=False):
    """Interactive display of animated transition between two honeycombs."""

    window_size = 800, 600
    background_color = 0.95, 0.9, 0.85
    ren, renWin, iren = makeVTKWindow(window_size, background_color, "Voronoi Honeycombs")
    colors = {}

    # Define the animation sequence
    num_frames = 30
    u_values = [iFrame / num_frames for iFrame in range(num_frames + 1)]
    u_values.extend(u_values[::-1])
    pauses = [num_frames]

    for iFrame, u in enumerate(u_values):

        # Make a list of 3D points
        unit_cell, size = lerpUnitCell(getUnitCell(a), getUnitCell(b), u)
        nx, ny, nz = (1, 1, 1) if showUnitCell else (2, 2, 2)
        genpt = lambda p, offset: [p[i] + offset[i]*size[i] for i in range(3)]
        internal_offsets = list(itertools.product(range(nx), range(ny), range(nz)))
        internal_pts = [genpt(p, offset) for p in unit_cell for offset in internal_offsets]
        all_offsets = list(itertools.product(range(-1, nx + 1), range(-1, ny + 1), range(-1, nz + 1)))
        external_pts = [genpt(p, offset) for p in unit_cell for offset in all_offsets if not genpt(p, offset) in internal_pts]

        # Compute a Voronoi structure from the list of 3D points
        v = scipy.spatial.Voronoi(internal_pts + external_pts)

        # Add the Voronoi cells to the scene
        for iVert, reg_num in enumerate(v.point_region[:len(internal_pts)]): # only interested in the internal cells, to avoid boundary effects
            indices = v.regions[reg_num]
            if -1 in indices: # external region including point at infinity
                continue
            verts = v.vertices[indices]
            faces = scipy.spatial.ConvexHull(verts).simplices
            addSurface(ren, verts, faces, temporallyConsistentRandomColor(iVert, colors), opacity=1, wireframe=showUnitCell)
            if showUnitCell:
                addSphere(ren, internal_pts[iVert], 0.05, temporallyConsistentRandomColor(iVert, colors))
                addUnitCube(ren, size)

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


def addParquetSlice(unit_cell_pts, size, pos, nx, ny, internal_pts, external_pts, coords=None):
    """Add points to internal_pts/external_pts for each unit cell in the X-Y range specified."""
    for ix,iy in itertools.product(range(-1, nx + 1), range(-1, ny + 1)):
        for unit_cell_pt in unit_cell_pts:
            offset_unit_cell_pt = (pos[0] + size[0]*ix + unit_cell_pt[0],
                                   pos[1] + size[1]*iy + unit_cell_pt[1],
                                   pos[2] + unit_cell_pt[2])
            if ix >= 0 and ix < nx and iy >= 0 and iy < ny:
                if not coords is None:
                    coords[len(internal_pts)] = ix, iy, pos[2]
                internal_pts.append(offset_unit_cell_pt)
            else:
                external_pts.append(offset_unit_cell_pt)


def animateParquetDeformation():
    """Interactive display of a parquet deformation between different honeycombs."""

    window_size = 720, 480
    background_color = 0.95, 0.9, 0.85
    ren, renWin, iren = makeVTKWindow(window_size, background_color, "Voronoi Honeycombs")
    label_color = 0, 0, 0
    label_scale = 0.5
    colors = {}

    ren.GetActiveCamera().SetPosition(-2, -3, -7)
    ren.GetActiveCamera().SetFocalPoint(0, 0, 0)
    ren.GetActiveCamera().SetViewUp(0, -1, 0)
    light = vtk.vtkLight()
    light.SetPosition(-5, -7, -5)
    light.SetAmbientColor(0.4, 0.4, 0.4)
    light.SetShadowAttenuation(0.3)
    ren.AddLight(light)
    first = True

    nx, ny, nz = 2, 1, 10
    sequence = ["cubic", "cubic", "laves", "laves", "diamond", "diamond", "a15", "a15", "cubic", "cubic", "cubic"]

    renWinToImage = vtk.vtkWindowToImageFilter()
    renWinToImage.SetInput(renWin)
    pngWriter = vtk.vtkPNGWriter()
    pngWriter.SetInputConnection(renWinToImage.GetOutputPort())

    # Collect points by stacking unit cells along the z-axis
    num_frames = 200
    for iFrame in range(0, num_frames):
        u_center = len(sequence[:-2]) * iFrame / num_frames
        internal_pts = [] # ones we will display the Voronoi cell for
        external_pts = [] # ones with possible boundary effects that we need to compute Voronoi but don't display
        coords = {}
        # iterate from the center slice in both directions
        for d in [1, -1]:
            pos = [0, 0, 0]
            for iz in range(0, nz):
                u = u_center + d * iz / nz
                u_pre = math.floor(u)
                u_post = u_pre + 1
                u_frac = u - u_pre
                #print(iFrame, u_center, d, iz, u, u_pre, u_post, u_frac, len(sequence))
                cell1 = getUnitCell(sequence[u_pre])
                cell2 = getUnitCell(sequence[u_post])
                unit_cell_pts, size = lerpUnitCell(cell1, cell2, u_frac)
                if iz > 0 or d == 1:
                    if iz < nz - 1:
                        addParquetSlice(unit_cell_pts, size, pos, nx, ny, internal_pts, external_pts, coords)
                    else:
                        # add hidden slice to avoid boundary effects
                        addParquetSlice(unit_cell_pts, size, pos, nx, ny, external_pts, external_pts, coords)
                pos[2] += d * size[2]
        for iVert, internal_pt in enumerate(internal_pts):
            addSphere(ren, internal_pt, 0.05, temporallyConsistentRandomColor(iVert, colors))

        # Run Voronoi over all the points
        v = scipy.spatial.Voronoi(internal_pts + external_pts)

        # Add the cells to the scene
        for iVert in range(len(internal_pts)): # only interested in the internal cells, to avoid boundary effects
            indices = v.regions[v.point_region[iVert]]
            if -1 in indices: # external region including point at infinity
                continue
            verts = v.vertices[indices]
            faces = scipy.spatial.ConvexHull(verts).simplices
            pt_coords = coords[iVert]
            wireframe = False # iVert > 0 and pt_coords[0] < 1 and math.fabs(pt_coords[2]) < 3
            addSurface(ren, verts, faces, temporallyConsistentRandomColor(iVert, colors), opacity=1, wireframe=wireframe)

        # Render
        renWin.Render()
        if first:
            print("Controls:")
            print("  mouse left drag: rotate the scene")
            print("  mouse right drag up/down, or mouse wheel: zoom in/out")
            print("  shift + mouse left drag: pan")
            print("\nPress 'q' to animate")
            first = False
            iren.Start()
            print("Animating...")
        renWin.Render()
        pngWriter.SetFileName(f"render_{iFrame:05d}.png")
        renWinToImage.Modified()
        pngWriter.Write()
        ren.RemoveAllViewProps()


if __name__ == "__main__":

    #animateTransitions(a="cubic", b="laves", showUnitCell=False)
    animateParquetDeformation()
