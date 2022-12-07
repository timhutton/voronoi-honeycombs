"""Renders a Voronoi diagram from a list of points."""
import itertools
import random
import scipy
import vtk


def makePolyData( verts, faces ):
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


def addSurface(renderer, verts, faces, red, green, blue):
    """Add the specified surface to the renderer scene with the specified color."""
    surface = makePolyData(verts, faces)
    surfaceMapper = vtk.vtkPolyDataMapper()
    surfaceMapper.SetInputData(surface)
    surfaceActor = vtk.vtkActor()
    surfaceActor.SetMapper(surfaceMapper)
    surfaceActor.GetProperty().SetColor(red, green, blue)
    renderer.AddActor(surfaceActor)


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


def main():
    ren, iren = makeVTKScene(800, 600, 0.95, 0.9, 0.85)

    # make a list of 3D points
    #unit_cell = [(0, 0, 0)] # cubic lattice
    unit_cell = [(0, 0, 0), (0.5, 0.5, 0.5)] # BCC
    #unit_cell = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)] # FCC
    nx, ny, nz = (5, 4, 3)
    internal_offsets = list(itertools.product(range(nx), range(ny), range(nz)))
    internal_pts = [(x + ox, y + oy, z + oz) for x,y,z in unit_cell for ox,oy,oz in internal_offsets]
    all_offsets = list(itertools.product(range(-1, nx + 1), range(-1, ny + 1), range(-1, nz + 1)))
    external_pts = [(x + ox, y + oy, z + oz) for x,y,z in unit_cell for ox,oy,oz in all_offsets if not (x + ox, y + oy, z + oz) in internal_pts]
    pts = internal_pts + external_pts

    # make a Voronoi structure from them
    print("Finding Voronoi...")
    v = scipy.spatial.Voronoi(pts)

    # add the finite cells to the scene
    print("Rendering...")
    for reg_num in v.point_region[:len(internal_pts)]: # only interested in the internal cells, to avoid boundary effects
        indices = v.regions[reg_num]
        if -1 in indices: # external region including point at infinity
            continue
        verts = v.vertices[indices]
        cell = scipy.spatial.ConvexHull(verts)
        faces = cell.simplices
        addSurface(ren, verts, faces, random.random(), random.random(), random.random())

    # render the scene and start the interaction loop
    ren.GetActiveCamera().SetPosition(0.5, 10.5, 20.5)
    ren.ResetCamera()
    iren.Start()


if __name__ == "__main__":
    main()