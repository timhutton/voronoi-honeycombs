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


def main():
    # setup a VTK scene
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    track = vtk.vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(track)
    ren.SetBackground(0.95, 0.9, 0.85)
    renWin.SetSize(800, 600)

    # make a list of 3D points
    nx, ny, nz = (5, 4, 3)
    pts = list(itertools.product(range(nx), range(ny), range(nz)))
    pts2 = [(x + 0.5, y + 0.5, z + 0.5) for x,y,z in pts] # extend cubic lattice to BCC
    pts.extend(pts2)

    # make a Voronoi structure from them
    v = scipy.spatial.Voronoi(pts)

    # add the finite cells to the scene
    for reg_num in v.point_region:
        indices = v.regions[reg_num]
        if -1 in indices: # external region including point at infinity
            continue
        verts = v.vertices[indices]
        cell = scipy.spatial.ConvexHull(verts)
        faces = cell.simplices
        addSurface(ren, verts, faces, random.random(), random.random(), random.random())


    iren.Start()


if __name__ == "__main__":
    main()