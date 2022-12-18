## How do I use it? ##

Run using Python 3:
```
python -m pip install scipy vtk
python parquet.py
```

## What does it do? ##

The script renders space-filling polyhedra using Voronoi cells. To animate the transitions between different patterns we move the points and run Voronoi on every frame. See inside the script for more options.

Discussion thread: https://mathstodon.xyz/@timhutton/109469821382910347

### Examples: ###

A tour of six different space-filling polyhedra: the cubic lattice, the A15 crystal (Weaire-Phelan), body-centered cubic (truncated octahedra), the Laves graph (triamond), face-centered cubic (rhombic dodecahedra), and diamond cubic (triakis truncated tetrahedra). In the lower-left is shown the unit cell, with 8 cell centers and two of the polyhedra shown in wireframe. On each frame we move the cell centers and run the Voronoi algorithm to obtain the polyhedra.

https://user-images.githubusercontent.com/647092/208315874-e75d5381-f20b-4d7d-96fa-6c587a575b78.mp4

The A15 crystal, which is the Voronoi approximation of the [Weaire-Phelan structure](https://en.wikipedia.org/wiki/Weaire%E2%80%93Phelan_structure):
![image](https://user-images.githubusercontent.com/647092/206327613-dec7f406-567a-4541-8f94-38f37e842a0c.png)

The [diamond-cubic](https://en.wikipedia.org/wiki/Diamond_cubic) structure, showing the unit cell on the left and more of the cells on the right:

https://user-images.githubusercontent.com/647092/206363589-336ff4e7-72a4-4b69-b06b-1f90c6ac1629.mp4
