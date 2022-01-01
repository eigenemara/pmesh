# pmesh
*Part of PRISM project.*

[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/EigenEmara/pmesh.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/EigenEmara/pmesh/context:python)
[![build](https://github.com/EigenEmara/pmesh/actions/workflows/python-app.yml/badge.svg)](https://github.com/EigenEmara/pmesh/actions/workflows/python-app.yml)
[![Pylint](https://github.com/EigenEmara/pmesh/actions/workflows/pylint.yml/badge.svg)](https://github.com/EigenEmara/pmesh/actions/workflows/pylint.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

**pmesh** is a compact, single zone, finite volume mesh format. pmesh stores points, faces, cell connectivity (owner-neighbor as inspired by OpenFOAM's polyMesh), boundary faces and arbitrary scalar/vector fields defined at cell center. Current implementation, as found in `pmesh.py` builds up over HDF5 file format.

Currently, pmesh only supports hexahedral, tetrahedral and wedge cells only, pyramid type support is planned to be added soon.

## pmesh outline
1. points: 
a list of 3D cartesian coordinates. Size of points is stored in `n_points` attribute.

2. faces: 
a list of list of points labels, with the first item in each row is the number of points/vertices that defines the face. Size of faces is stored in `n_faces` attribute.

For example: an entry in the form `3, 1, 3, 2` indicates a triangular face (first value in the list = 3), remaining labels are the indices of the points that defines such face.

3. owner: 
a list of of labels of cell that owns which faces.

For example: at index 23 of `owner` we might see the value `13`, this indicates that cell `13` owns face `23`. Size of owner cells is stored in `n_cells` attribute.

4. neighbor:
a list of labels that defines neighbor cells to interior faces.
For example: at index 23 of `neighbor` we might see the value `14`, this indicates that cell `14` is also sharing face `23`. Size of neighbor cells is stored in `n_neighbor` attribute.


Example of the outline of a pmesh file:

        pmesh
        ├── n_points (attribute: int)
        ├── n_faces (attribute: int)
        ├── n_cells (attribute: int)
        ├── n_neighbor(attribute: int)
        ├── points
        │   ├── 0.1, 0.0, 0.0
        │   ├── 0.1, 0.1, 0.0
        │   ├── ...
        ├── faces
        │   ├── 3,1,3,2
        │   ├── 3,2,4,1
        │   ├── ...
        ├── owner
        │   ├── 0
        │   ├── 2
        │   └── ...
        ├── neighbor
        │   ├── 1
        │   ├── 1
        │   └── ...
        ├── boundary
        │   ├── inlet
        │   │   ├── start_face (attribute)
        │   │   └── end_face (attribute)
        │   ├── outlet
        │   │   ├── start_face (attribute)
        │   │   └── end_face (attribute)
        │   └── ...
        └── time
            ├── 0
            │   ├── v
            │   │   ├── 2.5, 1.0, 0.0
            │   │   └── ...
            │   └── p
            │       ├── 1.103
            │       ├── 1.133
            │       └── ...
            └── 1
                └── ...

## `pmesh.py`
pmesh.py should handle all supported mesh formats by `meshio`, however, the tool was mostly tested on SALOME's MED format, and MED was the main target of conversion during development.

### Usage
        pmesh.py [input] [output]

### requirements
        meshio, h5py, numpy, rich
