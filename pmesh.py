# pylint: disable=missing-module-docstring

#     ____   ____   ____ _____  __  ___
#    / __ \ / __ \ /  _// ___/ /  |/  /
#   / /_/ // /_/ / / /  \__ \ / /|_/ /
#  / ____// _, _/_/ /  ___/ // /  / /
# /_/    /_/ |_|/___/ /____//_/  /_/

# Part of PRISM CFD project.

# pmesh.py: A converter from meshio supported formats to pmesh.

# MIT License

# Copyright (c) 2021 - M. Emara et al.
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# pylint: disable=too-few-public-methods
# pylint: disable=too-many-instance-attributes
import argparse
from typing import Callable, List

import h5py
import meshio
import numpy as np
import rich


class InvalidMeshException(Exception):
    """Simple Exception wrapper to avoid a catch-all 'except'"""


class UnsupportedElementsException(Exception):
    """Simple Exception wrapper for cases of meshio unsupported elemetes"""


class FaceType:
    """meshio face types"""

    Quad = "quad"
    Triangle = "triangle"


class CellType:
    """meshio cell types"""

    Hex = "hexahedron"
    Tetra = "tetra"
    Wedge = "wedge"


class FromMeshio:
    """A wrapper around meshio.mesh object

    FromMeshio provides access to mesh points, faces, boundary zones
    and owner/neighbor cell connectivity.

    Args:
            mesh_file_path (str): path to mesh file.
    """

    def __init__(self, mesh_file_path: str) -> None:
        try:
            self.mesh_raw = meshio.read(mesh_file_path)
        except meshio.ReadError as exception:
            error = "meshio error, could not open mesh file.\n"
            error += f"{exception}"
            raise InvalidMeshException(error) from exception

        self.points = self.mesh_raw.points
        self.n_points = len(self.points)
        self.n_cells = 0

        # list of points labels of processed faces (all types)
        self.processed_faces = []

        # a list of tuples of ordered labels of processed faces
        # we use `sorted_faces` mainly to check if a face has been processed before or not.
        # to allow for face checking using python's `in`, we store the face points labels as
        # an ordered tuple. check usage example in  face_exists()
        self.sorted_faces = []

        # maps sorted face points labels tuple to the face index in `processed_faces`
        self.face_to_faceid = {}

        # maps face id to a list of cells id. A face is shared by max. 2 cells.
        self.faceid_to_cellid = {}

        # keeps track of the face if to be processed.
        self.current_faceid = 0

        # maps face id to the boundary number defined in meshio.mesh.cell_tags
        self.boundary_faces = {}

        # maps boundary zone id to its list of faces
        self.zoneid_to_faces = {}

        self.check_mesh_elements()

        if not self.has_single_cell_zone():
            raise InvalidMeshException("pmesh do not support multiple zones mesh.")

        for cell_type, faces_fn in (
            (CellType.Hex, self.hex_cell_faces),
            (CellType.Tetra, self.tetra_cell_faces),
            (CellType.Wedge, self.wedge_cell_faces),
        ):
            self.process_cells(cell_type, faces_fn)

        # Allow only 3D meshes, 2D meshes are represented by 3D mesh with one cell in 3rd axis.
        if len(self.sorted_faces) == 0:
            raise InvalidMeshException("Invalid mesh. Found no cells in mesh.")

        self.link_boundary_faces_to_zones_ids()

        self.n_faces = len(self.sorted_faces)

    def check_mesh_elements(self) -> None:
        """Check if mesh contains any unsupported element type

        Raises:
            InvalidMeshException
        """
        supported_elements = set(
            [
                CellType.Hex,
                CellType.Tetra,
                CellType.Wedge,
                FaceType.Quad,
                FaceType.Triangle,
                "line",
            ]
        )
        for cell_block in self.mesh_raw.cells:
            if cell_block.type not in supported_elements:
                raise UnsupportedElementsException(
                    f"pmesh does not support element type '{cell_block.type}'"
                )

    def face_exists(self, face_labels: List) -> bool:
        """Checks if a list of face labels (aka a face) exists or not

        Args:
            face_labels (list): list of face points labels

        Returns:
            bool: True if face exists, False otherwise
        """
        return tuple(sorted(face_labels)) in self.sorted_faces

    def link_face_to_cell(self, face: List, cellid: int) -> None:
        """Associates a face (list of points labels) to an owner or
        a neighbor cell given the cell id.

        Args:
            face (List): list of face points labels
            cellid (int): owner/neighbor cell id
        """
        assert self.face_exists(
            face
        ), "Trying to link face that does not exist to a cell"
        face_id = self.face_to_faceid[tuple(sorted(face))]

        if face_id in self.faceid_to_cellid:
            self.faceid_to_cellid[face_id].append(cellid)
        else:
            self.faceid_to_cellid[face_id] = [cellid]

    def register_face(self, face: List) -> int:
        """Adds a face to list of processed faces and assign an id for it.

        Args:
            face (List): list of face points labels

        Returns:
            int: newly added face id
        """
        sorted_labels = tuple(sorted(face))
        self.face_to_faceid[sorted_labels] = self.current_faceid
        self.sorted_faces.append(sorted_labels)

        # add face points labels to `processed_faces`
        self.processed_faces.append(face)

        self.current_faceid += 1
        return self.current_faceid - 1

    def process_cells(self, cell_type: str, faces_list_fn: Callable) -> None:
        """Given a cell type and function for cell faces coordinates,
        loop over each cell, extract faces
        and construct owner-neighbor connectivity.

        Args:
            cell_type (str): CellType.Hex, CellType.Tetra or CellType.Wedge
            faces_list_fn (Callable): a function that returns a list of faces
                                      (list of points labels)for the type of cell given.
        """
        cell_block = list(
            filter(lambda cell_block: cell_block.type == cell_type, self.mesh_raw.cells)
        )

        if len(cell_block) == 0:
            # mesh has no cells with the given type, nothing to do here
            return

        cells = cell_block[0].data

        for cell_id, cell in enumerate(cells):
            self.n_cells += 1
            faces = faces_list_fn(cell)
            for face in faces:
                # have we met `face` before?
                if not self.face_exists(face):
                    self.register_face(face)

                # link the face to the cell who owns it
                self.link_face_to_cell(face, cell_id)

    def link_boundary_faces_to_zones_ids(self):
        """Associates boundary faces to the id of their boundary zone"""
        boundary_elements_idx = [
            i
            for i in range(len(self.mesh_raw.cells))
            if self.mesh_raw.cells[i].type in (FaceType.Quad, FaceType.Triangle)
        ]

        # All the faces stored in mesh.cells (with types: quad or tri) are boundary faces
        for face_type_idx in boundary_elements_idx:
            for j, face in enumerate(self.mesh_raw.cells[face_type_idx].data):
                faceid = self.face_to_faceid[tuple(sorted(face))]
                boundary_zone_id = self.mesh_raw.cell_data["cell_tags"][face_type_idx][
                    j
                ]
                self.boundary_faces[faceid] = boundary_zone_id

                self.link_face_to_boundary_zone(boundary_zone_id, faceid)

    def link_face_to_boundary_zone(self, boundary_zone_id: int, faceid: int) -> None:
        """A wrapper over zoneid_to_faces map, given a boundary zone id and a face id,
        appends the face id to the list of ids owned by the boundary zone

        Args:
            boundary_zone_id (int): boundary zone id
            faceid (int): face id
        """
        if boundary_zone_id not in self.zoneid_to_faces:
            self.zoneid_to_faces[boundary_zone_id] = [faceid]
        else:
            self.zoneid_to_faces[boundary_zone_id].append(faceid)

    def zoneid_to_name(self, boundary_zone_id: int) -> str:
        """Maps a boundary zone id to boundary zone name

        Args:
            boundary_zone_id (int): boundary zone id

        Returns:
            str: boundary zone name
        """
        if boundary_zone_id == 0:
            return "defaultBoundary"

        tag = self.mesh_raw.cell_tags[boundary_zone_id]
        assert len(tag) == 1, "found a zone id that holds more than one name!"

        return tag[0]

    def is_boundary_face(self, faceid: int) -> bool:
        """Check if faceid is an id of boundary face

        Args:
            faceid (int): face id

        Returns:
            bool: True if boundary face
        """
        return faceid in self.boundary_faces

    def has_single_cell_zone(self) -> bool:
        """Checks if mesh contains more than one zone.

        Returns:
            bool: True if multiple zones found.
        """
        cell_zones = set()
        for i, block in enumerate(self.mesh_raw.cells):
            if block.type in [CellType.Hex, CellType.Tetra, CellType.Wedge]:
                block_zones = np.unique(self.mesh_raw.cell_data["cell_tags"][i])
                for zone in block_zones:
                    cell_zones.add(zone)
        if len(cell_zones) > 1:
            return False
        return True

    def default_faces_count(self) -> int:
        """Returns number of boundary faces with no defined name (default boundary)

        Returns:
            int: number of default faces
        """
        return len(self.zoneid_to_faces[0])

    @staticmethod
    def hex_cell_faces(cell_points: List) -> List[List]:
        """Returns coordinates of 6 faces of a hexahedron cell, using meshio nodes ordering

        Args:
            cell_points (List): list of points defining the cell

        Returns:
            List[List]: list of list of faces points labels
        """
        faces = [
            [cell_points[0], cell_points[3], cell_points[2], cell_points[1]],
            [cell_points[4], cell_points[5], cell_points[6], cell_points[7]],
            [cell_points[0], cell_points[1], cell_points[5], cell_points[4]],
            [cell_points[2], cell_points[3], cell_points[7], cell_points[6]],
            [cell_points[0], cell_points[4], cell_points[7], cell_points[3]],
            [cell_points[1], cell_points[2], cell_points[6], cell_points[5]],
        ]
        return faces

    @staticmethod
    def wedge_cell_faces(cell_points: List) -> List[List]:
        """Returns coordinates of 5 faces of a wedge cell, using meshio nodes ordering

        Args:
            cell_points (List): list of points defining the cell

        Returns:
            List[List]: list of list of faces points labels
        """
        faces = [
            [cell_points[0], cell_points[2], cell_points[1]],
            [cell_points[3], cell_points[4], cell_points[5]],
            [cell_points[3], cell_points[0], cell_points[1], cell_points[4]],
            [cell_points[0], cell_points[3], cell_points[5], cell_points[2]],
            [cell_points[1], cell_points[2], cell_points[5], cell_points[4]],
        ]
        return faces

    @staticmethod
    def tetra_cell_faces(cell_points: List) -> List[List]:
        """Returns coordinates of 4 faces of a tetrahedral cell, using meshio nodes ordering

        Args:
            cell_points (List): list of points defining the cell

        Returns:
            List[List]: list of list of faces points labels
        """
        faces = [
            [cell_points[0], cell_points[2], cell_points[1]],
            [cell_points[1], cell_points[2], cell_points[3]],
            [cell_points[0], cell_points[1], cell_points[3]],
            [cell_points[0], cell_points[3], cell_points[2]],
        ]
        return faces


class PMeshConverter:
    """Converts a meshio.mesh wrapped by FromMeshio to a pmesh file.
    Args:
            mesh_file_path (str): path to input mesh file.
            output_filename (str): name of output pmesh file.
    """

    def __init__(self, input_mesh_path: str, output_filename: str) -> None:
        info(f"loading mesh file: {input_mesh_path}")
        self._in_mesh = FromMeshio(mesh_file_path=input_mesh_path)
        self.output = h5py.File(output_filename, "w")

        self.written_faces_ids = np.zeros((self._in_mesh.n_faces,))
        self.written_faces_count = 0
        self.boundary_zone_data = {}
        self.n_interior_faces = 0

        info("writing points...")
        self.write_points()
        info(f"\t{self._in_mesh.n_points} have been written.")

        info("writing faces...")
        self.write_faces()
        info(f"\t{self.written_faces_count} have been written.")

        info("writing owner/neighbor connectivity relations...")
        self.write_owner()
        self.write_neighbor()

        info("writing boundary faces...")
        self.write_boundary_zones()

        self.output.close()

    def write_points(self):
        """Writes mesh points as an HDF dataset"""
        self.output.create_dataset("points", data=self._in_mesh.points)
        self.output.attrs["n_points"] = self._in_mesh.points.shape[0]

    def write_faces(self) -> None:
        """Writes mesh faces as a variable length HDF dataset"""
        dtype = h5py.vlen_dtype(np.dtype("int32"))
        self.output.create_dataset("faces", (self._in_mesh.n_faces,), dtype=dtype)
        self.output.attrs["n_faces"] = len(self._in_mesh.processed_faces)
        self.write_interior_faces()
        self.write_boundary_faces()

    def write_face(self, faceid: int) -> None:
        """Writes a single face given the face id to faces dataset
        Args:
            faceid (int): face id
        """
        face = self._in_mesh.processed_faces[faceid]
        face_vertex_count = len(face)

        self.output["faces"][self.written_faces_count] = [face_vertex_count, *face]
        self.written_faces_ids[self.written_faces_count] = faceid
        self.written_faces_count += 1

    def write_interior_faces(self) -> None:
        """Writes interior faces as a list of list points to 'faces' dataset"""
        for faceid in self._in_mesh.faceid_to_cellid:
            if self._in_mesh.is_boundary_face(faceid):
                continue
            self.write_face(faceid)
            self.n_interior_faces += 1

    def write_boundary_faces(self) -> None:
        """Writes boundary faces as a list of list points to 'faces' dataset"""
        for zoneid, faces in self._in_mesh.zoneid_to_faces.items():
            zone_name = self._in_mesh.zoneid_to_name(zoneid)

            if zoneid == 0:
                # meshio writes cell_tags for boundary faces with no defined boundary name = 0
                warning = f"found {len(faces)} faces with no defined boundary name."
                warning += ' adding faces to boundary "defaultBoundary"'
                warn(warning)

            start_face_idx = self.written_faces_count

            for faceid in faces:
                self.write_face(faceid)

            end_face_idx = self.written_faces_count - 1

            self.boundary_zone_data[zoneid] = [zone_name, start_face_idx, end_face_idx]

    def write_owner(self) -> None:
        """Writes 'owner' cells dataset, indexed by face id"""
        owner_d = self.output.create_dataset(
            "owner", (len(self._in_mesh.sorted_faces),), dtype="int"
        )
        self.output.attrs['n_cells'] = self._in_mesh.n_cells
        for i, faceid in enumerate(self.written_faces_ids):
            owner_d[i] = self._in_mesh.faceid_to_cellid[faceid][0]

    def write_neighbor(self) -> None:
        """Writes 'neighbor' cells dataset, indexed by face id"""
        neighbor_d = self.output.create_dataset(
            "neighbor", (self.n_interior_faces,), dtype="int"
        )
        self.output.attrs["n_neighbor"] = self.n_interior_faces

        for i in range(self.n_interior_faces):
            faceid = self.written_faces_ids[i]
            neighbor_d[i] = self._in_mesh.faceid_to_cellid[faceid][1]

    def write_boundary_zones(self) -> None:
        """Creates the 'boundary' main group, and sub-group for each boundary face"""
        boundary_g = self.output.create_group("boundary")
        boundary_g.attrs["n_boundary"] = len(self.boundary_zone_data)
        for _, zone_data in self.boundary_zone_data.items():
            zone_name, start_face, end_face = zone_data
            boundary_g.create_group(zone_name)
            boundary_g[zone_name].attrs["start_face"] = start_face
            boundary_g[zone_name].attrs["end_face"] = end_face


def warn(msg: str) -> None:
    """Shows user a warning message"""
    rich.print("[bold yellow]Warning:[/bold yellow]", end="")
    print(msg)


def fatal(msg: str) -> None:
    """Shows user a fatal error message"""
    rich.print("[bold red]Fatal error: [/bold red]", msg)


def info(msg: str) -> None:
    """Ordinary information message"""
    print(msg)


def show_description() -> None:
    """Shows app PRISM logo and app description"""
    desc = r"""
    ____   ____   ____ _____  __  ___
   / __ \ / __ \ /  _// ___/ /  |/  /
  / /_/ // /_/ / / /  \__ \ / /|_/ /
 / ____// _, _/_/ /  ___/ // /  / /
/_/    /_/ |_|/___/ /____//_/  /_/

Part of PRISM CFD project

pmesh.py: A converter from meshio supported formats to pmesh
------------------------------------------------------------
"""
    print(desc)


if __name__ == "__main__":
    show_description()

    parser = argparse.ArgumentParser(
        description="",
    )
    parser.add_argument(
        "input_file",
        metavar="input",
        type=str,
        nargs=1,
        help="Input mesh file name",
    )
    parser.add_argument(
        "output_file",
        metavar="output",
        type=str,
        nargs=1,
        help="Output pmesh file name",
    )

    args = parser.parse_args()
    input_file = args.input_file[0]
    output_file = args.output_file[0]

    try:
        converter = PMeshConverter(input_file, output_file)
    except (InvalidMeshException, UnsupportedElementsException) as e:
        fatal(str(e))
