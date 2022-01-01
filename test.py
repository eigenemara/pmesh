import pmesh
import pytest
import h5py

def test_meshio_bad_mesh():
    with pytest.raises(pmesh.InvalidMeshException):
        _ = pmesh.FromMeshio('bad_mesh.vtk')

def test_multi_zones():
    with pytest.raises(pmesh.InvalidMeshException):
        _ = pmesh.FromMeshio(r'test_meshes/multi_cell_zones.med')

def test_plain_2d_mesh():
    with pytest.raises(pmesh.InvalidMeshException):
        _ = pmesh.FromMeshio(r'test_meshes/plain_2d_mesh.med')

def test_unsupported_cell_type():
    with pytest.raises(pmesh.UnsupportedElementsException):
        _ = pmesh.FromMeshio(r'test_meshes/unsupported_cell_type.med')

def test_missing_boundary_names():
    mesh = pmesh.FromMeshio(r'test_meshes/square_missing_boundary_names.med')
    assert mesh.default_faces_count() == 8

def test_boundary_names():
    mesh = pmesh.FromMeshio('test_meshes/square_complete_boundary_name.med')
    names = set(['front', 'bottom', 'top', 'right', 'left', 'back'])
    
    zones_ids = mesh.zoneid_to_faces.keys()
    mesh_boundary_names = set([mesh.zoneid_to_name(i) for i in zones_ids])
    
    assert names == mesh_boundary_names

def test_writing():
    _ = pmesh.PMeshConverter(
        'test_meshes/square_complete_boundary_name.med',
        'writing_test.pmesh'
    )
    
    file = h5py.File('writing_test.pmesh', 'r')
    assert file.attrs['n_faces'] == 20
    assert file.attrs['n_points'] == 18

def test_tetra_wedge():
    _ = pmesh.PMeshConverter(
        'test_meshes/tetra_wedge.med',
        'writing_test.pmesh'
    )
    
    file = h5py.File('writing_test.pmesh', 'r')
    boundary_faces_count = 0
    
    for boundary in file['boundary']:
        start_face_idx = file[f'boundary/{boundary}'].attrs['start_face']
        end_face_idx = file[f'boundary/{boundary}'].attrs['end_face']
        boundary_faces_count += end_face_idx - start_face_idx + 1
    
    assert boundary_faces_count == 272
    assert file.attrs['n_points'] == 213
