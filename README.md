# MeColDe
Mesh Collusion Detection for Final Project

This implementation will detect the collision information of 2 triangular meshes in ply format.

To build, simply invoke make command.

To use, pass the names of the 2 meshes you wish to check. It HAS to be TWO meshes of ply format. Optionally, you can pass a 3rd argument, an integer, as an instruction telling the code how many threads should be used for the OMP part. 

The following two command are valid examples:

./MeColDe bun_zipper.ply test_mesh.ply
./MeColDe bun_zipper.ply test_mesh.ply 16

In the repo I already included 2 ply meshes of Stanford bunny should you need them. 
