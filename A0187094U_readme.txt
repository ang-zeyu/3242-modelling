Additional notes:
- The red edge indicating boundary edges may disappear when selecting the triangle(s)
  connecting those edges due to z-fighting as the selected triangles have dark green edge highlights.
- A0187094U_mug.3ds included for testing reading .3ds files
- Usage instructions are in cli

Custom mesh files for testing final boss:

Mesh decimation:
- A0187094U_linkCheckShouldFail.obj replicates the shark-fin scenario in lecture, can be used to test link check
- Boundary check can be tested with the cow mesh: ctrl-alt select the top part separating the two halves of the cow;
  the boundary edges will always stay the same.
- Flipped triangle check can be tested with any mesh by repeatedly decimating.

Mesh subdivision, Mesh relaxation, Laplacian smoothing:
- any mesh would do

Laplacian Deformation:
- any mesh would do, but best tested on meshes with fewer CCs to see the results of shape preservation post deformation,
  because selected triangles of non-connected components that are connected to the selected vertex will stay the same
  in preserving shape by the optimization objective.
