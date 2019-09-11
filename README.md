# A h-adaptive scaled boundary finite element method based on maximisation of the error decrease rate

## Full report
https://github.com/congthanh184/sbfem-quad-lagrange/blob/master/ENS6127%20JO%20171%20Final%20Report%20by%20BUI%20Cong%20Thanh.pdf

## Contents

- Folder `quad_lagrange/` : Contains script of SBFEM using lagrange shape function

- Folder `projection-vectorize-h-adaptive/` : Contains script to solve h-adaptive method for square hole problem

- Folder `projection-vectorize-h-adaptive -square-plate/` : Contains script to solve h-adaptive method for square plate under vertical load

## Problems

- Square hole

![square_hole](https://github.com/congthanh184/sbfem-quad-lagrange/blob/master/Docs/example_square_hole.PNG)

- Square plate

![square_plate](https://github.com/congthanh184/sbfem-quad-lagrange/blob/master/Docs/square_plate_edited.png)

## Results

- Quad-Lagrange smooth recovery

![smooth_strain](https://github.com/congthanh184/sbfem-quad-lagrange/blob/master/Docs/compare_raw_smooth_strain.PNG)

![smooth_stress](https://github.com/congthanh184/sbfem-quad-lagrange/blob/master/Docs/compare_raw_smooth_stress.PNG)

- Method Errors

![patch_project](https://github.com/congthanh184/sbfem-quad-lagrange/blob/master/Docs/patch_vs_project_modified_error.png)

![patch_vs_projection](https://github.com/congthanh184/sbfem-quad-lagrange/blob/master/Docs/stress_patch_vs_projection.png)

- Refined mesh

![strain_mesh_stress](https://github.com/congthanh184/sbfem-quad-lagrange/blob/master/Docs/strain_mesh_stress_5.PNG)
