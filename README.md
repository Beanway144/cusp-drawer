# Cusp Drawer
This is a script which takes in an isosig of a geometric triangulation of a cusped hyperbolic 3-manifold, and outputs a 
picture of the cusp, i.e., a cross-section of the tetrahedra near the torus boundary.

# To use:
Replace `ISO_SIG` to your isosig in draw_cusp.py. Then run
``` sage draw_cusp.py ```

## TODO:
- Pretty up comments + code
- Make it easier to turn on and off vertex/edge coloring
- Fix large edge case failure
- Properly parse input data


## Large Edge Case Failure:
The method used currently fails on triangulations which see the following around an edge: ..., v, w, v, ...\
In other words, some edge class sees the same vertex appear twice with only one other vertex separating them.
This actually seems common among census triangulations, including many figure-8 triangulations.
