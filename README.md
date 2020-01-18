# Crease Aligned Cross Fields

## Introduction
This code includes algorithms for computing crease-aligned cross fields on triangle meshes:

## External Dependencies
- [Arff](https://github.com/dpa1mer/arff) 
- [Mosek](https://www.mosek.com) 9.0 ([C++ Fusion API](https://docs.mosek.com/9.0/cxxfusion/index.html#))
- Intel [TBB](https://github.com/intel/tbb)
- [gptoolbox](https://github.com/alecjacobson/gptoolbox)

## Installation
First, remember to build and install the Mosek Fusion API as described
[here](https://docs.mosek.com/9.0/cxxfusion/install-interface.html).

## Usage
The main command for computing fields is `SolveLpCrossField`. See `example` for building and run instructions.

### Loading Models
Some triangle meshes in `obj` format are included in the `Meshes` directory for convenience. See `example` for how to load.

### Computing Cross Fields
The following command computes a cross field on the triangle mesh whose vertices are X and triangles are T.
```crossField = SolveLpCrossField(X, T, '', 0, 2, True)
```

### Visualization
The sixth argument is a flag for toggling visualization. Leave as true to get a figure of the cross field. Tools are available online for computing streamlines of the cross field or finding singular locations. 
(https://github.com/avaxman/Directional/blob/master/docs/tutorial.md)
