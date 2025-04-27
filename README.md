# ooc-meshing
Out-Of-Core Tetrahedal Interface Remesher


# Requirements:
##### Install external Dependencies
```
apt-get install libtiff-dev
```

##### Bundled dependencies, downloaded and built during cmake
- geogram: [https://github.com/BrunoLevy/geogram](https://github.com/BrunoLevy/geogram)
- cli11: [https://github.com/CLIUtils/CLI11](https://github.com/CLIUtils/CLI11)

# Building
```
cd build
cmake ..
make
```

# Useful Software
To view .geogram meshes, install vorpaview from the geogram repository:
```
mkdir vorpaview
cd vorpaview
git clone https://github.com/BrunoLevy/geogram.git .
git submodule init
mkdir build
cd build
cmake ..
make [install] vorpaview
```


### Options
I still have to figure out how to do the entire CMAKE-Option-Configuration-Thingy...
- ```UNROLL_LOOPS```: Unrolls most loops in geometric predicates
- ```PARALLEL_LOCAL_OPERATIONS```: Parallelizes most element-local operations. For more information see [geogram's process library](https://github.com/BrunoLevy/geogram/blob/main/src/lib/geogram/basic/process.h)

Planned:
- Option to enable/disable timers
- Option to enable/disable building specific subprograms

# Usage
The project is split into multiple sub-programs, which can be incorporated in any workflow-framework.

### 1. Surface-Extraction
CGAL Surface-Extraction wrapper, that takes a stack of TIFF-Data and generates a surface mesh for use in meshers like TetWild. Also generates a sizing-field mesh if desired.

<b>Implementation pending.</b>

### 1.1. Tetrahedral Mesh Generation (maybe)
Will maybe implement this, probably just a wrapper around fTetWild.

<b>Implementation pending.</b>

### 2. Interface-Extraction
Takes two arbitrary tetrahedral meshes and a definition of an axis aligned plane, and generates a mesh containing the triangulated
union of vertices from both meshes that lie on the given plane.

Example usage:
```
./build/src/interface-generator/OutOfCore_Meshing_InterfaceGenerator
    -a "test/cylinder.mesh"           # first mesh
    -b "test/cube.mesh"               # second mesh
    -o "test/generated_interface.obj" # output
    -p "-1.0"                         # extent of the plane along the axis
    -x "Z"                            # axis the plane aligns to
    -e "0.00001"                      # epsilon around the extent of the axis
```
### 3. Interface-Insertion
Takes an arbitrary tetrahedral mesh and an interface-triangulation, and inserts the triangulation into the mesh via element-local operations. Outputs the new mesh, as well as metadata about the inserted interface for merging.

Example usage: (<b>Pending after refactoring.</b>)
```
./build/src/interface-inserter/OutOfCore_Meshing_InterfaceInserter
    -m "test/cube.mesh"               # mesh
    -i "test/generated_interface.obj" # interface triangulation
    -o "test/cube_inserted.mesh"      # output
    -p "-1.0"                         # extent of the plane along the axis
    -x "Z"                            # axis the plane aligns to
    -e "0.00001"                      # epsilon around the extent of the axis
```


### 4. Merging
Takes two arbitrary tetrahedral meshes, as well as their interface-metadata, and merges them into the first mesh.

<b>Implementation pending.</b>
