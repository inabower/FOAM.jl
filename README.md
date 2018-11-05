# OpenFOAM handler with using Julia

- Julia : v1.0.1
- OpenFOAM : v1806
- Python : 3.6.4 (anaconda 4.4.10)

```
pip install PyFoam

julia

using Pkg

Pkg.add("IJulia")
Pkg.add("PyPlot")
#ENV["PYTHON"] = /path/to/python/which/installed/PyFoam
#ENV["PYTHON"] = /home/user/.conda/env/py3/bin/python
#ENV["PYTHON"] = /home/user/Applications/anaconda3/bin/python
Pkg.add("PyCall")
Pkg.add("DataStructures")

exit()

export JULIA_NUM_THREADS=4

jupyter lab
```

# Example

![Example](src/example.png)

## Tasks

- __done__ read case dir
- __done__ read OpenFOAM dictFile
- __done__ read time
  - __done__ time list
  - field
	- __done__ dimensions
	- __done__ internalField
	- __ok__ boundaryField
- __ok__ read constant
  - __ok__ mesh
    - __done__ checkMesh
  - __done__ read regionProperties
- __ok__ read system
- __ok__ boundaryField
- __ok?__ parallel thread
- reconstruct and rm reconstructed times of decDir
- parseDict without pyFoam
