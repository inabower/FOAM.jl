# OpenFOAM handler with using Julia

- Julia : v1.0.1
- OpenFOAM : v1806

```
using Pkg

Pkg.add("IJulia")
Pkg.add("PyPlot")
Pkg.add("Gadfly")
Pkg.add("DataStructures")
```

## Tasks

- __done__ read case dir
- __done__ read time
  - __done__ time list
  - field
	- __done__ dimensions
	- __done__ internalField
	- boundaryField
- read constant
  - mesh
    - __done__ checkMesh
  - __done__ read regionProperties
- read system
