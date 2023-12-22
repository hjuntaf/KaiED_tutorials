# KaiED_tutorials

-------------
- Contact : *Hyeong Jun LEE*   < hjuntaf _at_ gmail.com >
-------------

This is a set of tutorial codes for KaiED Package. (The `KaiED` Package will be released soon.)

## 0. Preparation

- Download `KaiED` and `Julia`. ([Go to `KaiED` repo.](https://github.com/KAIST-ELST/KaiED)) ([Download Julia here.](https://julialang.org/downloads/ "official webpage"))
  >  Version-test successfully done for
  > * julia >= 1.8.5

## 1. Files

- `src/` contains source codes.
- `tutorials/` contains tutorial codes.
- `envs/` is a temporary directory for Julia-project-environment TOML files. `setup.jl` will use this location and generate `envs/KED/*.toml` files.
- `setup.jl` will activate a 'project' of julia and install the required packages.

## 2. Usage

```
$ julia tutorials/01_manybody_binary_basis.jl
```
Every tutorial code here contains `include("../src/mybase.jl")` with relative path to include `mybase.jl` and use `KaiED`. You will need to modify this path (`../src/`) when you want to run in a different directory/location.


### Jupyter notebook

We also provide some tutorials within `jupyter notebook` (.ipynb).
You can use the jupyter notebooks after an installation of `IJulia` as follows.

#### Installation

```
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.9.3 (2023-08-24)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia> using Pkg

julia> Pkg.add("IJulia")
```
or enter the Pkg REPL by pressing `]` and install it :
```
(@v1.9) pkg> add IJulia
```

#### Run

```
$ jupyter notebook
```


## 3. Tutorial contents
The following concepts are included.
- Quantum many-body basis state manipulations
- Quantum many-body wavefunction manipulations
- Green functions
- Impurity models
- Tight-binding models
- Exact diagonalization method
- Dynamical mean-field theory
