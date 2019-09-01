# AlternatingProjections.jl

A package implementing Alternating Projections (AP) methods in Julia.

## Package features

- types for the frequently used feasible sets and projections on them
- a general AP algorithm and examples of its adaptation to popular AP algorithms, including:
    - Gerschberg-Saxton
    - Vector Gerschberg-Saxton
    - Gerschberg-Papoulis
    - TIP
    - Douglas-Rachfor
    - DRAP
    

## Manual outline

```@contents
```

### Functions

```@docs
project(x, feasset::FeasibleSet)
```

### Types

```@docs 
FeasibleSet
```

```@docs 
ConvexSet
```

```@docs 
AmplitudeConstraint
```

```@docs 
SupportConstraint
```

## Index

```@index
```
