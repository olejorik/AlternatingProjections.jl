# AlternatingProjections.jl

A package implementing Alternating Projections methods in Julia.

## Alternating projections (background)

Alternating Projections (AP) is a simple method of solving a feasibility problem
```math
\text{for sets } A, B    \text{ find } x \in A \cap B
``` 
if ``A\cap B \neq \varnothing . `` 
In case of inconsistent problem ``A\cap B = \varnothing``, AP finds ``x \in A`` 
closest to ``B`` in some sense.  

The algorithm starts with some initial value ``x_0`` and proceeds by projecting
``x^k`` alternatively on ``A`` and ``B``
:
```math
y^{k} = P_B(x^k).
```
and
```math
x^{k+1} = P_A(y^k).
```

The method was originally proposeby von Neumann in 1949 for ``A`` and ``B`` being
linear subspaces, and was later generalised for any number of convex sets ``A_1,
\ldots,A_N`` as
```math
x^{k+1} = P_{A_{n(k)}}(x^k). 
```

It was also recently shown that the method can also work for not convex sets,
if the condition of _transversality_ is satisfied, (see for instance the Gerchberg-Saxton
method).

A lot of popular algorithms can be explained in AP framework, 
even if the original algorithm was invented base on other principles.
See, for instance, Gerchberg-Saxton algorithm for phase retrieval.

### Projections

Projection operator is often more easily described in math formula than in a computer program.
It is easy to project on a (multidimensional) sphere, for instance, or to a polyhedron.
The projection on a general convex set can be much more difficult (consider, for instance, projection on an ellipse).
This explains why in this package we limit ourselves to special types of sets that allows relatively easiness
 of defining a projection operator.
 
## Forward and backward operators
In some case the projection can be more easily computed in a transformed space (e.g. in the Fourier domain).
This can be generalised further by introducing some abstract (even not necassarily linear) forward and backward transforms
```math
x^{k+1} = b_{n(k)}(P_{A_{n(k)}}(f_{n(k)}(x^k))), 
```  
where ``f_n(x)`` and ``b_n(x)`` are the forward and backward transforms correspondingly.
For instance, in TIP algorithm ``f(x) = b(x) = i /_* x``, where ``/_*`` is ''deconvolution'' operator,
such that 
`` x /_* x = \delta_0``.

## Package features

- types for the frequently used feasible sets and projections on them
- a general AP algorithm and examples of its adaptation to popular AP algorithms, including:
    - Gerschberg-Saxton
    - Vector Gerschberg-Saxton
    - Gerschberg-Papoulis
    - TIP
    - Douglas–Rachford
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
