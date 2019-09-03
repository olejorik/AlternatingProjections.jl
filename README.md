# AlternatingProjections.jl

_A small project started as test of whether it is ineed as easy as advertised to translate algorithms from their mathematical description to Julia language._

Several simple alternating projection algorithms are planned for implementation:
- Gerchberg-Saxton (GS)
- TIP
- Vector GS
- DRAP

The goals are:
 - learn Julia
 -  try to keep the Julia implementation as close to the pseudocode description as possible
 - to unify several algorithms in one framework.
 
 ## Usage
 
 1. Download the project in a folder named AlternatingProjections.
 2. Start Julia in this folder.
 3. Activate the package in Julia:
     ```
     ] activate .
    ```
4. Now you can use the package by 
    ```
    using AlternatingProjections
   ```
## Documentation

To build documentation, cd to `docs` folder and run `julia make.jl`

## Info
### Version
Permanently alpha?

### Owner
Oleg Soloviev, o.a.soloviev@tudelft.nl