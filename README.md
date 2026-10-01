# AlternatingProjections.jl

Part of the [Phase.jl](https://github.com/JuliaPhase/Phase.jl) ecosystem.

<!-- DOI badge: add after first Zenodo release -->

| **Documentation**                                                               |                                                                          
|:-------------------------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | 

_A small project started as test of whether it is indeed as easy as advertised to translate algorithms from their mathematical description to Julia language._

Several simple alternating projection algorithms are planned for implementation:
- Gerchberg-Saxton (GS)
- TIP
- Vector GS
- DRAP

The goals are:
 - to learn Julia
 - to try to keep the Julia implementation as close to the pseudocode description as possible
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
0.1

### Owner
Oleg Soloviev, o.a.soloviev@tudelft.nl

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://juliaphase.github.io/AlternatingProjections.jl/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://juliaphase.github.io/AlternatingProjections.jl/latest/

## Funding

This work has received funding from the ECSEL Joint Undertaking (JU) under grant agreement No 826589 (MADEin4). The JU receives support from the European Union's Horizon 2020 research and innovation programme and France, Germany, Austria, Italy, Sweden, Netherlands, Belgium, Hungary, Romania and Israel.

This work is part of the 14AMI project (grant agreement No 101111948). The project is supported by the Chips Joint Undertaking and its members including the top-up funding by RVO (The Netherlands Enterprise Agency).

<img src="docs/src/assets/funding/EU-flag.svg" alt="European Union flag" height="60">
<img src="docs/src/assets/funding/ECSEL-JU.jpg" alt="ECSEL Joint Undertaking" height="60">
<img src="docs/src/assets/funding/Chips-JU.png" alt="Chips Joint Undertaking, co-funded by the European Union" height="60">
