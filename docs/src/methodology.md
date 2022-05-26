This is a short transcript of a [talk](https://olejorik.github.io/talk/generic-programming-in-julia/) given by me to my colleagues, explaining why I have worked on this package.
It's kept as a blog item.

# AlternatingProjections.jl (methodological background)

The overall auto-pedagogical and methodological goal of the package is to learn Julia and to check whether it is indeed as suitable for mathematical programming as it is advertised. 
So the package is written with the idea to keep its implementation as close to the mathematical or pseudocode description of the algorithms as possible.

In mathematics, the abstract concepts given by some definition and restrictive properties are often used.
Similar concepts are often united in a single more general concept and so on.
This idea of concepts hierarchy is reflected in Julia's *type system.*
The package, on example of the alternating projections (AP) framework, introduces the main concepts of AP as hierarchy of types.
Some of these implementations are obvious, some look artificial, and might be improved later, with the overall goal not the efficiency, but closeness to the mathematical text.

Consider, for instance, how you would explain to someone a method of alternating projections (see the previous chapter).
By explaining it, you would mention concepts of sets, then you would mark some of them as feasible sets by requiring them to have some properties, for instance, being convex.
Then you would explain what is the projection operation and how is it related to finding the closest point in the convex set.
Then you would introduce the broad class of the alternating projections, specifying that here you will use its variant sometimes named as POCS (projections on the convex sets), and explain that the method proceeds iteratively and converges either to a feasible point or to a stable cycle of length 2 (for 2 sets).
You would, most probably, mention that in practice the method is allowed to run for some fixed number of iterations or until the points get close enough to each other.

In this informal explanation one can extract however several important formal or abstract concepts. 
These are of the feasible sets, convex sets as a particular case of feasible sets, concept of projection operator, abstract concept of a feasibility problem, and a concept of a method of alternating projections to solve it, with a particular subtype POCS and maybe a supertype of abstract "algorithm".
These abstract and particular concepts can be implemented in Julia as abstract and concrete types, but you might think why does ne need to do it, as the algorithm itself is quite straightforward to program in most languages.
Well, the idea to do implement all the concepts in an abstract way has first of all the goal of generalisation.
When you write these concepts on paper, it is helpful to see analogies and/or refer one problem to another class of the problems, so in the proof of some properties you might concentrate only on the proof of some small details specific to this type of problem.

The same should work also in the programming languages, I hope.
For instance, suppose you want now to explain to someone Gerchberg-Saxton algorithm for phase retrieval problem.
You might simply describe four simple steps of the algorithm, and it would remain some magic for your listener, or you can formulate it in the framework of feasibility problem, with slightly different type of the feasible sets and slightly different projection operations.
With the same easiness it should be possible to be programmed -- we need just to  introduce another subtypes of the feasible sets and define the projection operation on them.
And this is exactly wat the multiple dispatch feature of Julia does, so that's why I think it might be interesting to try to implement this in Julia.
One additional reason, with several algorithms already implemented in one package, it would be easier to make quick experiments, like what if we use this algorithm to that sort of problem? 
With proper hierarchies of abstract types for the feasible sets and the methods for the problem solutions, its should be easy.  

## Example of type hierarchy: `FeasibleSet`

This concept is more or less clear. There is an abstract type FeasibleSet with subtypes of different types, like `ConvexSet` and further like `SupportConstrainedSet`, `AmplitudeConstrainedSet`, *etc.*
All these should be abstract types.

Then we should be able to make some realisation  of these abstract sets. 
The realisations should be implemented as concrete types, which cannot have any child types (the abstract types cannot have instantiations).
This means that for an abstract type AmplitudeConstrainedSet we should have a concrete subtype, which should have some other name and which should contain some data related to a concrete implementation of the abstract concept.
For instance, when in phase retrieval problem one specifies a set of all possible pupil fields with a *given* field amplitude ``a`` as 
`` 𝓐_a  = \{ x \in ℂ^{N\times N} \ | \ |x| = a \} ``
we can say that this is a structure named, for instance `AmplitudeConstrain(a)`, which is substantiation of a concrete type `AmplitudeConstrain`, which is, in turn, a subtype of abstract type `AmplitudeConstrainedSet`.
This particular concrete type reflects just one of many possible ways of implementation of the amplitude constrain by specifying the point-wise absolute value of of a finite dimensional complex vector.


To better distinguish between the abstract and concrete types, in this package there is a convention to name the concrete types with other word order, namely a concrete subtype of abstract type `AmplitudeConstrainedSet` is named as `ConstrainedByAmplitude(a)`.

Then we can define a general phase retrieval problem as
```
X = F x and 
x \in \aset <: AmplitudeConstrainedSet,
X \in \Aset <: AmplitudeConstrainedSet,
```
and some concrete problem we want to solve as
```
X = F x and 
x \in \aset <: ConstrainedByAmplitude(a),
X \in \Aset <: ConstrainedByAmplitude(A).
```

Now the phase retrieval problem can be seen as a feasibility problem
```math
\text{find } x \in A \cap B, \ A = 𝓐_A, \ B= F 𝓐_a
```



## `Problem`  and `Algorithm` abstract types


The package also introduces a special abstract type `Problem`, with concrete types as `FeasibilityProblem` and so on. 
This introduction is done not so much for reflecting an abstract concept (althought this seems to be captured by this type as well) as for easiness of use.
For instance, to solve a feasibility problem of finding ``x \in A\cap F^{-1} B`` with an AP method starting from some point ``x^0``, with stopping criteria `maxit` and `maxϵ`
one can write a function 
`apsolve(A, B, F, F⁻¹; x⁰=zeros(size(A)), maxit = 20, maxϵ =0.01)`,
but it seems to be much more elegant to group the first four arguments in one concept of *problem* and the last two being parameters of AP methods:
```julia
struct FeasibilityProblem <: Problem
    A::FeasibleSet
    B::FeasibleSet
    forward
    backward
end
```
and
```julia
struct AP <: APMethod
    maxit
    maxϵ
end
```
and the function looks both neater and more general now:
```julia
function solve(p::Problem,x⁰,alg::APMethod)
```

Thinking in the same logic as used when creating types for the feasible sets, we can create abstract types for the methods, something like
`APMethods :> POCS`
`APMethods :> GS`
`APMethods : DR`
and so on.

It is not however 100% clear what should belong to the concrete types.
For instance, whether the backward and forward operations should belong to the method or to the set definition?
On one hand, the definitions like `ConstrainedBySupport(a)` are valid only in some fixed basis, and thus the operations like forward and backward are only use to bring us to the basis where these operations are easier to perform (consider, for instacne, variant of TIP algorithm, where projection of the PSFs can be performed directly in the frequency domain, thus saving two `fft` operations).
On the other hand, the operations like forward and inverse transforms in the GS algorithm, or, in much more extent, component-wise inversion of TIP algorithm, are binded to their names, and same GS without the Fourier transforms and TIP without inversions would be named just alternating projections.
In addition, this simplifies the implementation of projections a little bit.

At the current moment, I have chosen to bind the forward and backward operations to the method definition, but the things might change with the future versions.

The initial and termination conditions (`initerm`) are not the part of the abstact method, but might be sought as part of a particular implementation.
Then it is much like as with the feasible sets, and we can introduce the concrete types containing `initerm`.

## Function `solve`

Now it is not difficult to write the implementation of the function `solve` for a 
feasibility problem using AP method:
```julia
function solve(p::FeasibilityProblem, x⁰, alg::AP)
    A = p.A
    B = p.B
    forward = p.forward
    backward = p.backward
    maxit = alg.maxit
    maxϵ =alg.maxϵ

    k = 0
    xᵏ = x⁰
    ϵ = Inf

    while k < maxit && ϵ > maxϵ
        ỹᵏ = forward(xᵏ)
        yᵏ = project(ỹᵏ, B)
        x̃ᵏ⁺¹ = backward(yᵏ)
        xᵏ⁺¹ = project(x̃ᵏ⁺¹, A)
        ϵ = LinearAlgebra.norm(xᵏ⁺¹ - xᵏ)
        xᵏ = xᵏ⁺¹
    #         println(ϵ)
        k += 1
    end

    println("To converge with $ϵ accuracy, it took me $k iterations")
    return xᵏ

end
```
Note that here the first 6 lines of the code are just to save me some typing; 
the code itself repeats literally the algorithm description.
So this actually proves the statement of easiness of transfer of the math algorithms 
in Julia. 


## Function `project`

Of course, to be able to use function `solve` from the previous section, we need
to introduce the projections.
And here, using the Julia's power of multiple dispatch, we can write teh implementations
of projections only for the sets for which we know how to do it.
This automatically implies the most general implementation as this:
```julia
function project(x, feasset::FeasibleSet)
    error("Don't know how to project on ", typeof(feasset))
end
```
This function does nothing but complains about its inabaility to do anything.
But we can teach it to project on the `ConstrainedByAmplitude`, for instance:
```julia
function project(x, feasset::ConstrainedByAmplitude)
    return feasset.amp .* exp.( im * angle.(x))
end
```
As you can see, the projection on this set only replaces the absolute value by that of
 the constraint and keeps the argument.
  
Now we can check our just implemented AP method
```julia
y = zeros(ComplexF32,10,10)
y[1:5,1:5] = randn(ComplexF32, 5,5)
Y = fft(y)
p= FeasibilityProblem(ConstrainedByAmplitude(abs.(y)),ConstrainedByAmplitude(abs.(Y)),fft, ifft)
gs=AP(3000,1e-10)
z= solve(p, zeros(size(y)),gs)
abs.(z) ≈ abs.(y)
```

## Workflow and package structure

As the second goal of this project is to become more familiar with the programming in Julia, I log here my steps to create the "workspace", like set up of IDE for the package and its documentation development and steps for creating and extending the package functionality.
This log represents my personal and unexperienced experience, so it might be far from the ideal and not optimised.

### Development workflow

I'm working in PyCharm with Julia plugin installed.
The tests are located in `test\runtests.jl` file, which I can rerun
using Shift+F10 shortcut.
In the terminal window I have opened several terminals, one for compiling the documentation, another for starting Jupyter lab, another for running Julia repl.

I start Julia REPL in the terminal by julia command. It works because I have created a special folder `D:\Documents\bin`, added it to the Windows PATH variable and created there `julia.bat` file containing single line

```commandline
"D:\Julia-1.2.0\bin\julia.exe" %*
```

Then in REPL (in package mode, press `]` to enter it) I activate the current environment with 

```julia
activate .
```

and start Revise typing 

```julia
using Revise
```

Then I try different things with the introduced by the package commands, modify the code, and save the best pieces in test files or in the demo notebooks or in the `runtests.jl`.

### Documentation workflow
in the second terminal window, in docs dirtectory, run 
```commandline
julia make.jl
```

In the third terminal window, also in the `docs` folder, issue
```commandline
python -m http.server --bind localhost
```

Now the documentation can be seen in web browser by  http://127.0.0.1:8000/build/ url.
Edit, save, run again `make.jl`, check the results in the browser.
More info in the [docs of Documenter.jl.](https://juliadocs.github.io/Documenter.jl/stable/man/guide/)

### Package structure

Here are the steps I used to create the package `AlternatingProjections`. 
I have decided to generate its initial structure with [`PkgTemplate` package](https://github.com/invenia/PkgTemplates.jl).


0. Check the git configuration for the presence of the required fields (see [doc page](https://invenia.github.io/PkgTemplates.jl/stable/#Usage-1) of `PkgTemplate`)
  
```commandline
shell> git config --global --list

...
user.name=Oleg Soloviev
user.email=oleg.soloviev@gmail.com
...
github.user=olejorik
```


1. open REPL, navigate to `Julia_learning` folder
```commandline
    julia> cd("..\\Documents\\Julia_learning\\")
```
 
2. Check installed packages (press `]` to enter package mode)
```commandline
(v1.1) pkg> st
```
I don't have `PkgTemplates`, so I add it
```commandline
(v1.1) pkg> add PkgTemplates
julia> using PkgTemplates
```

According to the description  https://github.com/invenia/PkgTemplates.jl, I make a template with the default values, but in the current directory current directory, and turning on SSH for git sync
```commandline
julia> t = Template(; dir = ".",ssh=true)
Template:
  → User: olejorik
  → Host: github.com
  → License: MIT (Oleg Soloviev <oleg.soloviev@gmail.com> 2019)
  → Package directory: D:\Documents\Julia_learning\
  → Minimum Julia version: v1.0
  → SSH remote: Yes
  → Add packages to main environment: Yes
  → Commit Manifest.toml: No
  → Plugins: None
```
```commandline
julia> generate("AlternatingProjections", t)
Generating project AlternatingProjections:
    D:\Documents\Julia_learning\AlternatingProjections\Project.toml
    D:\Documents\Julia_learning\AlternatingProjections\src/AlternatingProjections.jl
[ Info: Initialized Git repo at D:\Documents\Julia_learning\AlternatingProjections
[ Info: Set remote origin to git@github.com:olejorik/AlternatingProjections.jl.git
Resolving package versions...
  Updating `D:\Documents\Julia_learning\AlternatingProjections\Project.toml`
  [8dfed614] + Test
  Updating `D:\Documents\Julia_learning\AlternatingProjections\Manifest.toml`
  [2a0f44e3] + Base64
  [8ba89e20] + Distributed
  [b77e0a4c] + InteractiveUtils
  [56ddb016] + Logging
  [d6f4376e] + Markdown
  [9a3f8284] + Random
  [9e88b42a] + Serialization
  [6462fe0b] + Sockets
  [8dfed614] + Test
  Updating registry at `C:\Users\Oleg\.julia\registries\General`
  Updating git-repo `https://github.com/JuliaRegistries/General.git`
Resolving package versions...
  Updating `D:\Documents\Julia_learning\AlternatingProjections\Project.toml`
[no changes]
  Updating `D:\Documents\Julia_learning\AlternatingProjections\Manifest.toml`
  [2a0f44e3] - Base64
  [8ba89e20] - Distributed
  [b77e0a4c] - InteractiveUtils
  [56ddb016] - Logging
  [d6f4376e] - Markdown
  [9a3f8284] - Random
  [9e88b42a] - Serialization
  [6462fe0b] - Sockets
  [8dfed614] - Test
[ Info: Committed 6 files/directories: src/, Project.toml, test/, README.md, LICENSE, .gitignore
Resolving package versions...
  Updating `C:\Users\Oleg\.julia\environments\v1.1\Project.toml`
  [37589ef0] + AlternatingProjections v0.1.0 [`D:\Documents\Julia_learning\AlternatingProjections`]
  Updating `C:\Users\Oleg\.julia\environments\v1.1\Manifest.toml`
  [37589ef0] + AlternatingProjections v0.1.0 [`D:\Documents\Julia_learning\AlternatingProjections`]
[ Info: New package is at D:\Documents\Julia_learning\AlternatingProjections

```

 
3. Now we can cd to the project directory and activate it
```commandline
julia> cd("AlternatingProjections\\")
(v1.1) pkg> activate .
```
The status package shows installed packages
```commandline
(AlternatingProjections) pkg> st
Project AlternatingProjections v0.1.0
    Status `D:\Documents\Julia_learning\AlternatingProjections\Project.toml`
  (no changes since last commit)
```

In principle, I could run now 
```commandline
] add FFTW
```
but it's installed in the default environment, so I'm not sure whether it is needed or not. I will do it for reproductibility and to add it to the required packages

```commandline
(AlternatingProjections) pkg> add FFTW

(AlternatingProjections) pkg> st
Project AlternatingProjections v0.1.0
    Status `D:\Documents\Julia_learning\AlternatingProjections\Project.toml`
  [7a1cc6ca] + FFTW v0.3.0
    Status `D:\Documents\Julia_learning\AlternatingProjections\Manifest.toml`
  [7a1cc6ca] + FFTW v0.3.0

```

 
4.  I will open the directory in Atom-Juno and/or in Pycharm. The former has nice integration with Julia, the latter has powerful editing features. In PyCharm, set the project interpreter to "none" (it's not a Python project)

In Juno's Repl I activate the environment, load Revise package and load the just created AlternatingProjection package
    
```commandline
] activate.

julia> using Revise
[ Info: Recompiling stale cache file C:\Users\Oleg\.julia\compiled\v1.1\Revise\M1Qoh.ji for Revise [295af30f-e4ad-537b-8983-00126c2a3abe]


julia> using AlternatingProjections
[ Info: Precompiling AlternatingProjections [37589ef0-cc9c-11e9-2e35-d5d8b7f5a2df]
```

I can check that the package is loaded by the greet function (not exported by default)
```commandline
julia> AlternatingProjections.greet()
Hello World!
```

To check the Revise package tracks the changes in the code, let's export the greet function by adding in `AlternatingProjections.jl` a line after the first line (module ..)

```julia
export greet
```

and save the file. Now REPL should see `greet` function
```commandline
julia> greet()
Hello World!
julia> @which greet
AlternatingProjections
```
The package is ready for the development.
