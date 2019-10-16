var documenterSearchIndex = {"docs":
[{"location":"#AlternatingProjections.jl-1","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"A package implementing Alternating Projections methods in Julia.","category":"page"},{"location":"#Alternating-projections-(theoretical-background)-1","page":"AlternatingProjections.jl","title":"Alternating projections (theoretical background)","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Alternating Projections (AP) is a simple method of solving a feasibility problem","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"textfor sets  A B    text find  x in A cap B","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"if Acap B neq varnothing   In case of inconsistent problem Acap B = varnothing, AP finds x in A  closest to B in some sense.  ","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"The algorithm starts with some initial value x^0 and proceeds by projecting x^k alternatively on A and B :","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"y^k = P_B(x^k)","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"and","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"x^k+1 = P_A(y^k)","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"The method was originally proposed by von Neumann in 1949 for A and B being linear subspaces, and was later generalised for any number of convex sets A_1 ldotsA_N as","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"x^k+1 = P_A_n(k)(x^k) ","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"It was also recently shown that the method can also work for not convex sets, if the condition of transversality is satisfied, (see for instance the Gerchberg-Saxton method).","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"A lot of popular algorithms can be explained in AP framework,  even if the original algorithm was invented based on other principles. See, for instance, Gerchberg-Saxton algorithm for phase retrieval.","category":"page"},{"location":"#Projections-1","page":"AlternatingProjections.jl","title":"Projections","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Projection operator is often more easily described in math formula than in a computer program. It is easy to project on a (multidimensional) sphere, for instance, or to a polyhedron. The projection on a general convex set can be much more difficult (consider, for instance, projection on an ellipse). This explains why in this package we limit ourselves to special types of sets that allows relatively easiness  of defining a projection operator.","category":"page"},{"location":"#Forward-and-backward-operators-1","page":"AlternatingProjections.jl","title":"Forward and backward operators","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"In some case the projection can be more easily computed in a transformed space (e.g. in the Fourier domain). This can be generalised further by introducing some abstract (even not necassarily linear) forward and backward transforms","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"x^k+1 = b_n(k)(P_A_n(k)(f_n(k)(x^k))) ","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"where f_n(k)(x) and b_n(k)(x) are the forward and backward transforms correspondingly. For instance, in TIP algorithm f(x) = b(x) = i _* x, where _* is ''deconvolution'' operator, such that  x _* x = delta_0.","category":"page"},{"location":"#AlternatingProjections.jl-(methodological-background)-1","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl (methodological background)","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"The overall auto-pedagogical and methodological goal of the package is to learn Julia and to check whether it is indeed as suitable for mathematical programming as it is advertised.  So the package is written with the idea to keep its implementation as close to the mathematical or pseudocode description of the algorithms as possible.","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"In mathematics, the abstract concepts given by some definition and restrictive properties are often used. Similar concepts are often united in a single more general concept and so on. This idea of concepts hierarchy is reflected in Julia's type system. The package, on example of the alternating projections (AP) framework, introduces the main concepts of AP as hierarchy of types. Some of these implementations are obvious, some look artificial, and might be improved later, with the overall goal not the efficiency, but closeness to the mathematical text.","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Consider, for instance, how you would explain to someone a method of alternating projections (see the previous chapter). By explaining it, you would mention concepts of sets, then you would mark some of them as feasible sets by requiring them to have some properties, for instance, being convex. Then you would explain what is the projection operation and how is it related to finding the closest point in the convex set. Then you would introduce the broad class of the alternating projections, specifying that here you will use its variant sometimes named as POCS (projections on the convex sets), and explain that the method proceeds iteratively and converges either to a feasible point or to a stable cycle of length 2 (for 2 sets). You would, most probably, mention that in practice the method is allowed to run for some fixed number of iterations or until the points get close enough to each other.","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"In this informal explanation one can extract however several important formal or abstract concepts.  These are of the feasible sets, convex sets as a particular case of feasible sets, concept of projection operator, abstract concept of a feasibility problem, and a concept of a method of alternating projections to solve it, with a particular subtype POCS and maybe a supertype of abstract \"algorithm\". These abstract and particular concepts can be implemented in Julia as abstract and concrete types, but you might think why does ne need to do it, as the algorithm itself is quite straightforward to program in most languages. Well, the idea to do implement all the concepts in an abstract way has first of all the goal of generalisation. When you write these concepts on paper, it is helpful to see analogies and/or refer one problem to another class of the problems, so in the proof of some properties you might concentrate only on the proof of some small details specific to this type of problem.","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"The same should work also in the programming languages, I hope. For instance, suppose you want now to explain to someone Gerchberg-Saxton algorithm for phase retrieval problem. You might simply describe four simple steps of the algorithm, and it would remain some magic for your listener, or you can formulate it in the framework of feasibility problem, with slightly different type of the feasible sets and slightly different projection operations. With the same easiness it should be possible to be programmed – we need just to  introduce another subtypes of the feasible sets and define the projection operation on them. And this is exactly wat the multiple dispatch feature of Julia does, so that's why I think it might be interesting to try to implement this in Julia. One additional reason, with several algorithms already implemented in one package, it would be easier to make quick experiments, like what if we use this algorithm to that sort of problem?  With proper hierarchies of abstract types for the feasible sets and the methods for the problem solutions, its should be easy.  ","category":"page"},{"location":"#Example-of-type-hierarchy:-FeasibleSet-1","page":"AlternatingProjections.jl","title":"Example of type hierarchy: FeasibleSet","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"This concept is more or less clear. There is an abstract type FeasibleSet with subtypes of different types, like ConvexSet and further like SupportConstrainedSet, AmplitudeConstrainedSet, etc. All these should be abstract types.","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Then we should be able to make some realisation  of these abstract sets.  The realisations should be implemented as concrete types, which cannot have any child types (the abstract types cannot have instantiations). This means that for an abstract type AmplitudeConstrainedSet we should have a concrete subtype, which should have some other name and which should contain some data related to a concrete implementation of the abstract concept. For instance, when in phase retrieval problem one specifies a set of all possible pupil fields with a given field amplitude a as  𝓐_a  =  x in ℂ^Ntimes N    x = a  we can say that this is a structure named, for instance AmplitudeConstrain(a), which is substantiation of a concrete type AmplitudeConstrain, which is, in turn, a subtype of abstract type AmplitudeConstrainedSet. This particular concrete type reflects just one of many possible ways of implementation of the amplitude constrain by specifying the point-wise absolute value of of a finite dimensional complex vector.","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"To better distinguish between the abstract and concrete types, in this package there is a convention to name the concrete types with other word order, namely a concrete subtype of abstract type AmplitudeConstrainedSet is named as ConstrainedByAmplitude(a).","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Then we can define a general phase retrieval problem as","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"X = F x and \r\nx \\in \\aset <: AmplitudeConstrainedSet,\r\nX \\in \\Aset <: AmplitudeConstrainedSet,","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"and some concrete problem we want to solve as","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"X = F x and \r\nx \\in \\aset <: ConstrainedByAmplitude(a),\r\nX \\in \\Aset <: ConstrainedByAmplitude(A).","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Now the phase retrieval problem can be seen as a feasibility problem","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"textfind  x in A cap B  A = 𝓐_A  B= F 𝓐_a","category":"page"},{"location":"#Problem-and-Algorithm-abstract-types-1","page":"AlternatingProjections.jl","title":"Problem  and Algorithm abstract types","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"The package also introduces a special abstract type Problem, with concrete types as FeasibilityProblem and so on.  This introduction is done not so much for reflecting an abstract concept (althought this seems to be captured by this type as well) as for easiness of use. For instance, to solve a feasibility problem of finding x in Acap F^-1 B with an AP method starting from some point x^0, with stopping criteria maxit and maxϵ one can write a function  apsolve(A, B, F, F⁻¹; x⁰=zeros(size(A)), maxit = 20, maxϵ =0.01), but it seems to be much more elegant to group the first four arguments in one concept of problem and the last two being parameters of AP methods:","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"struct FeasibilityProblem <: Problem\r\n    A::FeasibleSet\r\n    B::FeasibleSet\r\n    forward\r\n    backward\r\nend","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"and","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"struct AP <: APMethod\r\n    maxit\r\n    maxϵ\r\nend","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"and the function looks both neater and more general now:","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"function solve(p::Problem,x⁰,alg::APMethod)","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Thinking in the same logic as used when creating types for the feasible sets, we can create abstract types for the methods, something like APMethods :> POCS APMethods :> GS APMethods : DR and so on.","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"It is not however 100% clear what should belong to the concrete types. For instance, whether the backward and forward operations should belong to the method or to the set definition? On one hand, the definitions like ConstrainedBySupport(a) are valid only in some fixed basis, and thus the operations like forward and backward are only use to bring us to the basis where these operations are easier to perform (consider, for instacne, variant of TIP algorithm, where projection of the PSFs can be performed directly in the frequency domain, thus saving two fft operations). On the other hand, the operations like forward and inverse transforms in the GS algorithm, or, in much more extent, component-wise inversion of TIP algorithm, are binded to their names, and same GS without the Fourier transforms and TIP without inversions would be named just alternating projections. In addition, this simplifies the implementation of projections a little bit.","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"At the current moment, I have chosen to bind the forward and backward operations to the method definition, but the things might change with the future versions.","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"The initial and termination conditions (initerm) are not the part of the abstact method, but might be sought as part of a particular implementation. Then it is much like as with the feasible sets, and we can introduce the concrete types containing initerm.","category":"page"},{"location":"#Function-solve-1","page":"AlternatingProjections.jl","title":"Function solve","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Now it is not difficult to write the implementation of the function solve for a  feasibility problem using AP method:","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"function solve(p::FeasibilityProblem, x⁰, alg::AP)\r\n    A = p.A\r\n    B = p.B\r\n    forward = p.forward\r\n    backward = p.backward\r\n    maxit = alg.maxit\r\n    maxϵ =alg.maxϵ\r\n\r\n    k = 0\r\n    xᵏ = x⁰\r\n    ϵ = Inf\r\n\r\n    while k < maxit && ϵ > maxϵ\r\n        ỹᵏ = forward(xᵏ)\r\n        yᵏ = project(ỹᵏ, B)\r\n        x̃ᵏ⁺¹ = backward(yᵏ)\r\n        xᵏ⁺¹ = project(x̃ᵏ⁺¹, A)\r\n        ϵ = LinearAlgebra.norm(xᵏ⁺¹ - xᵏ)\r\n        xᵏ = xᵏ⁺¹\r\n    #         println(ϵ)\r\n        k += 1\r\n    end\r\n\r\n    println(\"To converge with $ϵ accuracy, it took me $k iterations\")\r\n    return xᵏ\r\n\r\nend","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Note that here the first 6 lines of the code are just to save me some typing;  the code itself repeats literally the algorithm description. So this actually proves the statement of easiness of transfer of the math algorithms  in Julia. ","category":"page"},{"location":"#Function-project-1","page":"AlternatingProjections.jl","title":"Function project","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Of course, to be able to use function solve from the previous section, we need to introduce the projections. And here, using the Julia's power of multiple dispatch, we can write teh implementations of projections only for the sets for which we know how to do it. This automatically implies the most general implementation as this:","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"function project(x, feasset::FeasibleSet)\r\n    error(\"Don't know how to project on \", typeof(feasset))\r\nend","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"This function does nothing but complains about its inabaility to do anything. But we can teach it to project on the ConstrainedByAmplitude, for instance:","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"function project(x, feasset::ConstrainedByAmplitude)\r\n    return feasset.amp .* exp.( im * angle.(x))\r\nend","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"As you can see, the projection on this set only replaces the absolute value by that of  the constraint and keeps the argument.","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Now we can check our just implemented AP method","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"y = zeros(ComplexF32,10,10)\r\ny[1:5,1:5] = randn(ComplexF32, 5,5)\r\nY = fft(y)\r\np= FeasibilityProblem(ConstrainedByAmplitude(abs.(y)),ConstrainedByAmplitude(abs.(Y)),fft, ifft)\r\ngs=AP(3000,1e-10)\r\nz= solve(p, zeros(size(y)),gs)\r\nabs.(z) ≈ abs.(y)","category":"page"},{"location":"#Workflow-and-package-structure-1","page":"AlternatingProjections.jl","title":"Workflow and package structure","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"As the second goal of this project is to become more familiar with the programming in Julia, I log here my steps to create the \"workspace\", like set up of IDE for the package and its documentation development and steps for creating and extending the package functionality. This log represents my personal and unexperienced experience, so it might be far from the ideal and not optimised.","category":"page"},{"location":"#Development-workflow-1","page":"AlternatingProjections.jl","title":"Development workflow","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"I'm working in PyCharm with Julia plugin installed. The tests are located in test\\runtests.jl file, which I can rerun using Shift+F10 shortcut. In the terminal window I have opened several terminals, one for compiling the documentation, another for starting Jupyter lab, another for running Julia repl.","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"I start Julia REPL in the terminal by julia command. It works because I have created a special folder D:\\Documents\\bin, added it to the Windows PATH variable and created there julia.bat file containing single line","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"\"D:\\Julia-1.2.0\\bin\\julia.exe\" %*","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Then in REPL (in package mode, press ] to enter it) I activate the current environment with ","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"activate .","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"and start Revise typing ","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"using Revise","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Then I try different things with the introduced by the package commands, modify the code, and save the best pieces in test files or in the demo notebooks or in the runtests.jl.","category":"page"},{"location":"#Documentation-workflow-1","page":"AlternatingProjections.jl","title":"Documentation workflow","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"in the second terminal window, in docs dirtectory, run ","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"julia make.jl","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"In the third terminal window, also in the docs folder, issue","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"python -m http.server --bind localhost","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Now the documentation can be seen in web browser by  http://127.0.0.1:8000/build/ url. Edit, save, run again make.jl, check the results in the browser. More info in the docs of Documenter.jl.","category":"page"},{"location":"#Package-structure-1","page":"AlternatingProjections.jl","title":"Package structure","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Here are the steps I used to create the package AlternatingProjections.  I have decided to generate its initial structure with PkgTemplate package.","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Check the git configuration for the presence of the required fields (see doc page of PkgTemplate)","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"shell> git config --global --list\r\n\r\n...\r\nuser.name=Oleg Soloviev\r\nuser.email=oleg.soloviev@gmail.com\r\n...\r\ngithub.user=olejorik","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"open REPL, navigate to Julia_learning folder","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"    julia> cd(\"..\\\\Documents\\\\Julia_learning\\\\\")","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Check installed packages (press ] to enter package mode)","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"(v1.1) pkg> st","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"I don't have PkgTemplates, so I add it","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"(v1.1) pkg> add PkgTemplates\r\njulia> using PkgTemplates","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"According to the description  https://github.com/invenia/PkgTemplates.jl, I make a template with the default values, but in the current directory current directory, and turning on SSH for git sync","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"julia> t = Template(; dir = \".\",ssh=true)\r\nTemplate:\r\n  → User: olejorik\r\n  → Host: github.com\r\n  → License: MIT (Oleg Soloviev <oleg.soloviev@gmail.com> 2019)\r\n  → Package directory: D:\\Documents\\Julia_learning\\\r\n  → Minimum Julia version: v1.0\r\n  → SSH remote: Yes\r\n  → Add packages to main environment: Yes\r\n  → Commit Manifest.toml: No\r\n  → Plugins: None","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"julia> generate(\"AlternatingProjections\", t)\r\nGenerating project AlternatingProjections:\r\n    D:\\Documents\\Julia_learning\\AlternatingProjections\\Project.toml\r\n    D:\\Documents\\Julia_learning\\AlternatingProjections\\src/AlternatingProjections.jl\r\n[ Info: Initialized Git repo at D:\\Documents\\Julia_learning\\AlternatingProjections\r\n[ Info: Set remote origin to git@github.com:olejorik/AlternatingProjections.jl.git\r\nResolving package versions...\r\n  Updating `D:\\Documents\\Julia_learning\\AlternatingProjections\\Project.toml`\r\n  [8dfed614] + Test\r\n  Updating `D:\\Documents\\Julia_learning\\AlternatingProjections\\Manifest.toml`\r\n  [2a0f44e3] + Base64\r\n  [8ba89e20] + Distributed\r\n  [b77e0a4c] + InteractiveUtils\r\n  [56ddb016] + Logging\r\n  [d6f4376e] + Markdown\r\n  [9a3f8284] + Random\r\n  [9e88b42a] + Serialization\r\n  [6462fe0b] + Sockets\r\n  [8dfed614] + Test\r\n  Updating registry at `C:\\Users\\Oleg\\.julia\\registries\\General`\r\n  Updating git-repo `https://github.com/JuliaRegistries/General.git`\r\nResolving package versions...\r\n  Updating `D:\\Documents\\Julia_learning\\AlternatingProjections\\Project.toml`\r\n[no changes]\r\n  Updating `D:\\Documents\\Julia_learning\\AlternatingProjections\\Manifest.toml`\r\n  [2a0f44e3] - Base64\r\n  [8ba89e20] - Distributed\r\n  [b77e0a4c] - InteractiveUtils\r\n  [56ddb016] - Logging\r\n  [d6f4376e] - Markdown\r\n  [9a3f8284] - Random\r\n  [9e88b42a] - Serialization\r\n  [6462fe0b] - Sockets\r\n  [8dfed614] - Test\r\n[ Info: Committed 6 files/directories: src/, Project.toml, test/, README.md, LICENSE, .gitignore\r\nResolving package versions...\r\n  Updating `C:\\Users\\Oleg\\.julia\\environments\\v1.1\\Project.toml`\r\n  [37589ef0] + AlternatingProjections v0.1.0 [`D:\\Documents\\Julia_learning\\AlternatingProjections`]\r\n  Updating `C:\\Users\\Oleg\\.julia\\environments\\v1.1\\Manifest.toml`\r\n  [37589ef0] + AlternatingProjections v0.1.0 [`D:\\Documents\\Julia_learning\\AlternatingProjections`]\r\n[ Info: New package is at D:\\Documents\\Julia_learning\\AlternatingProjections\r\n","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"Now we can cd to the project directory and activate it","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"julia> cd(\"AlternatingProjections\\\\\")\r\n(v1.1) pkg> activate .","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"The status package shows installed packages","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"(AlternatingProjections) pkg> st\r\nProject AlternatingProjections v0.1.0\r\n    Status `D:\\Documents\\Julia_learning\\AlternatingProjections\\Project.toml`\r\n  (no changes since last commit)","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"In principle, I could run now ","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"] add FFTW","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"but it's installed in the default environment, so I'm not sure whether it is needed or not. I will do it for reproductibility and to add it to the required packages","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"(AlternatingProjections) pkg> add FFTW\r\n\r\n(AlternatingProjections) pkg> st\r\nProject AlternatingProjections v0.1.0\r\n    Status `D:\\Documents\\Julia_learning\\AlternatingProjections\\Project.toml`\r\n  [7a1cc6ca] + FFTW v0.3.0\r\n    Status `D:\\Documents\\Julia_learning\\AlternatingProjections\\Manifest.toml`\r\n  [7a1cc6ca] + FFTW v0.3.0\r\n","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"I will open the directory in Atom-Juno and/or in Pycharm. The former has nice integration with Julia, the latter has powerful editing features. In PyCharm, set the project interpreter to \"none\" (it's not a Python project)","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"In Juno's Repl I activate the environment, load Revise package and load the just created AlternatingProjection package","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"] activate.\r\n\r\njulia> using Revise\r\n[ Info: Recompiling stale cache file C:\\Users\\Oleg\\.julia\\compiled\\v1.1\\Revise\\M1Qoh.ji for Revise [295af30f-e4ad-537b-8983-00126c2a3abe]\r\n\r\n\r\njulia> using AlternatingProjections\r\n[ Info: Precompiling AlternatingProjections [37589ef0-cc9c-11e9-2e35-d5d8b7f5a2df]","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"I can check that the package is loaded by the greet function (not exported by default)","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"julia> AlternatingProjections.greet()\r\nHello World!","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"To check the Revise package tracks the changes in the code, let's export the greet function by adding in AlternatingProjections.jl a line after the first line (module ..)","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"export greet","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"and save the file. Now REPL should see greet function","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"julia> greet()\r\nHello World!\r\njulia> @which greet\r\nAlternatingProjections","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"The package is ready for the development.","category":"page"},{"location":"#Package-features-1","page":"AlternatingProjections.jl","title":"Package features","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"types for the frequently used feasible sets and projections on them\na general AP algorithm and examples of its adaptation to popular AP algorithms, including:[1]\nGerchberg-Saxton\nVector Gerchberg-Saxton\nGerchberg-Papoulis\nTIP\nDouglas–Rachford\nDRAP","category":"page"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"[1]: not implemented features are shown with italics","category":"page"},{"location":"#Manual-outline-1","page":"AlternatingProjections.jl","title":"Manual outline","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"","category":"page"},{"location":"#Functions-1","page":"AlternatingProjections.jl","title":"Functions","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"project(x, feasset::FeasibleSet)","category":"page"},{"location":"#AlternatingProjections.project-Tuple{Any,FeasibleSet}","page":"AlternatingProjections.jl","title":"AlternatingProjections.project","text":"project(x, A)\n\nProject x on set A.\n\n\n\n\n\n","category":"method"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"solve( p::FeasibilityProblem, x⁰, alg::AP)","category":"page"},{"location":"#AlternatingProjections.solve-Tuple{FeasibilityProblem,Any,AP}","page":"AlternatingProjections.jl","title":"AlternatingProjections.solve","text":"solve(p::Problem,x⁰,alg::APMethod)\n\nExamples\n\njulia>\n\n\n\n\n\n","category":"method"},{"location":"#Types-1","page":"AlternatingProjections.jl","title":"Types","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"FeasibleSet","category":"page"},{"location":"#AlternatingProjections.FeasibleSet","page":"AlternatingProjections.jl","title":"AlternatingProjections.FeasibleSet","text":"FeasibleSet\n\nAbstract type representing a set for feasibility problem.\n\n\n\n\n\n","category":"type"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"ConvexSet","category":"page"},{"location":"#AlternatingProjections.ConvexSet","page":"AlternatingProjections.jl","title":"AlternatingProjections.ConvexSet","text":"ConvexSet\n\nGeneral type, no projection method is specified.\n\nExamples\n\n\n\n\n\n\n\n","category":"type"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"AmplitudeConstrainedSet","category":"page"},{"location":"#AlternatingProjections.AmplitudeConstrainedSet","page":"AlternatingProjections.jl","title":"AlternatingProjections.AmplitudeConstrainedSet","text":"type AmplitudeConstraint\n\nExamples\n\njulia>\n\n\n\n\n\n","category":"type"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"ConstrainedBySupport","category":"page"},{"location":"#AlternatingProjections.ConstrainedBySupport","page":"AlternatingProjections.jl","title":"AlternatingProjections.ConstrainedBySupport","text":"ConstrainedBySupport(support)\n\nSpecial type of convex set.\n\nFor continuous case: consists of all functions that equals zero outside some fixed area called support: 𝓐_S = f  f(x) = 0 x  S  \n\nFor discrete case: all arrays that equals zero for indexes outside some index set: 𝓐_S = x  xi = 0 i  S  \n\nCurrently supports only discrete case, with the support defined as a boolean array.\n\nJulia version: \nAuthor: Oleg Soloviev\nDate: 2019-09-01\n\nExamples\n\n\njulia> S =  ConstrainedBySupport([true, false,true])\nConstrainedBySupport(Bool[true, false, true])\n\njulia> x = [1, 2, 3]; project(x, S)\n3-element Array{Int64,1}:\n 1\n 0\n 3\n\njulia> S = ConstrainedBySupport([x^2 + y^2 <=1  for x in -2:.2:2, y in -2:.2:2]);\njulia> x = ones(size(S.support));\njulia> project(x,S)\n21×21 Array{Float64,2}:\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0\n\n\n\n\n\n\n","category":"type"},{"location":"#Index-1","page":"AlternatingProjections.jl","title":"Index","text":"","category":"section"},{"location":"#","page":"AlternatingProjections.jl","title":"AlternatingProjections.jl","text":"","category":"page"}]
}