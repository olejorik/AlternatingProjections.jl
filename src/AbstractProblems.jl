"""
A collection of abstract types representing problems and methods for their solution
"""
module AbstractProblems

export Problem, Algorithm, IterativeAlgorithm, solve

"""
Abstract type representing a problem to solve
"""
abstract type Problem end

"""
Abstract type containing any algorithm
"""
abstract type Algorithm end



"""
Iterative methods form an important class of algorithms with iteratively adjust solution.

Concrete types of the `IterativeAlgorithm` should contain the initial value, tolerance and maximum number of iterations and 
are used more for convenience. These can be obtained by fucntions `initial`, `tolerance`, `maxit` fuctions. If applied the
abstract types, these fucntions return `missing` and can trigger using of the default values.

In addition, the concrete types contain instructions on whether or not to keeep the history of the convergence and 
some snapshots of the inner state of an iterator. These values are given by functions `keephistory` and `snapshots`.

"""
abstract type IterativeAlgorithm <: Algorithm end

initial(alg::IterativeAlgorithm) = missing
tolerance(alg::IterativeAlgorithm) = missing
maxit(alg::IterativeAlgorithm) = missing
keephistory(alg::IterativeAlgorithm) = false
snapshots(alg::IterativeAlgorithm) = Int64[]


"""
    solve(p::Problem,x⁰,alg::Algorithm)

    solve(p::Problem,alg::IterativeAlgorithm; x⁰, ϵ, maxit, keephistory, snapshshots)


Solve problem `p`, using method `alg`. For iterative algorithms the arguments may be specified separately.
Optionally keep the error history and the iteration snapshots.

"""
solve(p::Problem, alg::Algorithm) = error("Don't know how to solve ", typeof(p), " with method ", typeof(alg))

#unpack the values of iterative algorithm
solve(p::Problem, alg::IterativeAlgorithm; x⁰ = initial(alg), ϵ =tolerance(alg), maxit = maxit(alg), keephistory = keephistory(alg), snapshots = snapshots(alg)) = solve(p, alg, x⁰, ϵ, maxit, keephistory, snapshots)
# solve(p::Problem, alg::IterativeAlgorithm) = solve(p, alg, initial(alg), tolerance(alg), maxit(alg), keephistory(alg), snapshots(alg))

# Now proceed with the unpacked default or specified by the user values
function solve(p::Problem, alg::IterativeAlgorithm,  args...) 
    error("Don't know how to solve ", typeof(p), " with method ", typeof(alg))
end
    
end

# these types and fucntions will be modified by including file, so we import them.
import .AbstractProblems: Problem, Algorithm, IterativeAlgorithm, solve, initial, tolerance, maxit, keephistory, snapshots 
