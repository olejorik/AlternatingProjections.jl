"""
A collection of abstract types representing problems and methods for their solution
"""
module AbstractProblems

export Problem, Algorithm, IterativeAlgorithm, solve,
    solution, history, lasty, errhist, disthist, distgthist, xhist, itersteps, combinehist

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

"""
    initial(alg::IterativeAlgorithm)

Get the inital value of  iterative algorithm `alg`.
"""
initial(alg::IterativeAlgorithm) = missing

"""
    tolerance(alg::IterativeAlgorithm)::Float64

Get the tolerance value used as stopping criterium of iterative algorithm `alg`.
"""
tolerance(alg::IterativeAlgorithm) = missing

"""
    maxit(alg::IterativeAlgorithm)::Int

Get the maximum number of iterations used as stopping criterium of iterative algorithm `alg`.
"""
maxit(alg::IterativeAlgorithm) = missing

"""
    keephistory(alg::IterativeAlgorithm)::Bool

If set, iterative algorithm `alg` will keep the history of the error values (distances between subsequent states).
"""
keephistory(alg::IterativeAlgorithm) = false

"""
    snapshots(alg::IterativeAlgorithm)::Array{Int64}

Get the list of the iteration numbers, for which iterative algorithm `alg` will keep the state.
"""
snapshots(alg::IterativeAlgorithm) = Int64[]


"""
    solve(p::Problem,x⁰,alg::Algorithm)

    solve(p::Problem,alg::IterativeAlgorithm; x⁰, ϵ, maxit, keephistory, snapshshots)
    
    solve(p::Problem,(alg1, alg2,...); x⁰, ϵ, maxit, keephistory, snapshshots)


Solve problem `p`, using method `alg`. For iterative algorithms the arguments may be specified separately.
Optionally keep the error history and the iteration snapshots.

If sequence of the iterative algorithms is given, they are run subsequently, using the last state of the previous algorith as the inital value.

"""
solve(p::Problem, alg::Algorithm; kwargs...) = error("Don't know how to solve ", typeof(p), " with method ", typeof(alg))

#unpack the values of iterative algorithm
solve(p::Problem, alg::IterativeAlgorithm; x⁰=initial(alg), ϵ=tolerance(alg), maxit=maxit(alg), keephistory=keephistory(alg), snapshots=snapshots(alg), kwargs...) = solve(p, alg, x⁰, ϵ, maxit, keephistory, snapshots; kwargs...)
# solve(p::Problem, alg::IterativeAlgorithm) = solve(p, alg, initial(alg), tolerance(alg), maxit(alg), keephistory(alg), snapshots(alg))
#  unpack for several algorithms
solve(p::Problem, algs::Tuple{Vararg{IterativeAlgorithm}}; x⁰=initial(algs[1]), ϵ=tolerance(algs[1]), maxit=maxit(algs[1]), keephistory=keephistory(algs[1]), snapshots=snapshots(algs[1]), kwargs...) = solve(p, algs, x⁰, ϵ, maxit, keephistory, snapshots; kwargs...)
# for one tuple just call algorithm
solve(p::Problem, algs::Tuple{Vararg{IterativeAlgorithm,1}}; x⁰=initial(algs[1]), ϵ=tolerance(algs[1]), maxit=maxit(algs[1]), keephistory=keephistory(algs[1]), snapshots=snapshots(algs[1]), kwargs...) = solve(p, algs[1], x⁰, ϵ, maxit, keephistory, snapshots; kwargs...)

# Now proceed with the unpacked default or specified by the user values
function solve(p::Problem, alg::IterativeAlgorithm, args...; kwargs...)
    error("Don't know how to solve ", typeof(p), " with method ", typeof(alg))
end

# several algorithms, unpacked defaults 
function solve(p::Problem, algs::Tuple{Vararg{IterativeAlgorithm}}, args...; kwargs...)
    # solve with the first algorithm
    sol1 = solve(p, algs[1], args...; kwargs...)
    # if the error is smaller than the tolerane of the following algorithm, do nothing
    if last(errhist(sol1)) <= tolerance(algs[2])
        return sol1
    else
        # else continue starting with the last value and combine the solutions
        sol2 = solve(p, algs[2:end]; x⁰=copy(solution(sol1)), kwargs...)
        sol = (solution(sol2),
            # (lasty = lasty(sol2), 
            # errhist = [errhist(sol1); errhsit(sol2)],  
            # disthist = [disthist(sol1); disthist(sol2)], 
            # distgthist = [sol1[2][:distgthist]; sol2[2][:distgthist]], 
            # xhist = [sol1[2][:xhist];sol2[2][:xhist]],
            # k = sol1[2][:k] + sol2[2][:k]))
            combinehist(sol1, sol2)
        )
        return sol
    end
end

# getters for the solution, based on the form (solution, history)

function solution(sol::Tuple{Any,NamedTuple})
    return sol[1]
end

function history(sol::Tuple{Any,NamedTuple})
    return sol[2]
end

function lasty(sol::Tuple{Any,NamedTuple})
    return get(history(sol), :lasty, nothing)
end


function errhist(sol::Tuple{Any,NamedTuple})
    return get(history(sol), :errhist, nothing)
end


function disthist(sol::Tuple{Any,NamedTuple})
    return get(history(sol), :disthist, nothing)
end


function distgthist(sol::Tuple{Any,NamedTuple})
    return get(history(sol), :distgthist, nothing)
end

function xhist(sol::Tuple{Any,NamedTuple})
    return get(history(sol), :xhist, nothing)
end

function itersteps(sol::Tuple{Any,NamedTuple})
    return get(history(sol), :k, nothing)
end

function combinehist(sol1::Tuple{Any,NamedTuple}, sol2::Tuple{Any,NamedTuple})
    allkeys = union(keys(history(sol1)), keys(history(sol2)))
    allvalues = []
    for k in allkeys
        if k == :lasty
            push!(allvalues, lasty(sol2))
        elseif k == :k
            push!(allvalues, itersteps(sol1) + itersteps(sol2))
        else
            # push(allvalues, @eval [$k(sol1); $k(sol2)]) # this doesn't work
            push!(allvalues, [history(sol1)[k]; history(sol2)[k]]) # this is not using the functions I've defined above TODO correct
        end
    end
    return (; zip(allkeys, allvalues)...)
end


end




# these types and fucntions will be modified by parent module after including file, so we import them in the parent module.
import .AbstractProblems: Problem, Algorithm, IterativeAlgorithm, solve, initial, tolerance, maxit, keephistory, snapshots,
    solution, history, lasty, errhist, disthist, distgthist, xhist, itersteps, combinehist
