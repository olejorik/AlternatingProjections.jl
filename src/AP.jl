"""
    AP <: ProjectionsMethod

Classical Alternating Projections method

Project on one set, than on another, stop after fixed number of iterations or if the desired accuracy is achieved.
"""
abstract type AP <: ProjectionsMethod end

Base.@kwdef struct APparam <: AP
    x⁰ = missing
    maxϵ::Union{Float64,Missing} = missing
    maxit::Union{Missing,Int64} = missing
    keephistory::Bool = false
    snapshots::Array{Int64} = Int64[]
end

initial(alg::APparam) = alg.x⁰
tolerance(alg::APparam) = alg.maxϵ
maxit(alg::APparam) = alg.maxit
keephistory(alg::APparam) = alg.keephistory
snapshots(alg::APparam) = alg.snapshots

function solve(
    p::TwoSetsFP,
    alg::APparam,
    x⁰,
    maxϵ,
    maxit,
    keephistory::Bool,
    snapshots::Vector{Int64};
    gtfun=nothing,
)
    A = p.A
    B = p.B

    # process the default parameters
    !ismissing(x⁰) || (x⁰ = getelement(A))
    !ismissing(maxϵ) || (maxϵ = 1e-15)
    !ismissing(maxit) || (maxit = 100)

    k = 0
    ϵ = Inf

    xᵏ = copy(x⁰)
    xᵏ⁺¹ = similar(xᵏ)
    x̃ᵏ⁺¹ = similar(xᵏ)
    ỹᵏ = similar(xᵏ)
    yᵏ = similar(xᵏ)

    #  preallocated error vector
    err = similar(xᵏ)

    if keephistory
        errhist = Vector{Float64}(undef, maxit)
        disthist = Vector{Float64}(undef, maxit)
        if !isnothing(gtfun)
            distgthist = Vector{Float64}(undef, maxit)
        else
            distgthist = Float64[]
        end
    else
        errhist = Float64[]
        disthist = Float64[]
        distgthist = Float64[]
    end

    if length(snapshots) != 0
        xhist = [copy(x⁰) for i in 1:length(snapshots)]
    else
        xhist = typeof(x⁰)[]
    end
    j = 1

    while k < maxit && ϵ > maxϵ
        project!(yᵏ, xᵏ, B)
        project!(xᵏ⁺¹, yᵏ, A)

        err .= xᵏ⁺¹ .- xᵏ # This doesn't say much in infeasible case, but is OK in case of binary aperture
        # dist .= xᵏ⁺¹ .- yᵏ # this calculates true error but can stay large in case of infeasible case
        ϵ = LinearAlgebra.norm(err) / LinearAlgebra.norm(xᵏ) #relative norm
        xᵏ .= xᵏ⁺¹
        k += 1

        #         println(ϵ)
        if keephistory
            errhist[k] = ϵ
            disthist[k] = LinearAlgebra.norm(xᵏ⁺¹ .- yᵏ) / LinearAlgebra.norm(xᵏ) #relative norm
            if !isnothing(gtfun)
                distgthist[k] = gtfun(xᵏ⁺¹)
            end
        end

        if k ∈ snapshots
            println("Saving snapshot # $j, iteration # $k")
            xhist[j] .= xᵏ
            j += 1
        end
    end

    println(
        "Using $(supertype(typeof(alg))):  to converge with $ϵ accuracy, it took me $k iterations",
    )
    if keephistory
        @info "The distance between the sets at the solution point is $(disthist[k])"
    end
    return xᵏ,
    (
        lasty=yᵏ,
        errhist=errhist,
        xhist=xhist[1:(j - 1)],
        disthist=disthist,
        distgthist=distgthist,
        k=k,
    )
end

export AP, APparam
