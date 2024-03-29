"""
    DRAP <: ProjectionsMethod

Douglas-Rachford and Alternating Projection method

Combination of DR and AP method (see [1]).
Stop after fixed number of iterations or if the desired accuracy is achieved

[1] Nguyen Hieu Thao, Oleg Soloviev, and Michel Verhaegen. Convex combination of alternating
    projection and Douglas-Rachford operators for phase retrieval. Adv. Comput. Math.
"""
abstract type DRAP <: ProjectionsMethod end

Base.@kwdef struct DRAPparam <: DRAP
    x⁰ = missing
    maxϵ::Union{Float64, Missing} = missing
    maxit::Union{Missing,Int64} = missing
    keephistory::Bool = false
    snapshots::Array{Int64} = Int64[]
    β::Union{Float64, Missing} = missing
end

# DRAPparam() = DRAPparam(missing,missing,missing, false, Int64[],missing)

initial(alg::DRAPparam) = alg.x⁰
tolerance(alg::DRAPparam) = alg.maxϵ
maxit(alg::DRAPparam) = alg.maxit
keephistory(alg::DRAPparam) = alg.keephistory
snapshots(alg::DRAPparam) = alg.snapshots



function solve(p::TwoSetsFP, alg::DRAPparam, x⁰, maxϵ, maxit, keephistory::Bool, snapshots::Vector{Int64})
    A = p.A
    B = p.B

    # quick fix to be compatible with scaling-free problems
    nA = sqrt(sum(abs2,amp(A)))


    # process the default parameters
    !ismissing(x⁰) || ( x⁰ = getelement(A) )
    !ismissing(maxϵ) || (maxϵ = 1e-15)
    !ismissing(maxit) || (maxit = 100)

    β = alg.β
    !ismissing(β) || (β = 0.9)

    k = 0
    ϵ = Inf
    
    xᵏ = copy(x⁰)
    xᵏ⁺¹ = similar(xᵏ)
    x̃ᵏ⁺¹ = similar(xᵏ)
    ỹᵏ = similar(xᵏ)
    yᵏ = similar(xᵏ)
    zᵏ = similar(xᵏ)

    #  preallocated error vector
    err = similar(xᵏ)

    if keephistory
        errhist = Vector{Float64}(undef, maxit)
        disthist = Vector{Float64}(undef, maxit)
    else
        errhist = Float64[]
        disthist = Float64[]
    end

    if length(snapshots) != 0
        xhist = [copy(x⁰) for i in 1:length(snapshots)]
    else
        xhist = typeof(x⁰)[]
    end
    j = 1

    while k < maxit && ϵ > maxϵ

        # Tdrap = Pa( (1+β)Pb - β Id) - β(Pb -Id)
        project!(yᵏ, xᵏ, B) #Pb

        nPb = sqrt(sum(abs2,yᵏ))
        # yᵏ .*= (nA/nPb) #quick fix
        
        @. zᵏ = (1 + β) * yᵏ - β * xᵏ # (1+β)Pb - β Id
        project!(xᵏ⁺¹, zᵏ, A) # Pa( (1+β)Pb - β Id)
        @. xᵏ⁺¹ = xᵏ⁺¹ - β * (yᵏ - xᵏ)  # Pa( (1+β)Pb - β Id) - β(Pb -Id)

        # quick fix, to be changed as representative
        # removepiston
        cpiston = sum(xᵏ⁺¹)
        # cpiston = xᵏ⁺¹[1]
        cpiston /= abs(cpiston)
        # @info "Cpiston = $cpiston, angle = $(angle(cpiston))"
        xᵏ⁺¹ ./= cpiston

        err .= xᵏ⁺¹ .- xᵏ # This doesn't say much in infeasible case, but is OK in case of binary aperture
        # dist .= xᵏ⁺¹ .- yᵏ # this calculates true error but can stay large in case of infeasible case
        ϵ = LinearAlgebra.norm(err)
        xᵏ .= xᵏ⁺¹
        # xᵏ .= xᵏ⁺¹ * nA/sqrt(sum(abs2,xᵏ⁺¹)) # quick fix 
        k += 1

    #         println(ϵ)
        if keephistory
            errhist[k] = ϵ
            disthist[k] = LinearAlgebra.norm(xᵏ⁺¹ .- yᵏ)
        end

        if k ∈ snapshots
            println("Saving snapshot # $j, iteration # $k")
            xhist[j] .= xᵏ
             j += 1
        end

    end

    println("Using $(supertype(typeof(alg))): to converge with $ϵ accuracy, it took me $k iterations")
    if keephistory
        @info "The distance between the sets at the solution point is $(disthist[k])"
    end
    # if keephistory
    #     if length(snapshots) != 0
    #         return xᵏ, errhist, xhist
    #     else
    #         return xᵏ, errhist, Vector{T, (0,)}
    #     end
    # else
    #     return xᵏ, Vector{Float64, (0,)}, Vector{T, (0,)}
    # end
    return xᵏ, (lasty = yᵏ, errhist = errhist, xhist = xhist[1:j - 1], disthist = disthist, k= k)

end

export DRAP, DRAPparam


