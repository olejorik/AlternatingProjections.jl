"""
DR <: ProjectionsMethod
Douglas-Rachford method

Reflect with respect one set, than to another and move there halfway: x → (x + R₂R₁x)/2
Stop after fixed number of iterations or if the desired accuracy is achieved
"""
abstract type DR <: ProjectionsMethod end

struct DRparam <: DR
x⁰
maxit::Union{Missing,Int64}
maxϵ::Union{Float64, Missing}
keephistory::Bool 
snapshots::Array{Int64}
end

DRparam() = DRparam(missing,missing,missing, false, Int64[])

initial(alg::DRparam) = alg.x⁰
tolerance(alg::DRparam) = alg.maxϵ
maxit(alg::DRparam) = alg.maxit
keephistory(alg::DRparam) = alg.keephistory
snapshots(alg::DRparam) = alg.snapshots



function solve(p::TwoSetsFP, alg::DRparam, x⁰, maxϵ, maxit, keephistory::Bool, snapshots::Vector{Int64})
    A = p.A
    B = p.B

    # process the default parameters
    !ismissing(x⁰) || ( x⁰ = getelement(A) )
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

        reflect!(yᵏ, xᵏ, B)
        reflect!(x̃ᵏ⁺¹, yᵏ, A)
        @. xᵏ⁺¹ = (x̃ᵏ⁺¹+ xᵏ) / 2 

        err .= xᵏ⁺¹ .- xᵏ # This doesn't say much in infeasible case, but is OK in case of binary aperture
        # dist .= xᵏ⁺¹ .- yᵏ # this calculates true error but can stay large in case of infeasible case
        ϵ = LinearAlgebra.norm(err)
        xᵏ .= xᵏ⁺¹
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

export DR, DRparam


