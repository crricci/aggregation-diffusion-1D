
@with_kw struct parameters
    # aggregation diffusion
    γa::Float64 = 0.1
    γd::Float64 = 0.04
    h::Float64 = 0.3
    
    # reallocation gain
    V_GRD::Float64 = 0.0
    V_GRA::Float64 = 0.0

    # domain
    L::Float64 = 4.0
    T_end::Float64 = 10.0
    b_L::Float64 = - L + 0.2
    b_R::Float64 = L - 0.2

    # numerical
    Δx::Float64 = 1e-2
    Nx::Int = Int(L/Δx)
    x::LinRange{Float64, Int64} = LinRange(-L,L,Nx)
    ∇V::Vector{Float64} = [ForwardDiff.derivative(y -> V(y,b_L,b_R), xi) for xi in x]

    # initial condition
    NBumps::Int = 1
    bumpCenter::Vector{Float64} = [0.0]
    bumpRadius::Vector{Float64} = [L/4]
    y₀::Vector{Float64} = initialCondition(NBumps,bumpCenter,bumpRadius,Δx,Nx,x)
    centerOfMass::Float64 = sum(y₀ .* x) * Δx


    # human capital 
    β::Float64 = 1.0
    γ::Float64 = 1.0
    δ::Float64 = 1.0
    ϕ::Float64 = 1.0
    ψ::Float64 = 1.0

    # checks
    @assert length(bumpCenter) == NBumps
    @assert length(bumpRadius) == NBumps
end


function initialCondition(NBumps,bumpCenter,bumpRadius,Δx,Nx,x; ϵ = 0.1)
    """ initialize initial condition for PDE by summing all the bumps """
    y₀ = zeros(Nx)
    for i in eachindex(x)
        xi = x[i]
        for j in 1:NBumps
            y₀[i] += bump(bumpCenter[j],x[i],bumpRadius[j],ϵ)
        end
    end 
    # return y₀
    return y₀/( sum(y₀) * Δx )
end


function W(x)
    """ smoothing kernel """
    return abs(x) <= 1 ? 1-abs(x) : 0.0
end

function W(x,h)
    """ rescaled kernel """
    return 1/h*W(x/h)
end

function V(x,p1,p2; r = 0.1)
    """ confining potential (logistic) between p1 and p2 
        midpoint at p1+r and p2-r 
    """
    if p1 + r/2 < x < p2 - r/2
        return 0.0
    elseif x < p1 - r/2 || x > p2 + r/2
        return 1 
    elseif p1 - r/2 ≤ x ≤ p1 + r/2
        return 1 - joint_smooth(-(x-(p1+r/2)),r)
    elseif p2 - r/2 ≤ x ≤ p2 + r/2
        return 1 - joint_smooth(x-(p2-r/2),r)
    end
end

function bump(center,x,R,ϵ)
    """ smooth bump centered at 'center' with radius R + ϵ/2 """
    return bump(x-center,R,ϵ)
end

function bump(x,R,ϵ)
    """ smooth bump centered at 0 with radius R + ϵ/2 """
    if -R+ϵ/2 < x < R-ϵ/2
        return 1
    elseif x > R+ϵ/2 || x < -R-ϵ/2  
        return 0
    elseif R-ϵ/2 ≤ x ≤ R+ϵ/2    
        return joint_smooth(x-(R-ϵ/2),ϵ)
    elseif -R-ϵ/2 ≤ x ≤ -R+ϵ/2    
        return joint_smooth(-(x-(-R+ϵ/2)),ϵ)
    end
end

function joint_smooth(x,ϵ)
    """ smooth joining from (0,1) to (ϵ,0) at point x """
    d = x/ϵ
    if d ≤ 1
        return 1 - 6 * d^2 + 8 * d^3 - 3 * d^4
    else 
        return 0
    end
end