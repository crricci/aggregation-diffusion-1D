
function solve_l(;Vₐ=0.1, β=1.0, γ = 1.0, δ = 1.0, L=3.0,
        h=0.1,ϕ=1.0,ψ=1.0,ν = 0.0, V_GRD = 0.0, initialBumps = 1, T_fin = 10.0)
        
    Δx = 1e-2
    Nx = Int(L/Δx)
    x = LinRange(0,L,Nx)


    T_span = (0.0,T_fin)
    l₀ = initialCondition(x,L,Δx,initialBumps) 
    p = (Vₐ,β,γ,δ,L,h,ϕ,ψ,ν,V_GRD,x,Nx,Δx)
    prob = ODEProblem(df!,l₀,T_span,p)
    sol = solve(prob,Rosenbrock23())

    return sol,p
end

function df!(dl,l,p,t)
    (Vₐ,β,γ,δ,L,h,ϕ,ψ,ν,V_GRD,x,Nx,Δx) = p
    p1 = 0.2
    p2 = L-0.2   

    w = conv_w(l,p)  
    ∇w = ∇(w,p)
    ∇l = ∇(l,p)
    Δl = Δ(l.^δ,p)
    ∇V = ∇([V(xi,p1,p2) for xi in x],p)
    dl .= -Vₐ * ∇(l.^γ.*∇w,p) + ν*Δl + ∇(l .* ∇V,p) + V_GRD * abs.(∇l).^2 
end

function ∇(f,p)
    """gradient dirchlet boundary
    """
    (Vₐ,β,γ,δ,L,h,ϕ,ψ,ν,V_GRD,x,Nx,Δx) = p
    ∇f = similar(f)
    for i=2:Nx-1
        ∇f[i] = (f[i+1] - f[i-1])/(2*Δx)
    end
    ∇f[1] = (f[2] - 0.0)/(2*Δx)
    ∇f[Nx] = (0.0 - f[Nx-1])/(2*Δx)
    return ∇f
end

function Δ(f,p)
    """laplacian dirichlet boundary 
    """
    (Vₐ,β,γ,δ,L,h,ϕ,ψ,ν,V_GRD,x,Nx,Δx) = p
    Δf = similar(f)
    for i=2:Nx-1
        Δf[i] = (f[i+1] - 2*f[i] + f[i-1])/Δx^2
    end
    Δf[1] = (f[2] - 2*f[1] + 0.0)/Δx^2
    Δf[Nx] = (0.0 - 2*f[Nx] + f[Nx-1])/Δx^2
    return Δf
end


function conv_w(l,p)
    """ lϕ = l^ϕ 
        Whlϕ = (Wh*l^ϕ) 
        Whlϕψ = (Wh*l^ϕ)^ψ
        w = (Wh*l^ϕ)^ψ l^(β-1) 
    """ 
    (Vₐ,β,γ,δ,L,h,ϕ,ψ,ν,V_GRD,x,Nx,Δx) = p
    Nih = ceil(Int,h/Δx)

    lϕ = abs.(l).^ϕ
    Whlϕ = zeros(eltype(l),size(l))
    for i in 1:Nx
        j_range = mod.(i-Nih-2:i+Nih+1,Nx).+1 
        for k in 1:length(j_range)  # TRAPEZ
        jZₙₓ = j_range[k]
        jZₙₓ_next = k < length(j_range) ? j_range[k+1] : j_range[1]
        Whlϕ[i] += Δx*(W(d(x[i],x[jZₙₓ],L),h)*lϕ[jZₙₓ] + 
                W(d(x[i],x[jZₙₓ_next],L),h)*lϕ[jZₙₓ_next])/2
        end
    end

    Whlϕψ = abs.(Whlϕ).^ψ 
    w = Whlϕψ .* abs.(l).^(β-1)

    return w
end
