

function df!(du,u,p,t)
    w = conv_w(u,p)
    ∇w = ∇(w,p)
    ∇u = ∇(u,p)
    Δu = Δ(u.^(p.δ),p)
    du .= - p.γa * ∇(u.^p.γ.*∇w,p) + p.γd * Δu + ∇(u.*p.∇V,p) + p.V_GRD * abs.(∇u).^2 
end

function ∇(f,p)
    """ gradient dirchlet boundary 0.0 """
    @unpack_parameters p
    ∇f = similar(f)
    for i=2:Nx-1
        ∇f[i] = (f[i+1] - f[i-1])/(2*Δx)
    end
    ∇f[1] = (f[2] - 0.0)/(2*Δx)
    ∇f[Nx] = (0.0 - f[Nx-1])/(2*Δx)
    return ∇f
end

function Δ(f,p)
    """ laplacian dirichlet boundary 0.0 """
    @unpack_parameters p
    Δf = similar(f)
    for i=2:Nx-1
        Δf[i] = (f[i+1] - 2*f[i] + f[i-1])/Δx^2
    end
    Δf[1] = (f[2] - 2*f[1] + 0.0)/Δx^2
    Δf[Nx] = (0.0 - 2*f[Nx] + f[Nx-1])/Δx^2
    return Δf
end

function conv_w(f,p)
    """ fϕ = f^ϕ 
        Whfϕ = (Wh*f^ϕ) 
        Whfϕψ = (Wh*f^ϕ)^ψ
        w = (Wh*f^ϕ)^ψ f^(β-1) 
    """ 
    @unpack_parameters p
    Nih = ceil(Int,h/Δx)

    fϕ = abs.(f).^ϕ
    Whfϕ = zeros(eltype(f),size(f))
    for i in 1:Nx
        j_range = unique( max.(1, min.(Nx, i-Nih-1:i+Nih+1)) )
        for k in 1:length(j_range)  # TRAPEZ
        jZₙₓ = j_range[k]
        jZₙₓ_next = k < length(j_range) ? j_range[k+1] : j_range[1]
        Whfϕ[i] += Δx*(W(abs(x[i]-x[jZₙₓ]),h)*fϕ[jZₙₓ] + 
                W(abs(x[i]-x[jZₙₓ_next]),h)*fϕ[jZₙₓ_next])/2
        end
    end

    Whfϕψ = abs.(Whfϕ).^ψ 
    w = Whfϕψ .* abs.(f).^(β-1)

    return w
end
