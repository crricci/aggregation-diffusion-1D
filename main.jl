
include("L_loadAll.jl")


function solveModel(p)
    T_span = (0.0,p.T_end)
    cb = ManifoldProjection(comManifold)
    prob = ODEProblem(df!, p.y₀, T_span, p)
    sol = solve(prob,Rosenbrock23(),save_everystep=false,callback = cb)
            
    return sol
end

function comManifold(resid,u,p,t)
    resid[1] = sum(u .* p.x) * p.Δx - p.centerOfMass
    resid[2:end] .= 0.0
end 




# isoutofdomain=(u,p,t) -> sum(u .* p.x) > com + p.Δx || sum(u .* p.x) < com - p.Δx)