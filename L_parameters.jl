

function initialCondition(x,L,Δx,initialBumps)

    if initialBumps == 1
        σ₀ = 2.0
        L₀ = Normal(L/2,σ₀)
        l₀ = pdf.(L₀,x) .+ 1.0
    elseif initialBumps == 2
        σ₀ = 0.1
        L₀¹ = Normal(L/4+0.1,σ₀)
        L₀² = Normal(3/4*L-0.1,σ₀)
        l₀ = pdf.(L₀¹,x) + pdf.(L₀²,x) .+ 1.0
    else
        R = 0.2
        l₀ = Float64.([ind(xi-L/2,R,R) for xi in x])
    end

    return l₀ / (sum(l₀) * Δx)
end


function W(x)
    return abs(x) <= 1 ? 1-abs(x) : 0.0
end

function W(x,h)
    return 1/h*W(x/h)
end

function V(x,p1,p2)
    """ confining potential between p1 and p2"""
    r = 0.1
    p1 = p1 + r
    p2 = p2 - r
    if p1 < x < p2 
        return 0.0
    elseif x ≤ p1
        return V₀(p1-x,r)
    elseif x ≥ p2
        return V₀(x-p2,r)
    else
        error("somethings wrong")
    end
end

function V₀(x,r)
    k = 100
    return 1/(1 + exp(-k*(x-r)))
end



function g_smooth(x,ϵ)
    return  abs(x) <= ϵ ? exp(1/(x^2-ϵ^2))/exp(-1/ϵ^2) : 0
end

function ind(x,R,ϵ)
    if -R+ϵ <= x <= R-ϵ
        return 1
    elseif x > R+ϵ || x < -R-ϵ  
        return 0
    elseif R-ϵ <= x <= R+ϵ    
        return g_smooth(x-(R-ϵ),2*ϵ)
    else 
        return g_smooth(-x-(R-ϵ),2*ϵ)
    end
end


function d(x,y,L)
    return min(abs(x-y),L-abs(x-y))
end