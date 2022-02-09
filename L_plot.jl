
function plotSol(sol,t_step,p)
    (Vₐ,β,γ,δ,L,h,ϕ,ψ,ν,V_GRD,x,Nx,Δx) = p

    T_end = sol.t[end] 
    lmin = minimum(sol)
    lmax = maximum(sol)

    fig = figure()
    show()
    for t in 0.0:t_step:T_end
        plot(x,sol(t))
        ylim(lmin,lmax)
        title("t = $t")
        sleep(0.1)
        fig.clear()
    end
    close()
end

function plotSlice(sol,times, p)
    (Vₐ,β,γ,δ,L,h,ϕ,ψ,ν,V_GRD,x,Nx,Δx) = p
    for t in times
        plot(x,sol(t),label = "t = $t")
        xlabel("Space")
        ylabel("Income")
        grid(true)
        legend()
    end
end

function plotSlice()
    
    # two bumps initial
    sol1,p = solve_l(ν = 0.04, h = 0.2, initialBumps = 2, T_fin = 100.0);
    sol2,p = solve_l(ν = 0.04, h = 0.3, initialBumps = 2, T_fin = 100.0);

    (Vₐ,β,γ,δ,L,h,ϕ,ψ,ν,V_GRD,x,Nx,Δx) = p

    times1 = [0.0,1.0,2.0,100.0]
    for t in times1
        subplot(1,2,1)
        plot(x,sol1(t),label = "t = $t")
        xlabel("Space")
        ylabel("Income")
        grid(true)
        legend()
    end

    times2 = [0.0,1.0,2.0,100.0]
    for t in times2
        subplot(1,2,2)
        plot(x,sol2(t),label = "t = $t")
        xlabel("Space")
        ylabel("Income")
        grid(true)
        legend()
    end
    tight_layout()

    # flat initial
    sol3,p = solve_l(ν = 0.04, h = 0.2, initialBumps = 1, T_fin = 100.0);
    sol4,p = solve_l(ν = 0.04, h = 0.3, initialBumps = 1, T_fin = 100.0);

    times3 = [0.0,40.0,50.0,100.0]
    for t in times3
        subplot(1,2,1)
        plot(x,sol3(t),label = "t = $t")
        xlabel("Space")
        ylabel("Income")
        grid(true)
        legend()
    end

    times4 = [0.0,60.0,70.0,100.0]
    for t in times4
        subplot(1,2,2)
        plot(x,sol4(t),label = "t = $t")
        xlabel("Space")
        ylabel("Income")
        grid(true)
        legend()
    end
    tight_layout()

    
end