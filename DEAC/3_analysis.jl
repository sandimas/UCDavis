using CairoMakie
using JLD2

function analysis()
    ########################################################
    # DOS analysis
    ########################################################
    # Load the data
    data_cov = load("results/DEAC_DOS_cov.jld2")
    data_std = load("results/DEAC_DOS_std.jld2")
    data_ana_cont = load("../ana_cont/results/julia_ana_cont_DOS.jld2")

    f = Figure(size=(2000,600),fontsize=40)
    ax = Axis(f[1,1],xlabel="ω",ylabel="A(ω)",title="DEAC (covariance)")
    for (fi, fitness) in enumerate(data_cov["fitness"])
        lines!(ax, data_cov["ωs"], data_cov["A"][:,fi] .+ ((fi - 1)*0.1), label="$fitness")
    end
    axislegend(ax)
    xlims!(ax,-8,12)

    ax = Axis(f[1,2],xlabel="ω",ylabel="A(ω)",title="DEAC (standard deviation)")
    for (fi, fitness) in enumerate(data_std["fitness"])
        lines!(ax, data_std["ωs"], data_std["A"][:,fi] .+ ((fi - 1)*0.1), label="$fitness")
    end
    axislegend(ax)
    xlims!(ax,-8,12)

    ax = Axis(f[1,3],xlabel="ω",ylabel="A(ω)",title="ana_cont")
    lines!(ax, data_ana_cont["ws"], data_ana_cont["sol"], label="AnaCont")
    axislegend(ax)
    xlims!(ax,-8,8)

    save("DOS.png", f)



    ########################################################
    # Current analysis
    ########################################################
    data_cov = load("results/DEAC_curr_cov.jld2")
    data_std = load("results/DEAC_curr_std.jld2")
    data_ana_cont = load("../ana_cont/results/julia_ana_cont_curr.jld2")


    f = Figure(size=(2000,600),fontsize=40)
    ax = Axis(f[1,1],xlabel="ω",ylabel="A(ω)",title="DEAC (covariance)")
    for (fi, fitness) in enumerate(data_cov["fitness"])
        lines!(ax, data_cov["ωs"], data_cov["A"][:,fi] .+ ((fi - 1)*0.1), label="$fitness")
    end
    axislegend(ax)
    xlims!(ax,0,8)

    ax = Axis(f[1,2],xlabel="ω",ylabel="A(ω)",title="DEAC (standard deviation)")
    for (fi, fitness) in enumerate(data_std["fitness"])
        lines!(ax, data_std["ωs"], data_std["A"][:,fi] .+ ((fi - 1)*0.1), label="$fitness")
    end
    axislegend(ax)
    xlims!(ax,0,8)

    ax = Axis(f[1,3],xlabel="ω",ylabel="A(ω)",title="ana_cont")
    lines!(ax, data_ana_cont["ws"], data_ana_cont["sol"], label="AnaCont")
    axislegend(ax)
    xlims!(ax,0,8)

    save("curr.png", f)

    ########################################################
    # Spectral function analysis
    ########################################################

    x_ticks = ([1,8,15,15+7*sqrt(2)],["Γ","X","M","Γ"])
	
	x_pos = collect(LinRange(1.0,22.0,22))
	for i in 1:7
	    x_pos[15+i] = 15.0 + (i*sqrt(2))
	end
    
    
    data_cov = load("results/DEAC_DOS_cov.jld2") # ["A"][:,4]
    
    data_ana_cont_dos = load("../ana_cont/results/julia_ana_cont_DOS.jld2")["sol"]
    data_ana_cont_ws = load("../ana_cont/results/julia_ana_cont_DOS.jld2")["ws"]
    data_deac = zeros(22,1001)
    data_ana_cont = zeros(22,1001)
    
    for k in 1:21
        data_deac[k,:] = load("results/DEAC_spec_$(k).jld2")["A"][:,end]
        data_ana_cont[k,:] = load("../ana_cont/results/julia_ana_cont_spec_k$(k).jld2")["sol"]
    end
    data_deac[22,:] = data_deac[1,:]
    data_ana_cont[22,:] = data_ana_cont[1,:]
    maxval = max(maximum(data_deac),maximum(data_ana_cont))

    f = Figure(size=(2000,600),fontsize=40)
    ga = GridLayout(f[1,1])

    ax1 = Axis(ga[1,1],ylabel="A(ω)",title="DEAC (covariance)")
    heatmap!(ax1, x_pos, data_cov["ωs"], data_deac,colormap=:terrain,colorrange=(0,maxval))#,colorscale=sqrt)
    ax2 = Axis(ga[1,2],)
    lines!(ax2, data_cov["A"][:,4], data_cov["ωs"], label="DEAC")
    xlims!(ax2,0,1.1*maximum(data_cov["A"][:,4]))
    ylims!(ax1,-8,9)
    ylims!(ax2,-8,9)
    colgap!(ga, 0)
    ax1.xticks = x_ticks
    ax2.xticklabelsvisible = false
    
    
    ax3 = Axis(ga[1,3],ylabel="A(ω)",title="ana_cont")
    heatmap!(ax3, x_pos, data_ana_cont_ws, data_ana_cont,colormap=:terrain,colorrange=(0,maxval))#,colorscale=sqrt)
    ax4 = Axis(ga[1,4])
    lines!(ax4, data_ana_cont_dos,data_ana_cont_ws, label="ana_cont")
    xlims!(ax4,0,1.1*maximum(data_cov["A"][:,4]))
    ylims!(ax3,-8,9)
    ylims!(ax4,-8,9)
    colgap!(ga,0)
    ax3.xticks = x_ticks
    ax4.xticklabelsvisible = false

    hideydecorations!(ax2)
    hideydecorations!(ax4)

    colsize!(ga,1, Relative(0.4))
    colsize!(ga,3, Relative(0.4))

    save("spec.png", f)





    return
    
end

analysis()