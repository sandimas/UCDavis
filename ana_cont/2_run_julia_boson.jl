using FileIO
using JLD2
using PyCall
using Statistics



# define gaussian function for the model
function gauss(x, x0, σ)
    return @. 2.0*exp(-(x - x0)^2 / (2 * σ^2)) / (sqrt(2 * π) * σ)
end


# main function
function run_AC()

    # import python modules
    np = pyimport("numpy")
    cont = pyimport("ana_cont.continuation")
    solvers = pyimport("ana_cont.solvers")

    # Output directory
    mkpath("results")

    # Simulation parameters and output omegas
    ws = collect(LinRange(0.0, 20.0, 1001))
    taus = collect(LinRange(0.0, 20.0, 201))[1:101]
    beta = 20.0
    

    data_dict = load(joinpath("../data_files/AC_data.jld2"))

    #############################################
    ## 1. Density of States for the real-space G(0,0)
    #############################################

    # Load data
    G = data_dict["curr"][:,1:101]

    # Calculate the covariance matrix and average Green's function
    covar = Statistics.cov(G) # (taus,taus)
    G_avg = mean(G, dims=1)[1,:] # (taus)

    # Create the problem
    probl = cont.AnalyticContinuationProblem( im_axis = taus, re_axis = ws, im_data = G_avg, kernel_mode="time_bosonic", beta = beta)

    # Solve the problem
    sol, _ = probl.solve(method="maxent_svd", alpha_determination="chi2kink", cov=covar, model=gauss(ws, 0.0, 3.0),alpha_start=1e12,alpha_end=1e-3)
    
    # Save the results
    # sol.A_opt is the Density of States (shape: (601))
    save(joinpath("results", "julia_ana_cont_curr.jld2"), "sol", sol.A_opt, "ws", ws)



end



# Run the script if called from the command line
if abspath(PROGRAM_FILE) == @__FILE__
    run_AC()
end