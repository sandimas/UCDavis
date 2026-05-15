using SmoQyDEAC
using JLD2
using Statistics

function run_example()


    # Run this multithreaded!
    if Threads.nthreads() == 1
        println("\nRunning on 1 thread. It is recommended to run on multiple threads for faster computation.")
        println("\t julia --threads=auto 1_fermion_example.jl\n")
    end

    # Output directory
    mkpath("results")

    # Load the data
    data_dict = load("../data_files/AC_data.jld2")
    G_r = data_dict["G_r"]
    G_k = data_dict["G_k"]
    
    # set up parameters
    taus = collect(LinRange(0.0, 20.0, 201)) # no need to truncate for linear independence, SmoQyDEAC will handle this
    ws = collect(LinRange(-15.0, 15.0, 1001))
    beta = 20.0
    kernel_mode = "time_fermionic"
    n_bins = 100
    n_runs_per_bin = 100


    # # DEAC on G(0,0) for DOS
    # # using covariance matrix
    # deac_dict = DEAC_Binned(
    #     G_r, # shape (bins,taus)
    #     beta,
    #     taus,
    #     ws,
    #     kernel_mode, # "time_fermionic" or "time_bosonic"
    #     n_bins, # number of bins to split the data into
    #     n_runs_per_bin, # number of runs to perform per bin
    #     "results/DEAC_DOS_cov.jld2", # output file
    #     "chk.jld2", # checkpoint file
    #     find_fitness_floor = false,
    #     number_of_generations = 20_000, #max number of generations within a run
    #     keep_bin_data = false, # keep the bin data and store it in the output file
    #     verbose = true, # print progress
    #     fitness = [1000.0, 100.0, 10.0, 1.0]#, 0.5], # target χ² fitnesses for the DEAC algorithm
    #                                                     # Typically, you would only want to use 1.0 unless that overfits.
    #                                                     # I am doing multiple fitnesses to show you how DEAC converges
    # )

    # # using standard deviation
    # deac_dict = DEAC_Std(
    #     mean(G_r,dims=1)[1,:], # shape (taus)
    #     std(G_r,dims=1)[1,:], # shape (taus)
    #     beta,
    #     taus,
    #     ws,
    #     kernel_mode, # "time_fermionic" or "time_bosonic"
    #     n_bins, # number of bins to split the data into
    #     n_runs_per_bin, # number of runs to perform per bin
    #     "results/DEAC_DOS_std.jld2", # output file
    #     "chk.jld2", # checkpoint file
    #     find_fitness_floor = false,
    #     number_of_generations = 20_000, #max number of generations within a run
    #     keep_bin_data = false, # keep the bin data and store it in the output file
    #     verbose = true, # print progress
    #     fitness =  [ 100.0, 10.0, 1.0, 0.1, 0.01, 0.001,  0.0005]#0.0001], # target χ² fitnesses for the DEAC algorithm
    #                                                     # Typically, you would only want to use 1.0 unless that overfits.
    #                                                     # I am doing multiple fitnesses to show you how DEAC converges
    # )

    # DEAC on G(k,0) for spectral function
    for k in 15:16
        deac_dict = DEAC_Binned(
            G_k[:,:,k], # shape (bins,taus)
            beta,
            taus,
            ws,
            kernel_mode, # "time_fermionic" or "time_bosonic"
            n_bins, # number of bins to split the data into
            n_runs_per_bin, # number of runs to perform per bin
            "results/DEAC_spec_$(k).jld2", # output file
            "chk.jld2", # checkpoint file
            find_fitness_floor = false,
            number_of_generations = 5_000, #max number of generations within a run
            keep_bin_data = false, # keep the bin data and store it in the output file
            verbose = true, # print progress
            fitness = 1.5, # target χ² fitnesses for the DEAC algorithm 
        )
    end

end


run_example()