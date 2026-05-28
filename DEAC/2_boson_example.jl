using SmoQyDEAC
using JLD2
using Statistics

function run_example()


    # Run this multithreaded!
    if Threads.nthreads() == 1
        println("\nRunning on 1 thread. It is recommended to run on multiple threads for faster computation.")
        println("\t julia --threads=auto 2_boson_example.jl\n")
    end

    # Output directory
    mkpath("results")

    # Load the data
    data_dict = load("../data_files/AC_data.jld2")
    curr = data_dict["curr"]
    
    # set up parameters
    taus = collect(LinRange(0.0, 20.0, 201)) # no need to truncate for linear independence, SmoQyDEAC will handle this
    ws = collect(LinRange(0.0, 20.0, 1001))
    beta = 20.0
    kernel_mode = "time_bosonic_symmetric_w"
    n_bins = 100
    n_runs_per_bin = 100


    # DEAC on G(0,0) for DOS
    # using covariance matrix
    deac_dict = DEAC_Binned(
        curr[:,1:101], # shape (bins,taus)
        beta,
        taus[1:101],
        ws,
        kernel_mode, # "time_fermionic" or "time_bosonic"
        n_bins, # number of bins to split the data into
        n_runs_per_bin, # number of runs to perform per bin
        "results/DEAC_curr_cov.jld2", # output file
        "chk.jld2", # checkpoint file
        find_fitness_floor = false,
        number_of_generations = 20_000, #max number of generations within a run
        keep_bin_data = false, # keep the bin data and store it in the output file
        verbose = true, # print progress
        fitness = [1000.0, 100.0, 10.0, 2.0, 1.3, 1.2,1.0], # target χ² fitnesses for the DEAC algorithm
                                                        # Typically, you would only want to use 1.0 unless that overfits.
                                                        # I am doing multiple fitnesses to show you how DEAC converges
    )

    # using standard deviation
    deac_dict = DEAC_Std(
        mean(curr,dims=1)[1,1:101], # shape (taus)
        std(curr,dims=1)[1,1:101], # shape (taus)
        beta,
        taus[1:101],
        ws,
        kernel_mode, # "time_fermionic" or "time_bosonic"
        n_bins, # number of bins to split the data into
        n_runs_per_bin, # number of runs to perform per bin
        "results/DEAC_curr_std.jld2", # output file
        "chk.jld2", # checkpoint file
        find_fitness_floor = false,
        number_of_generations = 20_000, #max number of generations within a run
        keep_bin_data = false, # keep the bin data and store it in the output file
        verbose = true, # print progress
        fitness = [ 1.0, 1e-1, 1e-2, 1e-3, 5e-4, 1e-4], # target χ² fitnesses for the DEAC algorithm
                                              # Fitness values are far harder to pin down using standard deviation.
                                              # I am doing multiple fitnesses to show you how DEAC converges
    )

end


run_example()