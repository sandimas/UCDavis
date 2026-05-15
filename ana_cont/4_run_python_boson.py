import numpy as np
import ana_cont.continuation as cont
import ana_cont.solvers as solvers


def gauss(x, x0, sigma):
    return 2.0 * np.exp(-(x - x0)**2 / (2 * sigma**2)) / (np.sqrt(2 * np.pi) * sigma)

def run_AC():
    # import data
    # Curr(tau) is only linearly independent in [0, beta/2] so we only use the first 101 time steps
    curr = np.load("../data_files/curr.npy")[:,0:100]
    
    # params
    taus = np.linspace(0.0, 20.0, 201)[0:100]
    ws = np.linspace(0.0, 20.0, 1001)
    beta = 20.0
    


    # DOS for the real-space Green's function
    print("Solving for the real-space Green's function")
    cov = np.cov(curr.T)
    curr_avg = np.mean(curr, axis=0)

    probl = cont.AnalyticContinuationProblem(im_axis=taus, re_axis=ws, im_data=curr_avg, kernel_mode="time_bosonic", beta=beta)
    sol, _ = probl.solve(method="maxent_svd", alpha_determination="chi2kink", cov=cov, model=gauss(ws, 0.0, 3.0), alpha_start=1e12, alpha_end=1e-3)
    A = sol.A_opt
    np.save("results/python_ana_cont_curr.npy", A)
    return

run_AC()