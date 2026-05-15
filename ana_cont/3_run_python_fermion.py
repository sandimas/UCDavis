import numpy as np
import ana_cont.continuation as cont
import ana_cont.solvers as solvers


def gauss(x, x0, sigma):
    return 2.0 * np.exp(-(x - x0)**2 / (2 * sigma**2)) / (np.sqrt(2 * np.pi) * sigma)

def run_AC():
    # import data
    # G(beta) = 1-G(0) so we need to limit to linearly independent time steps
    G_r = np.load("../data_files/G_r.npy")[:,0:199]
    G_k = np.load("../data_files/G_k.npy")[:,0:199,:]
    
    # params
    taus = np.linspace(0.0, 20.0, 201)[0:199]
    ws = np.linspace(-15.0, 15.0, 1001)
    beta = 20.0
    


    # DOS for the real-space Green's function
    print("Solving for the real-space Green's function")
    cov_r = np.cov(G_r.T)
    G_r_avg = np.mean(G_r, axis=0)
    model = gauss(ws, 0.0, 3.0)

    probl = cont.AnalyticContinuationProblem(im_axis=taus, re_axis=ws, im_data=G_r_avg, kernel_mode="time_fermionic", beta=beta)
    sol, _ = probl.solve(method="maxent_svd", alpha_determination="chi2kink", cov=cov_r, model=model, alpha_start=1e12, alpha_end=1e-3)
    A_r = sol.A_opt
    np.save("results/python_ana_cont_DOS.npy", A_r)

    # Spectral Function for the momentum-space Green's function
    # k = 1: Γ, k = 8: X, k = 15: M
    for k in range(21):
        print("Solving for G_k = {}".format(k))
        G_k_avg = np.mean(G_k[:,:,k], axis=0)
        cov_k = np.cov(G_k[:,:,k].T)
        probl = cont.AnalyticContinuationProblem(im_axis=taus, re_axis=ws, im_data=G_k_avg, kernel_mode="time_fermionic", beta=beta)
        sol, _ = probl.solve(method="maxent_svd", alpha_determination="chi2kink", cov=cov_k, model=gauss(ws, 0.0, 3.0), alpha_start=1e12, alpha_end=1e-3)
        A_k = sol.A_opt
        np.save("results/python_ana_cont_spec_k{}.npy".format(k), A_k)
    return

run_AC()