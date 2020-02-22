using Statistics
using Plots
using StatsBase

ppal_dir = pwd()

function init_lattice(N, M, prob)
   
    lattice = zeros((N,M))
    
    for i in 1:N
        for j in 1:M
            
            dice = rand()
            
            if dice < prob
                
                lattice[i,j] = 1
                
            else
                lattice[i,j] = -1
                
            end
            
        end
    end
        
    return lattice
    
end

function energy(i, j, lattice, N, M)
    
    u = 1 + mod(j, M)   # up
    d = 1 + mod(j-2,M) # down
    r = 1 + mod(i, N)   # right
    l = 1 + mod(i-2, N) # left
    
    Sn = lattice[i,u] + lattice[i,d] + lattice[l,j] + lattice[r,j]
    
    E = -Sn * lattice[i,j]
    
    return E
    
end

function magnetisation(lattice, N, M)
    
    return abs(sum(lattice))# / (N * M)
    
end

function total_energy(lattice, N, M)
    
    E = 0.
    
    for i in 1 : N
        for j in 1 : M
            
            E += energy(i, j, lattice, N, M)
            
        end
    end
            
    return E / 4
    
end

function metropolis(lattice, N, M, beta, i, j)
    
    E_0 = energy(i, j, lattice, N, M)

    lattice[i,j] *= -1

    E_f = energy(i, j, lattice, N, M)

    delta_E = E_f - E_0

    prob = exp(-beta * delta_E)

    #Just keep the flipped value if the energy decreases, else keep only with probability e^(-\beta(delta_E))
    if !(delta_E < 0 || (rand() < prob) )
        
        lattice[i,j] *= -1
    end

    return lattice
    
end

function simulation(N, M, T, t, eq_time, K, lattice)
   
    magnetization_time = zeros(t)
    energy_time = zeros(t)
    
    magnetization_squared_time = zeros(t)
    energy_squared_time = zeros(t)
    
    beta = 1. / T
    N_dipols = N * M
    
    #Thermalisation
    for MCS in 1 : eq_time

        for i in 1 : N
            for j in 1 : M
            
                lattice = metropolis(lattice, N, M, beta, i, j)
                
            end
        end
            
    end
    
    #Calculation
    for MCS in 1 : t
        
        for k in 1 : K

            for i in 1 : N
                for j in 1 : M

                    lattice = metropolis(lattice, N, M, beta, i, j)

                end
            end
            
        end
        
        mag = magnetisation(lattice, N, M)
        ene = total_energy(lattice, N, M)

        magnetization_time[MCS] = mag
        energy_time[MCS] = ene
        
        magnetization_squared_time[MCS] = mag^2
        energy_squared_time[MCS] = ene^2

    end
    
    return magnetization_time, energy_time, magnetization_squared_time, energy_squared_time
        
end

function exp_corr(array)
    
    length = size(array)[1]
    corr = 0.
    
    @fastmath @inbounds for i in 1 : length - 1
        
        corr += array[i] * array[i + 1]
        
    end
        
    x = corr / length
    
    rho = (x - mean(array)^2) / var(array)
    
    final_corr = rho / (1 - rho)
    
    return final_corr
    
end

function avg_error(array)
    
    avg = mean(array)
    
    tau = exp_corr(array)
    
    L = length(array)
    
    if tau > 0
    
        error = sqrt(var(array) / L) * sqrt(2 * tau + 1)
        
    else
        
        error = sqrt(var(array) / L)
        
    end
    
    return (avg, error)
    
end

function simulate(N, M, T, prob_up, t, eq_time, K)
   
    lattice = init_lattice(N, M, prob_up)
    
    M_arr, E_arr, M_squared_arr, E_squared_arr = simulation(N, M, T, t, eq_time, K, lattice)
    
    #<M> and <E>
    avg_M = avg_error(M_arr)
    avg_E = avg_error(E_arr)
    
    #Susceptibility and Specific heat
    
    avg_M_squared = avg_error(M_squared_arr)
    avg_E_squared = avg_error(E_squared_arr)
    
    beta = 1. / T
    
    #Errors computed from error propagation
    err_xi = beta * sqrt(avg_M_squared[2]^2 + 4 * avg_M[1]^2 * avg_M[2]^2)
    err_Cv = beta^2 * sqrt(avg_E_squared[2]^2 + 4 * avg_E[1]^2 * avg_E[2]^2)
    
    xi = (beta * (avg_M_squared[1] - avg_M[1]^2), err_xi)
    Cv = (beta^2 * (avg_E_squared[1] - avg_E[1]^2), err_Cv)
    
    return avg_M, avg_E, xi, Cv
    
end

function simulation_T(N, M, eq_time, t, K, Ts, folder="Data")
    
    M_T = zeros(length(Ts))
    E_T = zeros(length(Ts))

    xi_T = zeros(length(Ts))
    Cv_T = zeros(length(Ts))

    M_T_errors = zeros(length(Ts))
    E_T_errors = zeros(length(Ts))
    
    xi_T_errors = zeros(length(Ts))
    Cv_T_errors = zeros(length(Ts))
    
    N_spins = N * M #Compute previously to gain in efficiency

    pos = 0
    
    if ! isdir(folder)
        mkdir(folder)
    end
    
    cd(folder)
    
    f_results = open("results.txt", "w")
    f_errors = open("errors.txt", "w")
    
    println(f_results, "#<M>\t<E>\txi\tCv\tT")
    println(f_errors, "#<M>\t<E>\txi\tCv")

    for T in Ts
        
        println("T: ", T)

        pos += 1

        if T < 2

            prob_up = 1

        else

            prob_up = 0.5

        end

        avg_M, avg_E, xi, Cv = simulate(N, M, T, prob_up, t, eq_time, K)

        M_T[pos] = avg_M[1] / N_spins
        E_T[pos] = avg_E[1] / N_spins

        xi_T[pos] = xi[1] / N_spins
        Cv_T[pos] = Cv[1] / N_spins

        M_T_errors[pos] = avg_M[2] / N_spins
        E_T_errors[pos] = avg_E[2] / N_spins
        
        xi_T_errors[pos] = xi[2] / N_spins
        Cv_T_errors[pos] = Cv[2] / N_spins

        println(f_results, M_T[pos] , "\t", E_T[pos], "\t", xi_T[pos], "\t",Cv_T[pos], "\t", T)
        println(f_errors, M_T_errors[pos], "\t", E_T_errors[pos], "\t", xi_T_errors[pos], "\t", Cv_T_errors[pos])

    end

    close(f_results)
    close(f_errors)
    
    cd(ppal_dir)
    
    return M_T, E_T, xi_T, Cv_T, M_T_errors, E_T_errors, xi_T_errors, Cv_T_errors
    
end

function check_corr_t(N, t)
    
    M = N

    T = 2.22

    prob_up = 0.5

    eq_time = 10^4

    lattice = init_lattice(N, M, prob_up)

    M_arr, E_arr, M_squared_arr, E_squared_arr = simulation(N, M, T, t, eq_time, lattice)
    
    lags = [i for i in 1 : t - 1]

    corr = autocor(M_arr, lags)

    return corr
    
end

function correlation_time(x)
    
    tau = 0.
    
    for val in x
        
        if val > 0
            
            tau += val
        
        else
            
            return tau
            break
            
        end
    end
    
    return tau
    
end

function corr_t_study(times, Ls, filename)

    f = open("$filename.txt", "w")

    println("#L\tcorrelation time")
    println(f, "#L\tcorrelation time")

    t = 10^5

    for L in Ls

        avg_tau = 0.

        for i in 1 : times

            corr = check_corr_t(L, t)

            avg_tau += correlation_time(corr)

        end

        avg_tau = avg_tau / times

        println(L, "\t", avg_tau)
        println(f, L, "\t", avg_tau)

    end

    close(f)
    
end

function simulation_K(N, M, T, t, eq_time, lattice, tau_M)
   
    magnetization_time = zeros(t)
    energy_time = zeros(t)
    
    magnetization_squared_time = zeros(t)
    energy_squared_time = zeros(t)
    
    beta = 1. / T
    N_dipols = N * M
    
    #Thermalisation
    for MCS in 1 : eq_time

        for i in 1 : N
            for j in 1 : M
            
                lattice = metropolis(lattice, N, M, beta, i, j) 
                
            end
        end
            
    end
    
    time_MC = 0.
    time_measure = 0.
    
    #Calculation
    for MCS in 1 : t
        
        time_start = time_ns()

        for i in 1 : N
            for j in 1 : M
            
                lattice = metropolis(lattice, N, M, beta, i, j)
                
            end
        end
        
        time_MC = time_ns() - time_start
        
        time_start = time_ns()
        
        mag = magnetisation(lattice, N, M)
        ene = total_energy(lattice, N, M)

        magnetization_time[MCS] = mag
        energy_time[MCS] = ene
        
        magnetization_squared_time[MCS] = mag^2
        energy_squared_time[MCS] = ene^2
        
        time_measure = time_ns() - time_start
        
        opt = sqrt(2*tau_M*time_measure/time_MC)
        
        println("Time MC: ", time_MC, "\t", "Time measure: ", time_measure, "\t", "Opt: ", opt)

    end
    
    return magnetization_time, energy_time, magnetization_squared_time, energy_squared_time
        
end

#Parameters
N = 64
M = N

eq_time = 10^4
t = 10^6
K = 5

Ts = range(2.0, stop=3.0, length=50)

folder = "L_$N"

#Start simulation

t0 = time_ns()

M_T, E_T, xi_T, Cv_T, M_T_errors, E_T_errors, xi_T_errors, Cv_T_errors = simulation_T(N, M, eq_time, t, K, Ts, folder)

tf = time_ns()

#Simulation ended

#Write parameters used and elapsed time to file
    
ET = round((tf - t0) * 1e-9, sigdigits=3)
    
f_parameters_used = open("$folder/parameters_used.txt", "w")

println(f_parameters_used, "N:", N)
println(f_parameters_used, "M:", M)
println(f_parameters_used, "Run_MCS:", t)
println(f_parameters_used, "Eq_MCS:", eq_time)
println(f_parameters_used, "T0:", Ts[1])
println(f_parameters_used, "Tf:", Ts[end])
println(f_parameters_used, "dT:", abs(Ts[2]-Ts[1]))

println(f_parameters_used, "\nElapsed time:", ET)

close(f_parameters_used)

println("Elapsed time: ", ET, " s")
