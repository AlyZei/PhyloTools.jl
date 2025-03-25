function sampling_ANC(p::Array{Float64,2}; n_seq = 100)
    L,q = size(p)
    res_msa = [sample(1:q, Weights(p[i,:]), n_seq) for i in 1:L]
    return Int.(hcat(res_msa...)')
end


function reshuffle(msa_in::Array{Int,2}, h::Array{T,2}, J::Array{T,4}; shuffle_units = 0.1, temp = 1., info = false) where {T}
    
    msa = deepcopy(msa_in)
    L, M = size(msa);
    unit = L*M*(M-1)/2
    N_shuffles = round(Int,shuffle_units * unit)
    
    acc_rate = 0
    ens = energy(msa_in, h, J)
    println(mean(ens))
    ens_evol = zeros(N_shuffles+1)
    ens_evol[1] = mean(ens)
    
    for iter in 1:N_shuffles
        i = rand(1:L)
        m, n = sample(1:M,2,replace = false)
        
        a = msa[i,m]
        b = msa[i,n]
        
        if a == b
            ens_evol[iter+1] = ens_evol[iter]
        else 
            dE_m = single_mut_dE(msa[:,m], h, J, b, i, L)
            dE_n = single_mut_dE(msa[:,n], h, J, a, i, L)
        
            if rand() < exp(-(dE_m + dE_n)/temp)
                msa[i,m] = b
                msa[i,n] = a
                acc_rate += 1
                ens_evol[iter+1] = ens_evol[iter] + ((dE_m + dE_n)/M)
            else
                ens_evol[iter+1] = ens_evol[iter]
            end
        end
    end
    
    println("Acceptance rate $(acc_rate) / $(N_shuffles) = $(acc_rate/N_shuffles)")
    
    if info == true
        return msa, ens_evol[1:1000:end], [i for i in 1:(N_shuffles+1)][1:1000:end] ./ unit
    else
        return msa
    end
end


function reshuffle_entr(msa_in::Array{Int,2}, h::Array{T,2}, J::Array{T,4}; times = 1., temp = 1., info = false) where {T}
    
    msa = deepcopy(msa_in)
    L, M = size(msa);
    unit = M*(M-1)/2
    
    axis_unit = L*M*(M-1)/2;
    
    
    f = compute_empirical_freqs(msa, 21, eps = 10^-5);
    #compute eff_number of aminos at every position, it's 2^CIE_i for evert site i
    eff_numb = round.(Int, times * unit .* (2 .^ sum(-f .* log.(f) ./ log(2), dims = 2)[:]));
    
    acc_rate = 0
    ens = energy(msa_in, h, J)
    println(mean(ens))
    #sample for the total effective number of aminoacids
    N_shuffles = Int.(sum(eff_numb))
    ens_evol = zeros(N_shuffles+1)
    ens_evol[1] = mean(ens)
    
    for iter in 1:round(Int,N_shuffles)
        #extract positions according to how many eff_number of aminos i have in that site
        i = sample(1:L, Weights(eff_numb))
        m, n = sample(1:M,2,replace = false)
        
        a = msa[i,m]
        b = msa[i,n]
        
        if a == b
            ens_evol[iter+1] = ens_evol[iter]
        else 
            dE_m = single_mut_dE(msa[:,m], h, J, b, i, L)
            dE_n = single_mut_dE(msa[:,n], h, J, a, i, L)
        
            if rand() < exp(-(dE_m + dE_n)/temp)
                msa[i,m] = b
                msa[i,n] = a
                acc_rate += 1
                ens_evol[iter+1] = ens_evol[iter] + ((dE_m + dE_n)/M)
            else
                ens_evol[iter+1] = ens_evol[iter]
            end
        end
    end
    
    println("Acceptance rate $(acc_rate) / $(N_shuffles) = $(acc_rate/N_shuffles)")
    
    if info == true
        return msa, ens_evol[1:1000:end], [i for i in 1:(N_shuffles+1)][1:1000:end] ./ axis_unit
    else
        return msa
    end
end


function reshuffle_entr_min(msa_in::Array{Int,2}, h::Array{T,2}, J::Array{T,4}; temp = 1, info = false) where {T}
    
    msa = deepcopy(msa_in)
    L, M = size(msa);
    unit = M*(M-1)/2
    
    axis_unit = L*M*(M-1)/2;
    
    
    f = compute_empirical_freqs(msa, 21, eps = 10^-5);
    eff_numb = round.(Int, unit .* (2 .^ sum(-f .* log.(f) ./ log(2), dims = 2)[:]));
    
    acc_rate = 0
    ens = energy(msa_in, h, J)
    println(mean(ens))
    N_shuffles = Int.(sum(eff_numb))
    ens_evol = zeros(N_shuffles+1)
    ens_evol[1] = mean(ens)
    
    for iter in 1:N_shuffles
        i = sample(1:L, Weights(eff_numb))
        m, n = sample(1:M,2,replace = false)
        
        a = msa[i,m]
        b = msa[i,n]
        
        if a == b
            ens_evol[iter+1] = ens_evol[iter]
        else 
            dE_m = single_mut_dE(msa[:,m], h, J, b, i, L)
            dE_n = single_mut_dE(msa[:,n], h, J, a, i, L)
        
            if rand() < exp(-(min(dE_m,dE_n))/temp)
                msa[i,m] = b
                msa[i,n] = a
                acc_rate += 1
                ens_evol[iter+1] = ens_evol[iter] + ((dE_m + dE_n)/M)
            else
                ens_evol[iter+1] = ens_evol[iter]
            end
        end
    end
    
    println("Acceptance rate $(acc_rate) / $(N_shuffles) = $(acc_rate/N_shuffles)")
    
    if info == true
        return msa, ens_evol[1:1000:end], [i for i in 1:(N_shuffles+1)][1:1000:end] ./ axis_unit
    else
        return msa
    end
end
    
            
            
        
        
    

    
    