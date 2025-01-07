function sampling_ANC(p::Array{Float64,2}; n_seq = 100)
    L,q = size(p)
    res_msa = [sample(1:q, Weights(p[i,:]), n_seq) for i in 1:L]
    return Int.(hcat(res_msa...)')
end


function reshuffle(msa_in::Array{Int,2}, h::Array{T,2}, J::Array{T,4}; shuffle_units = 0.1, temp = 1, info = false) where {T}
    
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
        m = rand(1:M)
        n = rand(1:M)
        a = msa[i,m]
        b = msa[i,n]
        
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
    
    println("Acceptance rate $(acc_rate) / $(N_shuffles) = $(acc_rate/N_shuffles)")
    
    if info == true
        return msa, ens_evol[1:1000:end], [i for i in 1:(N_shuffles+1)][1:1000:end] ./ unit
    else
        return msa
    end
end
    
            
            
        
        
    

    
    