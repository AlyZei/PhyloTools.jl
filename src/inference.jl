
function extract_distances(T, start_seq; gap_option = false)
    ds = [];
    ts = [distance(T.root, a) for a in leaves(T)];
    for a in keys(T.lleaves)
        push!(ds, ham_dist(start_seq, data(T[a]).seq))
    end
    
    return ts, ds
end

function extract_distances(FileTree::String, FileSequences::String, start_seq; gap_option = false)
    Z = read_fasta_dict(FileSequences)
    T = read_tree(FileTree, node_data_type = Seqontree)
    for leaf in leaves(T)
              data(T[label(leaf)]).seq = Z[label(leaf)] 
    end
       
    return extract_distances(T, start_seq, gap_option = gap_option)
end

function infer_mu(T; gap_option = true)
    ds = [];ts = [];
    for a in keys(T.lleaves)
        for b in keys(T.lleaves)
            push!(ts, distance(T, a, b));
            if gap_option == true
                mask = (data(T[a]).seq .!= 21) .& (data(T[b]).seq .!= 21)
                push!(ds, ham_dist(data(T[a]).seq[mask], data(T[b]).seq[mask]))
            else
                push!(ds, ham_dist(data(T[a]).seq, data(T[b]).seq))
            end
        end 
    end
    
    model(x, p) = p[1] *(1 .- exp.(-p[2]*x)) 
    p0 = [0.5, 0.5]
    fiti = curve_fit(model, ts, ds, p0)
    
    return fiti.param[2]
end


function infer_mu(FileTree::String, FileSequences::String; gap_option = true)
    Z = read_fasta_dict(FileSequences)
    T = read_tree(FileTree, node_data_type = Seqontree)
    for leaf in leaves(T)
              data(T[label(leaf)]).seq = Z[label(leaf)] 
    end
       
    return infer_mu(T, gap_option = gap_option)
end

function get_anc_loglike(Z, AncestorProb::Array{Array{Float64,1},1}, TransitionProb::Array{Array{Float64,2},1}, W::Array{Float64,2}, Ts, n::String, m::Float64, L::Int; q::Int = 21)
    
    #clean tree prob and put deltas on leaves
    @tasks for s in 1:L
        #define function for both for cycles
        for node in nodes(Ts[s])
            for a in 1:q
                data(Ts[s][label(node)]).prob[a] = 0.
            end
        end
    
        for leaf in leaves(Ts[s])
            for a in 1:q
                data(Ts[s][label(leaf)]).prob[a] = 0.
            end 
            data(Ts[s][label(leaf)]).prob[Z[label(leaf)][s]] = 1.
        end
    end
    
    like = zeros(Float64, L)
    #W_new = hcat([softmax(W[:,s]) for s in 1:L]...) 
    W_new = W
    
    @tasks for s in 1:L
        for a in 1:q
            AncestorProb[s][a] = 0.
        end
        logZn = 0.
        #AncestorProb[s], logZn = FelsensteinSingle2!(logZn, TransitionProb[s], Ts[s], n, m, W_new[:,s], q)
        a,b = FelsensteinSingle2!(logZn, TransitionProb[s], Ts[s], n, m, view(W_new,:,s), q)
        normal = sum(a .* W_new[:,s])
        like[s] += log(normal)+b 
    end
    
    #println(like)
       
    return -sum(like)
end


function tree_inizialization(Z, T, s::Int; q::Int = 21)
   
    for node in nodes(T)
        # Create a new array for the probability distribution
        data(T[label(node)]).prob = zeros(Float64, q)
    end
    
    for leaf in leaves(T)
        prob = [i == Z[label(leaf)][s] ? 1.0 : 0.0 for i in 1:q]
        data(T[label(leaf)]).prob = prob
    end
end


function get_loglike_singlesite(Z, AncestorProb::Array{Float64,1}, TransitionProb::Array{Float64,2}, W::Array{Float64,1}, T, s::Int, n::String, m::Float64, L::Int; q::Int = 21)
    
    tree_inizialization(Z, T, s; q = q)
    
    like = 0.
    W_new = softmax(W) 
    
       
    logZn = 0.
    a,b = FelsensteinSingle2!(logZn, TransitionProb, T, n, m, W_new, q)
    normal = sum(a .* W_new)
    like += log(normal)+b 

    return -like 
    
end





function optimizer(Z, AncestorProb::Array{Float64,1}, TransitionProb::Array{Float64,2}, T, s::Int, n::String, m::Float64, L::Int; 
    q::Int = 21,
    iterations::Int=1000, 
    tol=1e-4,
    optimizer=LBFGS())
    
    
    #W0 = ones(q) ./ q
    W0 = rand(q) 
    W0 ./= sum(W0)
    wraplike! = W -> get_loglike_singlesite(Z, AncestorProb, TransitionProb, W, T, s, n, m, L; q = q)
    res = Optim.optimize(wraplike!, 
        W0,
        optimizer, 
        Optim.Options(iterations=iterations, f_tol=tol))
    
    return res.minimizer    
end





function run_inference(FileNat::String, FileTree::String, FileSequences::String, m::Float64; 
        verbose::Bool = false, 
        eps::Float64 = 1e-5,
        η::Float64 = 1e-1,
        tol = 1e-4,
        opt_iter = 1000)
    
     _, NatMSA, _, _, q = read_fasta(FileNat,1.0,0.2,false)
    L,M = size(NatMSA)
    
    
    
    StationaryProb = compute_empirical_freqs(Int.(NatMSA), q; eps = eps)
    W = ones(q,L) ./ q 
    
    #W = hcat([StationaryProb[i,:] for i in 1:L]...);
    W_new = deepcopy(W);
    
    LeavesSequences = read_fasta_dict(FileSequences)
    
    T = read_tree(FileTree, node_data_type = ProbabilityOnTree)
    Ts = [copy(T) for _ in 1:L]
    
    AncestorProb = [zeros(Float64,q) for _ in 1:L]
    TransitionProb = [zeros(Float64,q,q) for _ in 1:L]
    
    Z = read_fasta_dict(FileSequences)
    
    like = zeros(L)
    
    n = label(root(T))
    
    @tasks for s in 1:L
        for leaf in leaves(Ts[s])
            ρ = zeros(Float64,q)
            ρ[Z[label(leaf)][s]] = 1.
            data!(leaf, ProbabilityOnTree(prob = ρ))
        end
    end
    
          
    @tasks for s in 1:L
        
        opt_W = optimizer(Z, AncestorProb[s], TransitionProb[s], Ts[s], s, n, m, L; q = q, iterations = opt_iter, tol = tol)
        for a in 1:q
            W[a,s] = opt_W[a]
        end
    end
        
    #W_new = W; 
    W_new = hcat([softmax(W[:,s]) for s in 1:L]...)
    InfProb = Float64.(W_new');
    
    for s in 1:L
        like[s] = -get_loglike_singlesite(Z, AncestorProb[s], TransitionProb[s], W[:,s], Ts[s], s, n, m, L; q = q)
    end
    
    if verbose == true
        println("Like $(round(sum(like), digits = 9)) Pearson $(round(cor(InfProb[:], StationaryProb[:]), digits = 5))")
    end
    
    return W_new
end

    



