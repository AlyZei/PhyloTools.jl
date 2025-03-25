Base.@kwdef mutable struct SeqAmino <: TreeNodeData # Create a custom data type
    seq::Array{Int,1} = [0,0]
end

Base.copy(s::SeqAmino) = SeqAmino(copy(s.seq))

struct StorageAmino{T}
    tree::Tree{SeqAmino}
    generator::Xoshiro
    log_prob::Vector{T}
end


function StorageAmino(start_seq::Array{Int,1}, tree_file::String, generator::Xoshiro; T::DataType = Float64)
    
    tree = read_tree(tree_file, node_data_type = SeqAmino)
    for l in keys(tree.lnodes)
        data!(tree[l], SeqAmino(seq = copy(start_seq)))
    end
    log_prob = T.([0.,0.]); 
    StorageAmino{T}(tree, generator, log_prob)
end



function log_prob_tree_amino!(log_prob::Array{T,1}, 
        seq::Array{Int,1}, 
        h::Array{T,2}, 
        J::Array{T,4},
        temp::T,
        seq_site::Int,
        L::Int;
        q::Int = 21) where {T}
    
    @inbounds for a in 1:q 
        log_prob[a] = h[a,seq_site] 
        @inbounds for j in 1:L
            log_prob[a] += J[seq[j], j, a, seq_site]
        end
    end
    
end
      

function loc_sample_amino!(rng::AbstractRNG, wv, dest_seq::Array{Ti, 1}, seq_site::Int) where {Ti<:Integer}
    t = rand(rng) * sum(wv)
    n = length(wv)
    i = one(Ti)
    cw = wv[1]
    while cw < t && i < n
        i += one(Ti)
        @inbounds cw += wv[i]
    end
        
    dest_seq[seq_site] = i

end


function get_prob_tree_amino!(seq::Array{Int,1}, 
        log_prob::Array{T,1},
        h::Array{T,2}, 
        J::Array{T,4},
        seq_site::Int,
        rng::Xoshiro,
        temp::T,
        L::Int;
        q::Int = 21) where {T} 
    
    log_prob_tree_amino!(log_prob, seq, h, J, temp, seq_site, L, q = q)
    log_prob ./= temp
    loc_softmax!(log_prob)
    loc_sample_amino!(rng, log_prob, seq, seq_site)   
end


function prob_cond_tree_amino!(chain::StorageAmino{Float64}, 
        seq::Array{Int,1}, 
        h::Array{T,2}, 
        J::Array{T,4},
        seq_site::Int,
        temp::T,
        L::Int;
        q::Int = 21) where {T}
    
    get_prob_tree_amino!(seq, chain.log_prob, h, J, seq_site, chain.generator, temp, L, q = q)    
end



function run_gibbs_sampling_tree_amino!(chain::StorageAmino{Float64}, 
        seq::Array{Int,1}, 
        h::Array{T,2}, 
        J::Array{T,4}, 
        temp::T,
        L::Int;
        q::Int = 21) where {T}

    seq_site = rand(1:L)    
    prob_cond_tree_amino!(chain, seq, h, J, seq_site, temp, L, q = q)
end


function assign_sequences_amino!(node::TreeNode{SeqAmino}, 
        chain::StorageAmino{Float64},
        h::Array{T,2}, 
        J::Array{T,4}, 
        temp::T,
        mu::Float64,
        p::Float64,
        L::Int;
        q::Int = 21) where {T}
    
    if isempty(node.child)
        return 0
    end
    
    for a in node.child
        
        for i in 1:L
            data(a).seq[i] = data(a.anc).seq[i] 
        end
                
        estimate_steps = L*mu*branch_length(a)
        steps = floor(Int, estimate_steps)
        if rand() < (estimate_steps - steps)
            steps += 1
        end
        
        for _ in 1:steps
            #sampling gibbs with probability p
            run_gibbs_sampling_tree_amino!(chain, data(a).seq, h, J, temp, L, q = q)
            #println("Node $(a.label) $(sum(data(a).seq .== data(a.anc).seq))")
        end
        
        assign_sequences_amino!(a, chain, h, J, temp, mu, p, L, q = q)
    end
end


function run_evolution_ontree_amino(start_seq::Union{Array{Int,1}, Array{Int,2}}, tree_file::String, h::Array{T,2}, J::Array{T,4};
        temp::Float64 = 1.0, 
        mu::Float64 = 1.0,
        p::Float64 = 0.5, 
        verbose = false,
        q::Int = 21) where {T}
    
    
    L = size(start_seq, 1)
    if (size(J,1) !== size(J,3)) || (size(J,2) !== size(J,4))
        error("Size of J should be (q,L,q,L)")
    elseif (size(J,2) !== L) || (size(h,2) !== L)
        error("Length of sequences different from length of parameters")
    end
    
    
    
    temp = T(temp)
        
    N_trees = size(start_seq,2)    
    
    if N_trees == 1
        rng = random_gens(N_trees+1)  
        chains = [StorageAmino(start_seq, tree_file, rng[n]) for n in 1:N_trees+1]
        assign_sequences_amino!(chains[1].tree.root, chains[1], h, J, temp, mu, p, L, q = q)
        return chains[1].tree
        
    elseif N_trees > 1
        
        trees = []
        rng = random_gens(N_trees)  
        chains = [StorageAmino(start_seq[:,n], tree_file, rng[n]) for n in 1:N_trees]
        @tasks for n in 1:N_trees
            assign_sequences_amino!(chains[n].tree.root, chains[n], h, J, 
                temp, mu, p, L, q = q)
        end 

        return [chains[i].tree for i in 1:N_trees]
        
    elseif (N_trees == 0) || (N_trees < 0)
        error("N_trees must be >= 1")
    end
end
    



