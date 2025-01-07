Base.@kwdef mutable struct Seqontree <: TreeNodeData # Create a custom data type
    seq::Array{Int,1} = [0,0]
end


function propagate_with_felse!(node::TreeNode{Seqontree}, 
        p, 
        prob::Array{Float64,1}, 
        m::Float64,
        site::Int,
        rng::Xoshiro;
        q = 21) 
    
    if isempty(node.child)
        return 0
    end
    
    for n in node.child
                      
        t = branch_length(n)
        w = exp(- m*t)
        for aa in 1:q
            prob[aa] = (1-w)*p[aa]
        end
        
        old_aa = data(node).seq[site]
        prob[old_aa] += w 
        
        loc_sample_phylo!(rng, prob, data(n).seq, site)
        
        propagate_with_felse!(n, p, prob, m, site, rng)
    end
end

function FelsensteinSampler(start_seq::Array{Int,1}, S::Array{Float64,2}, T, m::Float64; q::Int = 21, verbose::Bool = true)
"""
Function to run the Felsenstein algorithm on a given tree. 
Input:
    
    - S: Stationary distribution of aminoacid on each site (NOT obtained with Z, but from all the natural sequences, that may be different from Z)
    - T: Tree we want to run over (using Pierre TreeTools package)
    - m: Value of μ (timescale of evolution)
    - q (optional): Alphabet length (21 in proteins)
    - verbose (optional): Prints each time a site is inferred (true)
Output: 
   
"""
    for l in keys(T.lnodes)
        data!(T[l], Seqontree(seq = copy(start_seq)))
    end

    L,q = size(S)    
    rng = random_generators(L)
   
    @tasks for site in 1:L
        prob = zeros(q)
        propagate_with_felse!(T.root, view(S,site,:), prob, m, site, rng[site], q = q)
        if verbose == true 
            println("Evolution of site $(site) completed")
        end
    end

    return T
end





function FelsensteinSampler(start_seq::Array{Int,1}, FileNat::String, FileTree::String, m::Float64; verbose::Bool = false, eps::Float64 = 10^-5)
"""
Function to run the Felsenstein algorithm given the files containing the sequences and tree. 
Input:
    - FileNat: File containing the natural sequences (from which we estimate the stationary probability)
    - FileTree: File containing the tree
    - m: Value of μ (timescale of evolution)
    - verbose (optional): Prints each time a site is inferred (true)
    - eps (optional): Pseudocount value (10^-5) 
Output: 
    
"""  
    _, NatMSA, _, _, q = read_fasta(FileNat,1.0,0.2,false)
    StationaryProb = compute_empirical_freqs(Int.(NatMSA), q; eps = eps)

    MyTree = read_tree(FileTree, node_data_type = Seqontree)

    return FelsensteinSampler(start_seq, StationaryProb, MyTree, m; q = q, verbose = verbose)

    
end