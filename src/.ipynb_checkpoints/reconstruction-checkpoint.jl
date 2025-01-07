#jinka

#################################################################
############### FELSENSTEIN ALGORITHM ###########################
#################################################################

### Define the type of data that I want to assign to the tree nodes. In my case it's just a vector of probabilities. Can it be done better?
Base.@kwdef mutable struct ProbabilityOnTree <: TreeNodeData # Create a custom data type
    prob :: Array{Float64,1} = zeros(Float64,21)
end

function FelsensteinSingle!(T, n, m, p, q)
    
    LogNodeProb = zeros(Float64,q)
    for child in children(T[n])
        if isleaf(child)
            ChildProb = data(child).prob
        else
            ChildProb = FelsensteinSingle!(T, label(child), m, p, q) # Iterate over the child
        end
        t = branch_length(child)
        if t == 0
            t += 10^-5 # ADD CORRECTION FOR 0 BRANCH LENGTH CASE. TO BE MODIFIED CORRECTLY
        end 
        W = exp(- m * t) 
        TransitionProb = repeat(p', q) * (1-W) + W * I(q) # This is equivalent to P(b|a,t) = δ_{a,b} * exp(- μt) + (1-exp(-μt)) w(b)
        LogNodeProb += log.(TransitionProb * ChildProb) # This is equivalent to log(\sum_b P(b|a,t) ρ(b))
    end
    NodeProb = exp.(LogNodeProb) ./ sum(exp.(LogNodeProb))
    data!(T[n], ProbabilityOnTree(prob = NodeProb))
    return NodeProb
end

#=function FelsensteinSingle2!(logZn::Float64, TransitionProb::Array{Float64,2}, T, n::String, m::Float64, p, q::Int)
   
    LogNodeProb = zeros(Float64,q)
    for a in 1:q
        for b in 1:q
            TransitionProb[a,b] = 0.
        end
    end 
    
        
    for child in children(T[n])
        if isleaf(child)
            ChildProb = data(child).prob
            logZn += 0.
        else
            ChildProb, logZn = FelsensteinSingle2!(logZn, TransitionProb, T, label(child), m, p, q)
        end
        
        t = branch_length(child)
        if t == 0
            t += 10^-5 # ADD CORRECTION FOR 0 BRANCH LENGTH CASE. TO BE MODIFIED CORRECTLY
        end 
        esp = exp(-m*t)
        for aa in 1:q
            for bb in 1:q
                TransitionProb[aa,bb] = (1-esp)*p[bb]
                if aa == bb
                    TransitionProb[aa,bb] += esp
                end
            end
        end
               
        LogNodeProb += log.(TransitionProb * ChildProb) # This is equivalent to log(\sum_b P(b|a,t) ρ(b))        
    end
        
    NodeProb = zeros(q); 
    conto = 0.
    for a in 1:q
        num = exp(LogNodeProb[a])
        conto += num
        #NodeProb[a] = num
    end
    for a in 1:q
        NodeProb[a] = NodeProb[a] / conto
    end
    
    logZn += log(conto)
    data!(T[n], ProbabilityOnTree(prob = NodeProb))
    
    return NodeProb, logZn
end =#


function FelsensteinSingle2!(logZn::Float64, TransitionProb::Array{Float64,2}, T, n::String, m::Float64, p, q::Int)
    # Initialize LogNodeProb without in-place modification later
    LogNodeProb = zeros(Float64, q)
    
    # Process children nodes
    for child in children(T[n])
        if isleaf(child)
            ChildProb = data(child).prob
            logZn += 0.0  # This doesn't change anything but retains structure
        else
            # Recursive call
            ChildProb, logZn = FelsensteinSingle2!(logZn, TransitionProb, T, label(child), m, p, q)
        end
        
        # Calculate the transition probabilities using a temporary array
        t = branch_length(child)
        if t == 0
            t += 1e-5  # Adjust for zero branch length
        end
        esp = exp(-m * t)
        
        # Compute TransitionProb temporarily without in-place updates
        newTransitionProb = [(1 - esp) * p[bb] + (aa == bb ? esp : 0.0) for aa in 1:q, bb in 1:q]
        #TransitionProb .= newTransitionProb  # Update the pre-allocated TransitionProb
        
        # Update LogNodeProb without mutating in-place
        LogNodeProb = LogNodeProb .+ log.(newTransitionProb * ChildProb)
    end
    
    # Compute normalization factor
    conto = sum(exp.(LogNodeProb))  # No in-place operations
    NodeProb = softmax(LogNodeProb)  # Softmax instead of manual normalization
    logZn += log(conto)
    
    # Update tree data
    data!(T[n], ProbabilityOnTree(prob = NodeProb))
    
    return NodeProb, logZn
end



function Felsenstein(Z::AbstractDict, S::Array{Float64,2}, T, n, m::Float64; q::Int = 21, verbose::Bool = true, eps::Float64 = 10^-5)

    @assert in(n,T) "The selected node is not present inside the tree"
    @assert abs(sum(S[1,:])-1) < 10^-8 "Stationary probability is not normalized, $(sum(S[1,:]))"
    @assert length(leaves(T)) == length(Z) "Number of leaves in the tree T does not match the number of sequences in Z"
    
    M = length(Z)
    L = length(collect(values(Z))[1])
    
    AncestorProb = zeros(Float64,(L,q))
    Ancestor = zeros(Int,L)
    
    @tasks for s in 1:L
        _T = copy(T)
        for leaf in leaves(_T)
            ρ = zeros(Float64,q)
            ρ[Z[label(leaf)][s]] = 1
            data!(leaf, ProbabilityOnTree(prob = ρ))
        end
        #AncestorProb[s,:] = FelsensteinSingle!(_T, n, m, S[s,:], q)
        AncestorProb[s,:] = FelsensteinSingle!(_T, n, m, S[s,:], q) .* S[s,:]
        #AncestorProb[s,:] ./= sum(AncestorProb[s,:])
        Ancestor[s] = argmax(AncestorProb[s,:])
        if verbose == true 
            println("Reconstruction of site $s completed, a[$s] = $(Ancestor[s])")
        end
    end

    return Ancestor, AncestorProb
end

function Felsenstein2(Z::AbstractDict, S::Array{Float64,2}, T, n::String, m::Float64; q::Int = 21, verbose::Bool = false, eps::Float64 = 10^-5)
    
    @assert in(n,T) "The selected node is not present inside the tree"
    @assert abs(sum(S[1,:]) - 1) < 1e-8 "Stationary probability is not normalized, $(sum(S[1,:]))"
    @assert length(leaves(T)) == length(Z) "Number of leaves in the tree T does not match the number of sequences in Z"
    
    M = length(Z)
    L = length(collect(values(Z))[1])
    
    AncestorProb = zeros(Float64,(L,q))
    Ancestor = zeros(Int,L)
    AncLike = zeros(Float64, L)
    
    @tasks for s in 1:L
        _T = copy(T)
        LogNodeProb = zeros(Float64,q)
        TransitionProb = zeros(Float64,q,q)
        for leaf in leaves(_T)
            ρ = zeros(Float64,q)
            ρ[Z[label(leaf)][s]] = 1
            data!(leaf, ProbabilityOnTree(prob = ρ))
        end
        #AncestorProb[s,:] = FelsensteinSingle2!(TransitionProb, _T, n, m, view(S,s,:), q)
        logZn = 0.
        a,b = FelsensteinSingle2!(logZn, TransitionProb, _T, n, m, view(S,s,:), q)
        AncestorProb[s,:] = a .* S[s,:]
        normal = sum(AncestorProb[s,:])
        AncestorProb[s,:] ./= normal
        Ancestor[s] = argmax(view(AncestorProb,s,:))
        if verbose == true 
            println("Reconstruction of site $s completed, a[$s] = $(Ancestor[s])")
        end
        AncLike[s] += normal*exp(b)
    end
    return Ancestor, AncestorProb, AncLike
end


function Felsenstein(FileNat::String, FileTree::String, FileSequences::String, m::Float64; verbose::Bool = false, eps::Float64 = 10^-5)

    println("Saverio version")
    _, NatMSA, _, _, q = read_fasta(FileNat,1.0,0.2,false)
    StationaryProb = compute_empirical_freqs(Int.(NatMSA), q; eps = eps)
    
    LeavesSequences = read_fasta_dict(FileSequences)
    MyTree = read_tree(FileTree, node_data_type = ProbabilityOnTree)

    AncestorSequence, AncestorProbability = Felsenstein(LeavesSequences, StationaryProb, MyTree, label(root(MyTree)), m; q = q, verbose = verbose, eps=eps)

    return AncestorSequence, AncestorProbability
end

function Felsenstein2(FileNat::String, FileTree::String, FileSequences::String, m::Float64; verbose::Bool = false, eps::Float64 = 10^-5)
  
    _, NatMSA, _, _, q = read_fasta(FileNat,1.0,0.2,false)
    StationaryProb = compute_empirical_freqs(Int.(NatMSA), q; eps = eps)
    
    LeavesSequences = read_fasta_dict(FileSequences)
    MyTree = read_tree(FileTree, node_data_type = ProbabilityOnTree)

    AncestorSequence, AncestorProbability, AncLike = Felsenstein2(LeavesSequences, StationaryProb, MyTree, label(root(MyTree)), m; q = q, verbose = verbose, eps=eps)

    return AncestorSequence, AncestorProbability, AncLike
end



