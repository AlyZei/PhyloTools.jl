aminoalphabet = Dict{String,Int64}()
nucleoalphabet = Dict{String,Int64}()
aminoalphabet = Dict('A' => 1, 'B' => 21, 'C' => 2, 'D' => 3, 'E' => 4, 'F' => 5, 'G' => 6, 'H' => 7, 'I' => 8, 'J' => 21, 'K' => 9, 'L' => 10, 'M' => 11, 'N' => 12, 'O' => 21, 'P' => 13, 'Q' => 14, 'R' => 15, 'S' => 16, 'T' => 17, 'U' => 21, 'V' => 18, 'W' => 19, 'X' => 21, 'Y' => 20, 'Z' => 21, '-' => 21)
nucleoalphabet = Dict('A' => 1, 'C' => 2, 'G' => 3, 'U' => 4, '-' => 5)

aminonum = Dict{Int64, String}()
nucleonum = Dict{Int64, String}()
aminonum = Dict(1 => 'A', 2 => 'C', 3 => 'D', 4 => 'E', 5 => 'F', 6 => 'G', 7 => 'H', 8 => 'I', 9 => 'K', 10 => 'L', 11 => 'M', 
                12 => 'N', 13 => 'P', 14 => 'Q', 15 => 'R', 16 => 'S', 17 => 'T', 18 => 'V', 19 => 'W', 20 => 'Y', 21 => '-')
nucleonum = Dict(1 => 'A', 2 => 'C', 3 => 'G', 4 => 'U', 5 => '-')


codontoamino = Dict("AAA" => 'K', "AAC" => 'N', "AAG" => 'K', "AAT" => 'N',
                    "ACA" => 'T', "ACC" => 'T', "ACG" => 'T', "ACT" => 'T',
                    "AGA" => 'R', "AGC" => 'S', "AGG" => 'R', "AGT" => 'S',
                    "ATA" => 'I', "ATC" => 'I', "ATG" => 'M', "ATT" => 'I',
                    "CAA" => 'Q', "CAC" => 'H', "CAG" => 'Q', "CAT" => 'H',
                    "CCA" => 'P', "CCC" => 'P', "CCG" => 'P', "CCT" => 'P',
                    "CGA" => 'R', "CGC" => 'R', "CGG" => 'R', "CGT" => 'R',
                    "CTA" => 'L', "CTC" => 'L', "CTG" => 'L', "CTT" => 'L',
                    "GAA" => 'E', "GAC" => 'D', "GAG" => 'E', "GAT" => 'D',
                    "GCA" => 'A', "GCC" => 'A', "GCG" => 'A', "GCT" => 'A',
                    "GGA" => 'G', "GGC" => 'G', "GGG" => 'G', "GGT" => 'G',
                    "GTA" => 'V', "GTC" => 'V', "GTG" => 'V', "GTT" => 'V',
                    "TAC" => 'Y', "TAT" => 'Y',
                    "TCA" => 'S', "TCC" => 'S', "TCG" => 'S', "TCT" => 'S',
                    "TGC" => 'C', "TGG" => 'W', "TGT" => 'C',
                    "TTA" => 'L', "TTC" => 'F', "TTG" => 'L', "TTT" => 'F',
                    "---" => '-')

function translate_codontoamino(nucleosequence::String)   
    L = length(nucleosequence)
    if L % 3 == 0
        AminoLength = L ÷ 3
    else
        println("Length is not a multiple of 3")
        return 0
    end
    
    codons = [nucleosequence[n*3+1:n*3+3] for n in 0:AminoLength-1]
    
    aminostring = ""
    for i in 1:AminoLength
        aminostring *= codontoamino[codons[i]]
    end

    return aminostring

end


function H_distance(seq1::Vector, seq2::Vector)
    ### Function computing the Hamming distance between two vectors

    L = length(seq1)
    @assert L == length(seq2) "Error: the two vectors do not have the same length"
    d = 0
    for i in 1:L
        if seq1[i] != seq2[i]
            d += 1
        end
    end
    return d
end


function H_distance_nogaps(seq1::Vector, seq2::Vector)
    ### Function computing the Hamming distance between two vectors without considering gaps

    L = length(seq1)
    @assert L == length(seq2) "Error: the two vectors do not have the same length"
    d = 0
    for i in 1:L
        if seq1[i] != 21 && seq2[i] != 21
            if seq1[i] != seq2[i]
                d += 1
            end
        end
    end
    return d
end


function read_fasta(filename::AbstractString, max_gap_fraction::Real, theta::Any, remove_dups::Bool)
    """
    filename: name of the fasta file 
    max_gap_fraction : maximum fraction of gaps per sequence
    theta : parameter of the reweigthing# function read_fasta_alignment(filename::AbstractString, max_gap_fraction::Real)
    remove_dups : remove duplicate sequences
    """
    Z = read_fasta_alignment(filename, max_gap_fraction) #reads file and traslate it immediately to numbers (but the gap is at 21 (?)) - function from the package DCAutils
    if remove_dups
        Z, _ = remove_duplicate_sequences(Z) #removes duplicate sequences - function from the package DCAutils
    end
    N, M = size(Z)
    q = round(Int, maximum(Z)) #number of symbols, I think
    W, Meff = compute_weights(Z, theta) #computes the weights of each sequence based on similarity with the others, and gives the effective number of sequences Meff
    return W, Z, N, M, q
end


function read_fasta_dict(FileSeq::String)
    """
    Reads the sequences in the input file and saves them in a dictionary with the name of the sequences as keys and the sequences (in numeric form) as values
    """
    seq_text = readlines(FileSeq)
    Z_text = Dict()
    seq = ""
    name = ""
    for i in 1:length(seq_text)
        if seq_text[i][1] == '>'
            name = split(seq_text[i],">")[2]
            seq = ""
        else
            seq *= seq_text[i]
        end
        Z_text[name] = seq
    end

    Z = Dict{}()
    for key in keys(Z_text)
        num_seq = [aminoalphabet[Z_text[key][i]] for i in 1:length(Z_text[key])]
        # if !(num_seq in values(Z))
        #     Z[key] = num_seq
        # end
        Z[key] = num_seq
    end
    return Z
end


function compute_empirical_freqs(Z::Array{Int64,2}, q::Int; eps::Float64 = 0.0) 
    """
    Function to compute the frequency of occurrence of each aminoacid in the MSA.

    parameters:
    - Z: MSA in format Length_of_sequences x Number_of_sequences
    - q: length of the alphabet
    - eps: pseudocount value

    output:
    - f: matrix Length_of_sequences x q with the frequence of each aminoacid at each site

    """
    L, M = size(Z)
    f = zeros(L, q)
    for i in 1:L
        for s in 1:M
            f[i, Z[i, s]] += 1
        end
    end
    f ./= M

    for s in 1:L
        for a in 1:q
            f[s,a] += eps
        end
        f[s,:] ./= sum(f[s,:])
    end
    return f
    
end




function loc_softmax_phylo!(out::AbstractArray{T}, x::AbstractArray{T}) where {T}
    max_ = T(maximum(x))
    if isfinite(max_)
        @fastmath out .= exp.(x .- max_)
    else
        _zero, _one, _inf = T(0), T(1), T(Inf)
        @fastmath @. out = ifelse(isequal(max_,_inf), ifelse(isequal(x,_inf), _one, _zero), exp(x - max_))
    end
    #tmp = dims isa Colon ? sum(out) : sum!(max_, out)
    out ./= sum(out)#tmp
end

loc_softmax_phylo!(x::AbstractArray) = loc_softmax_phylo!(x, x)

function loc_sample_phylo!(rng, wv, dest_seq::Array{Ti, 1}, seq_site::Int) where {Ti<:Integer}
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


function random_generators(num_generators::Int) 
    rng_array = []
    for seed in 1:num_generators
        push!(rng_array, Random.Xoshiro(seed))
    end
    return rng_array
end

function seqs_from_leaves(tree)
    msa = []; for a in keys(tree.lleaves)
    push!(msa, data(tree[a]).seq) end; msa = hcat(msa...);
    return msa
end

function seqs_from_nodes(tree)
    msa = []; for a in keys(tree.lnodes)
    push!(msa, data(tree[a]).seq) end; msa = hcat(msa...);
    return msa
end


function leavestofasta(path, tree)
    n_seq = length(tree.lleaves)
    FastaWriter(path, "w") do file
        for a in keys(tree.lleaves) 
            writeentry(file, "$(a)", vec2string(tree["$(a)"].data.seq[:]))
        end
    end
end


"""
    vec2string(v)

    Takes as input a vector of integers (representing an amino acids)
    and returns a list of characters, the corresponding amino acids. 
    In this case with the convention "1 2 .. 21" == "A C .. -"
"""


function vec2string(v)
    s = ""
    for i in v
        s = s*num2letter(i)
    end
    return s
end


function moving_average(data, window_size)
    n = length(data)
    ma = Vector{Float64}(undef, n - window_size + 1)
    for i in 1:(n - window_size + 1)
        ma[i] = mean(data[i:i + window_size - 1])
    end
    return ma
end