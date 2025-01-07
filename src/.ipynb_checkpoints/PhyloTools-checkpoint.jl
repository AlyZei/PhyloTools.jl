module PhyloTools

using Revise
using OhMyThreads
using LinearAlgebra
using PyPlot
using TreeTools
using DCAUtils
using Random
using Statistics
using FastaIO
using LsqFit
using Zygote
using StatsBase
using Optim

import Flux: DataLoader, Adam, gradient, softmax
import Flux.Optimise: update! 
import Flux.Optimisers: setup

import KitMSA: fasta2matrix, letter2num, num2letter, extract_params, read_par_BM, set_max_field_to_0

include("reconstruction.jl")
include("inference.jl")
include("sampler.jl")
include("utils.jl")
include("hamming.jl")
include("run_ontree.jl")
include("utils2.jl")
include("energy.jl")
include("reshuffling.jl")

export Felsenstein, Felsenstein2, FelsensteinSampler, seqs_from_leaves, seqs_from_nodes, read_fasta, compute_empirical_freqs, run_inference, leavestofasta, infer_mu, ham_dist, run_evolution_ontree, pairwise_ham_dist, Seq, msa_from_leafs, sampling_ANC, energy, reshuffle, moving_average


end
