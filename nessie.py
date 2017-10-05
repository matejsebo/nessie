#!/usr/bin/env python

######################## 
######## NESSIE ########
######################## 
# Author: Matej Sebo   #
# Language: Python 2.7 #
########################
VERSION = 3.8

# Standard modules
import os, sys, bisect
import random
import math as m
import argparse
import pickle
import itertools as it
import collections as cl
from multiprocessing import Pool

# Nonstandard modules
import numpy as np

# My modules
from slyce import SlyceList, Slyce

# Call "python nessie.py -x -i Neurospora_crassa_OR74A_FungiDB-3.1.fasta -o Neurospora_crassa_OR74A_isolated.fasta -s 100000"
# to generate the fasta file with isolated scaffolds. 

# Call "python nessie.py -i mygenomes/Neurospora_crassa_OR74A_FungiDB-3_isolated.fasta -tm thallic_mat.csv -sm sub_mat.csv -pt confirmed_rearrangement_rates.tsv"
# to run the simulation. 

# Default values for our parameters -----------------------
OUT_DIRECTORY = 'nessie_out/'
NUM_PROCESSORS = 8 
DEPTH = 3 #log_2(number of evolved output genomes we want)
MIN_SCAFFOLD_LENGTH = 100000
TIME_LAMBDA = 1.0 # rate of speciation (leave at 1.0)

SUB_MAT_FILE = "sub_mat.csv"
THALLIC_MAT_FILE = "thallic_mat.csv"
PARTITION_FILE = "confirmed_rearrangement_rates.tsv"

NUM_PARTITIONS = 2 # 2 partition model assumed as a standard...

MEAN_SUB_RATE = 0.004 # just a dummy value

FREQ_REARR_HETERO = 1.596e-6 # VARY this parameter
FREQ_REARR_HOMO = 2.612e-6 # VARY this parameter

# MEAN_FREQ_INV_HOMO = 2.612e-6 # empirically estimated
# SD_FREQ_INV_HOMO = 1.384e-6 # empirically estimated
# MEAN_FREQ_INV_HETERO = 1.596e-6 # empirically estimated
# SD_FREQ_INV_HETERO = 1.496e-6 # empirically estimated

MEAN_INV_SIZE_HETERO = 1438 # empirically estimated
SD_INV_SIZE_HETERO = 603 # empirically estimated
MEAN_INV_SIZE_HOMO = 2468 # empirically estimated
SD_INV_SIZE_HOMO = 746 # empirically estimated

LIKELIHOOD_INV = 0.267            #i #NO PARAMETERS #TODO compute these from rearr rate
LIKELIHOOD_IO_REARR = 0.266       #q #NO PARAMETERS #TODO compute from rearr rate
# does NOT vary with evolutionary time; they depend on the basic o/e rearr rates
# is computable at start of modle

LIKELIHOOD_ASSIST_TRANSLOC = 0.2 #a
LIKELIHOOD_TANDEM = 0.70         #t 

# BASE INFORMATION / CONSTANTS ----------------------------
BLANK = "-"
A = "A"
C = "C"
G = "G"
T = "T"
BASES = [A, C, G, T]
NUM_BASES = len(BASES)
BASE_MAP = {A:0, C:1, G:2, T:3}
acgt = [A, C, G, T]

# UTILITY METHODS -----------------------------------------

def fasta_to_dict(fasta_file):
    f_fasta = open(fasta_file, 'r')
    data = f_fasta.read()

    scaffolds = data[1:].split('>')
    scaffolds = [scaf.split('\n', 1) for scaf in scaffolds]
    try:
        scaffolds = {scaf[0] : scaf[1].replace('\n', '') for scaf in scaffolds}
    except IndexError:
        print fasta_file
        raise
    f_fasta.close()
    return scaffolds

# Return a genome triplist that contains only scaffolds where 
# bp length > MIN_SCAFFOLD_LENGTH. Run this once and keep the output
# file for use in nessie. 
def isolate_chroms(genome_doublist): 
    new_genome = []
    for metadata, seq in genome_doublist.items(): 
        if len(seq) >= MIN_SCAFFOLD_LENGTH:
            print "    Scaffold " + metadata.split(' ')[0] + ' will be included . . .'
            new_genome += [[metadata, seq]]
    genome_len = 0
    processed_genome = []
    for metadata, seq in new_genome: 
        genome_len += len(seq)
        seq = np.array(list(seq))
        seq[seq == 'N'] = BASES[random.randint(0, 3)] # replace all 'N' bases with random bases
        processed_genome += [[metadata, seq, None]]
    return processed_genome

# Export a genome triplist to a file.
def export_genome(genome_triplist, output_file):
    # export to fasta file
    f = open(output_file, 'w')
    print "Writing to " + output_file + " . . ."
    output_text = ""
    for metadata, seq, bins in genome_triplist:
        print bins
        if not bins: 
            output_text = ">" + metadata + '\n'
            f.write(output_text)
            seq.tofile(f)
            f.write('\n')
            print "    Written scaffold " + metadata.split(' ')[0] + " to file . . ."
        else:
            # TODO how to export a triplist WITH BIN INFORMATION to a file
            pass
            
# Define command-line arguments.
def parse_args():
    parser = argparse.ArgumentParser(description= \
        'Nessie v' + str(VERSION) + ': an evolution simulator for Neurospora.')

    parser.add_argument('-p','--processors', \
        help='Number of processors on this machine (for scaffold point mutation ' + \
        'parallelization).', default=NUM_PROCESSORS, required=False)

    parser.add_argument('-i','--in', \
        help='Input file (.fasta)', required=True)
    parser.add_argument('-o','--out', default=OUT_DIRECTORY, \
        help='Output directory.', required=False)

    parser.add_argument('-d','--depth', default=DEPTH, \
        help='Depth of generated phylogenetic tree (defaults to 3).', required=False)

    parser.add_argument('-x', \
        help='Isolate and process (remove N-bases) scaffolds of size -s or greater ' + \
        'from the input file to the output file.', required=False, dest='isolate', action='store_true')
    parser.add_argument('-s','--size', default=MIN_SCAFFOLD_LENGTH, \
        help='Minimum scaffold size for isolation. See -x.', required=False)

    parser.add_argument('-tm','--thallic_mat', \
        help='Source file for the thallic transition matrix  (in *.csv format).', required=False)
    parser.add_argument('-sm','--sub_mat', \
        help='Source file for the base substitution transition matrix (in *.csv format).', required=False)
    parser.add_argument('-pt','--partitions', \
        help='Source file for the partition indices and rearrangement rates (in *.tsv format).', required=False)

    parser.add_argument('-linv','--likelihood_inv', default=LIKELIHOOD_INV, \
        help='Probability(inversion | rearrangement) = "i" in our model.', required=False)
    parser.add_argument('-lio','--likelihood_io_rearr', default=LIKELIHOOD_IO_REARR, \
        help='Pr(duplication | rearrangement) = Pr(deletion | rearrangement) = "q" in our model.', required=False)
    parser.add_argument('-ltan','--likelihood_tandem', default=LIKELIHOOD_TANDEM, \
        help='Probability(tandem duplication | duplication) = "t" in our model.', required=False)
    parser.add_argument('-lat','--likelihood_assist_transloc', default=LIKELIHOOD_ASSIST_TRANSLOC, \
        help='Probability(assisted translocation | rearrangement) = "a" in our model.', required=False)

    parser.add_argument('-ros','--speciation_rate', default=TIME_LAMBDA, \
        help='Rate of speciation. Used to compute evolutionary time (t).', required=False)

    parser.add_argument('-fre','--freq_rearr_hetero', default=FREQ_REARR_HETERO, \
        help='Frequency of rearrangement for heterothallic individuals.', required=False)
    parser.add_argument('-fro','--freq_rearr_homo', default=FREQ_REARR_HOMO, \
        help='Frequency of rearrangement for homothallic individuals.', required=False)

    # parser.add_argument('-mfie','--mean_freq_inv_hetero', default=MEAN_FREQ_INV_HETERO, \
    #     help='Mean frequency of inversion for heterothallic individuals.', required=False)
    # parser.add_argument('-sdfie','--sd_freq_inv_hetero', \
    #     help='Standard deviation of the frequency of inversion for heterothallic '+\
    #     'individuals.', required=False, default=SD_FREQ_INV_HETERO)
    # parser.add_argument('-mfio','--mean_freq_inv_homo', default=MEAN_FREQ_INV_HOMO, \
    #     help='Mean frequency of inversion for homothallic individuals.', required=False)
    # parser.add_argument('-sdfio','--sd_freq_inv_homo', \
    #     help='Standard deviation of the frequency of inversion for homothallic '+\
    #     'individuals.', required=False, default=SD_FREQ_INV_HOMO)

    parser.add_argument('-mise','--mean_inv_size_hetero', default=MEAN_INV_SIZE_HETERO, \
        help='Mean inversion size for heterothallic individuals.', required=False)
    parser.add_argument('-sdise','--sd_inv_size_hetero', \
        help='Standard deviation of the inversion size for heterothallic '\
        'individuals.', required=False, default=SD_INV_SIZE_HETERO)
    parser.add_argument('-mfiso','--mean_inv_size_homo', default=MEAN_INV_SIZE_HOMO, \
        help='Mean inversion size for homothallic individuals.', required=False)
    parser.add_argument('-sdiso','--sd_inv_size_homo', \
        help='Standard deviation of the inversion size for homothallic '\
        'individuals.', required=False, default=SD_INV_SIZE_HOMO)


    args = vars(parser.parse_args())
    args['processors'] = int(args['processors'])
    args['depth'] = int(args['depth'])
    args['size'] = int(args['size'])
    args['likelihood_tandem'] = float(args['likelihood_tandem'])
    args['likelihood_inv'] = float(args['likelihood_inv'])
    args['likelihood_io_rearr'] = float(args['likelihood_io_rearr'])
    args['likelihood_assist_transloc'] = float(args['likelihood_assist_transloc'])
    args['speciation_rate'] = float(args['speciation_rate'])
    args['freq_rearr_hetero'] = float(args['freq_rearr_hetero'])
    args['freq_rearr_homo'] = float(args['freq_rearr_homo'])

    # args['mean_freq_inv_hetero'] = float(args['mean_freq_inv_hetero'])
    # args['sd_freq_inv_hetero'] = float(args['sd_freq_inv_hetero'])
    # args['mean_freq_inv_homo'] = float(args['mean_freq_inv_homo'])
    # args['sd_freq_inv_homo'] = float(args['sd_freq_inv_homo'])

    args['mean_inv_size_hetero'] = float(args['mean_inv_size_hetero'])
    args['sd_inv_size_hetero'] = float(args['sd_inv_size_hetero'])
    args['mean_inv_size_homo'] = float(args['mean_inv_size_homo'])
    args['sd_inv_size_homo'] = float(args['sd_inv_size_homo'])

    if args['processors'] < 1:
        sys.exit("ERROR: Must have >= 1 processors.")

    if args['depth'] < 1:
        sys.exit("ERROR: Must have depth of >= 1.")

    if 'isolate' in args and args['isolate']: 
        if not 'in' in args or not args['in']:
            sys.exit("ERROR: Input file required.")
        if not 'out' in args or not args['out']:
            sys.exit("ERROR: Output file required.")

        anc_genome = fasta_to_dict(args['in'])
    
        # This handles scaffold isolation. 
        print "Isolating/processing scaffolds from original genome . . ."
        out_file = args['out']
        export_genome(isolate_chroms(anc_genome), out_file)
        print "Done! Generated " + out_file
        print "Rerun nessie.py with this file as the input."
        sys.exit(0)
    
    # Otherwise, we will run the evolution simulator.

    if args['likelihood_io_rearr'] > 1 or args['likelihood_io_rearr'] < 0 or \
        args['likelihood_inv'] > 1 or args['likelihood_inv'] < 0 or \
        args['likelihood_tandem'] > 1 or args['likelihood_tandem'] < 0:
        sys.exit("ERROR: All probabilities must be >= 0 and <= 1.")


    if args['likelihood_inv'] < args['likelihood_io_rearr']:
        sys.exit("ERROR: likelihood_inv < likelihood_io_rearr, leading to a negative probability in our model.\n" \
            + "(i < q, thus P(random translocation | rearrangement) = i - q < 0)")
    if abs(args['likelihood_inv']*2 + args['likelihood_io_rearr'] + args['likelihood_assist_transloc'] - 1) > 0.01:
        sys.exit("ERROR: likelihood_inv*2 + likelihood_io_rearr + \n    likelihood_assist_transloc ~= 1.0. \n" + 
            "Please check your rearrangement modeling parameters and try again.")

    if args['size'] < 0:
        sys.exit("ERROR: Must have minimum scaffold size of >= 0.")

    if args['speciation_rate'] < 0:
        sys.exit("ERROR: Speciation rate must be positive.")

    # if args['mean_freq_inv_hetero'] < 0 or args['sd_freq_inv_hetero'] < 0 or \
    #     args['mean_freq_inv_homo'] < 0 or args['sd_freq_inv_homo'] < 0:
    #     sys.exit("ERROR: Frequencies/standard deviations must be positive.")

    if args['freq_rearr_hetero'] < 0 or args['freq_rearr_homo'] < 0 or \
        args['mean_inv_size_hetero'] < 0 or args['sd_inv_size_hetero'] < 0 or \
        args['mean_inv_size_homo'] < 0 or args['sd_inv_size_homo'] < 0:
        sys.exit("ERROR: Frequencies/standard deviations must be positive.")

    return args

args = parse_args()

# required args
IN_FILE = args['in']
THALLIC_MAT_FILE = args['thallic_mat']

SUB_MAT_FILE = args['sub_mat'] # either this (1-partition), or the following two (for multipartition)
THALLIC_MAT_FILE = args['thallic_mat']
PARTITION_FILE = args['partitions']

# optional args
OUT_DIRECTORY = args['out']
NUM_PROCESSORS = args['processors']
DEPTH = args['depth']
MIN_SCAFFOLD_LENGTH = args['size']
TIME_LAMBDA = args['speciation_rate']
IS_ISOLATE = args['isolate']
FREQ_REARR_HETERO = args['freq_rearr_hetero']
FREQ_REARR_HOMO = args['freq_rearr_homo']

# MEAN_FREQ_INV_HOMO = args['mean_freq_inv_homo']
# SD_FREQ_INV_HOMO = args['sd_freq_inv_homo']
# MEAN_FREQ_INV_HETERO = args['mean_freq_inv_hetero']
# SD_FREQ_INV_HETERO = args['sd_freq_inv_hetero']

MEAN_INV_SIZE_HETERO = args['mean_inv_size_hetero']
SD_INV_SIZE_HETERO = args['sd_inv_size_hetero']
MEAN_INV_SIZE_HOMO = args['mean_inv_size_homo']
SD_INV_SIZE_HOMO = args['sd_inv_size_homo']

LIKELIHOOD_INV = args['likelihood_inv']
LIKELIHOOD_TANDEM = args['likelihood_tandem']
LIKELIHOOD_IO_REARR = args['likelihood_io_rearr']
LIKELIHOOD_ASSIST_TRANSLOC = args['likelihood_assist_transloc']


# homo- vs hetero- thallic constants ----------------------
try:
    THALLIC_MAT = np.loadtxt(open(THALLIC_MAT_FILE, "rb"), delimiter=",") # the substitution matrix
except:
    sys.exit("ERROR: Thallic matrix file not found at " + THALLIC_MAT_FILE)
#TODO add dimension checker
# DEFAULT PROVIDED AS thallic_mat.csv: 
# # [[ 0.75  0.25]
# #  [ 0.25  0.75]]
HOMO_MAT = np.matrix([[1.0], [0.0]])
HETERO_MAT = np.matrix([[0.0], [1.0]])
# generate the multiplied vectors. 
OUTPUT_HOMO = THALLIC_MAT.dot(HOMO_MAT)
OUTPUT_HETERO = THALLIC_MAT.dot(HETERO_MAT)

# rearrangement rate partition file import ----------------

try:
    partition_base = np.genfromtxt(PARTITION_FILE, dtype=None, delimiter='\t')[1:] # lop off header
except:
    sys.exit("ERROR: Partition file not found at " + PARTITION_FILE)
#TODO add dimension checker
# indices within partition_base
start = 0
end = 1
min_rate = 2
max_rate = 3

TOP_RATE = 0.002

# indices within BIN_INDICES
from_i = 0 # from which base does this bin refer to?
up_to = 1 # up-to which base does this bin refer to?
partition_index = 2 # which partition does this bin correspond to?

# smooth the rates, determine intervals
num_bins = partition_base.shape[0]
BIN_INDICES = np.zeros(shape=(num_bins,3), dtype=np.int32) # reference library of all the partition bins
PARTITIONS = np.zeros(shape=(NUM_PARTITIONS, )) # mean rearrangement probabilities associated with the bins
BIN_LEN = 7000 #TODO compute or parametrize this!
sum_num_partitions = np.zeros(shape=(NUM_PARTITIONS, ))

#TODO zap all the rates that are clearly bad (centromeres, ends of chroms) (we are currently thresholding only really high rates away)

averaged_rates = [(float(partition_base[i][max_rate]) + \
    float(partition_base[i][min_rate]))/2.0 for i in range(num_bins)]
averaged_rates = [min(a, TOP_RATE) for a in averaged_rates]
sorted_rates = np.array(list(sorted(averaged_rates)))
threshold_indices = [(i_bin+1)*num_bins // NUM_PARTITIONS for i_bin in range(NUM_PARTITIONS-1)]
# print threshold_indices, len(averaged_rates)
# print sorted_rates
thresholds = np.array([sorted_rates[ti] for ti in threshold_indices])
# print thresholds

for i, big_interval in enumerate(partition_base):
    big_interval_len = int(big_interval[end]) - int(big_interval[start])
    av_rearr_rate = (float(partition_base[i][max_rate]) + float(partition_base[i][min_rate]))/2.0
    BIN_INDICES[i][from_i] = 0 if i==0 else BIN_INDICES[i-1][up_to]
    BIN_INDICES[i][up_to] = int(big_interval[start]) + 3 * big_interval_len / 5
    BIN_INDICES[i][partition_index] = np.searchsorted(thresholds, av_rearr_rate)
    #print av_rearr_rate, thresholds, BIN_INDICES[i][partition_index]
    #PARTITIONS[BIN_INDICES[i][partition_index]] += av_rearr_rate
    #sum_num_partitions[BIN_INDICES[i][partition_index]] += 1

PARTITION_LEN = len(partition_base) // NUM_PARTITIONS
# PARTITIONS = np.array(([BIN_INDICES[i*PARTITION_LEN+PARTITION_LEN//2][from_i] + \
#     [BIN_INDICES[i*PARTITION_LEN+PARTITION_LEN//2][to_i])/2.0 for i in range(NUM_PARTITIONS)])
partition_indices = [i*PARTITION_LEN+PARTITION_LEN//2 for i in range(NUM_PARTITIONS)] # inaccurate rearr rate values here will matter; need to strip bad (centromeric/edge) values from the original dataset
# print len(sorted_rates)
# print partition_indices
PARTITIONS = np.array([sorted_rates[pi] for pi in partition_indices])
PARTITION_PROBS = PARTITIONS / sum(PARTITIONS)
# print BIN_INDICES
# print PARTITIONS
# print PARTITION_PROBS

#print PARTITIONS
#sys.exit(0)

# partition handling --------------------------------------
# The following is for point mutation partitioning and can be ignored for now

# # A partition holder "struct"
# class Partition(object):
#     def __init__(self, num_original_bases, alpha, tree_len, \
#         rac, rag, rat, rcg, rct, rgt, pia, pic, pig, pit):
#         #self.base_index_list = base_index_list
#         self.num_original_bases = num_original_bases
#         self.alpha = alpha
#         self.tree_len = tree_len
#         self.rac = rac
#         self.rag = rag
#         self.rat = rat
#         self.rcg = rcg
#         self.rct = rct
#         self.rgt = rgt
#         self.pia = pia
#         self.pic = pic
#         self.pig = pig
#         self.pit = pit
#         self.slyce_list = None # computed later...

# # reads partition data from two files into a list of Partition objects
# # takes in a *.partition index file and an ExaML output file (with partition parameters)
# def read_partitions(index_file, param_file):
#     print "Reading partitions . . ."

#     try:
#         f_indices = open(index_file, 'r')
#     except:
#         sys.exit("ERROR: Multipartition index file (*.partition) not found.\n" + \
#             "Use -ri to provide it if you have not already done so.")
#     try:
#         f_params = open(param_file, 'r')
#     except:
#         sys.exit("ERROR: Multipartition parameter file (ExaML output) not found.\n" + \
#             "Use -rp to provide it if you have not already done so.")

#     indices = f_indices.readlines()
#     params = [l.strip() for l in f_params.readlines() if l.strip() != ""]

#     partitions = []
#     i_partition = 0
#     for i in range(len(indices)):
#         num_original_bases = len([indices[i_partition].split(",")]) - 1
#         # -1 because the first comma-separated entry in the *.partition file is not a base location

#         # we have 13 param-lines for every one index line in the param file...
#         base_index = i_partition * 13

#         # partition header is mod 0
#         alpha = float(params[base_index + 1].split(" ")[-1]) # mod 1
#         tree_len = float(params[base_index + 2].split(" ")[-1]) # mod 2
#         rac = float(params[base_index + 3].split(" ")[-1]) # mod 3
#         rag = float(params[base_index + 4].split(" ")[-1]) # mod 4
#         rat = float(params[base_index + 5].split(" ")[-1]) # mod 5
#         rcg = float(params[base_index + 6].split(" ")[-1]) # mod 6
#         rct = float(params[base_index + 7].split(" ")[-1]) # mod 7
#         rgt = float(params[base_index + 8].split(" ")[-1]) # mod 8
#         pia = float(params[base_index + 9].split(" ")[-1]) # mod 9
#         pic = float(params[base_index + 10].split(" ")[-1]) # mod 10
#         pig = float(params[base_index + 11].split(" ")[-1]) # mod 11
#         pit = float(params[base_index + 12].split(" ")[-1]) # mod 12

#         partitions += [Partition(num_original_bases, alpha, tree_len, \
#             rac, rag, rat, rcg, rct, rgt, pia, pic, pig, pit)]
#         i_partition += 1

#     f_indices.close()
#     f_params.close()
#     num_partitions = i_partition + 1
#     print "Finished reading " + str(num_partitions) + " partitions."
#     return num_partitions, partitions

# base transition matrix (for single partition model) -----
# DEFAULT PROVIDED AS sub_mat.csv: 
# # [[ 0.996  0.001  0.001  0.001  0.001]
# #  [ 0.001  0.996  0.001  0.001  0.001]
# #  [ 0.001  0.001  0.996  0.001  0.001]
# #  [ 0.001  0.001  0.001  0.996  0.001]
# #  [ 0.001  0.001  0.001  0.001  0.996]]

try: # Try reading the substitution matrix
    sm = np.loadtxt(open(SUB_MAT_FILE, "rb"), delimiter=",") # the substitution matrix
    NUM_PARTITIONS = 1
    print "Using single partition model for point mutations . . ."
    a_msr = sm[0][1] + sm[0][2] + sm[0][3]
    c_msr = sm[1][0] + sm[1][2] + sm[1][3]
    g_msr = sm[2][0] + sm[2][1] + sm[2][3]
    t_msr = sm[3][0] + sm[3][1] + sm[3][2]
    MEAN_SUB_RATE = (a_msr + c_msr + g_msr + t_msr) / NUM_BASES
    MEAN_INS_RATE = (sm[4][0] + sm[4][1] + sm[4][2] + sm[4][3]) / NUM_BASES
    MEAN_DEL_RATE = (sm[0][4] + sm[1][4] + sm[2][4] + sm[3][4]) / NUM_BASES

    OUTPUT_A   = np.array([0.0,      sm[0][1], sm[0][2], sm[0][3]]) / a_msr
    OUTPUT_C   = np.array([sm[1][0], 0.0,      sm[1][2], sm[1][3]]) / c_msr
    OUTPUT_G   = np.array([sm[2][0], sm[2][1], 0.0     , sm[2][3]]) / g_msr
    OUTPUT_T   = np.array([sm[3][0], sm[3][1], sm[3][2], 0.0     ]) / t_msr
    OUTPUT_INS = np.array([sm[4][0], sm[4][1], sm[4][2], sm[4][3]]) / MEAN_INS_RATE / NUM_BASES
    OUTPUT_DEL = np.array([sm[0][4], sm[1][4], sm[2][4], sm[3][4]]) / MEAN_DEL_RATE / NUM_BASES

    CUM_OUTPUT_A = np.cumsum(OUTPUT_A)
    CUM_OUTPUT_C = np.cumsum(OUTPUT_C)
    CUM_OUTPUT_G = np.cumsum(OUTPUT_G)
    CUM_OUTPUT_T = np.cumsum(OUTPUT_T)
    CUM_OUTPUT_INS = np.cumsum(OUTPUT_INS)
    CUM_OUTPUT_DEL = np.cumsum(OUTPUT_DEL)

    SUB_OUTPUTS = [OUTPUT_A, OUTPUT_C, OUTPUT_G, OUTPUT_T]
    CUM_OUT_MAP = {A:CUM_OUTPUT_A, C:CUM_OUTPUT_C, G:CUM_OUTPUT_G, T:CUM_OUTPUT_T}
except: # If substitution matrix is not found . . .
    NUM_PARTITIONS, PARTITIONS = \
        read_partitions(PARTITION_INDEX_FILE, PARTITION_PARAMS_FILE)

# conventions/constants -----------------------------------

META_LOG_NAME = "metadata.pickle"

# the values of these don't matter as long as they are different
# example use: if thallic == HOMOTHALLIC: (action)
HOMOTHALLIC = True
HETEROTHALLIC = False

# Number of mutations "on-top-of" each other that we can model
TOP_POISSON = 5

# Core information ----------------------------------------
# change the value of -p based on what machine you are running on...
# PHOENIX: 2 core, 4 logical processors (2.40 GHz) -> 4
# Mycelium: 8 processors (2.66 GHz) -> 8

class TreeNode(object): # An object that represents a node in a phylogenetic tree
    # Store all metadata in this class & then pickle it.
    def __init__(self, genome, branch_name, state):
        self.genome = genome # genome of the node
        self.children = [] # list of children nodes (length 2 or 0)
        self.branch_name = branch_name # name of this node
        self.state = state # Homo or Hetero -thallic
        s = 0
        for meta, seq, bins in genome:
            s += seq.size
        self.size = s

        # self.base_nums = {}
        # # print "begin"
        # # for meta, seq in genome:
        # #     base, counts = np.unique(seq, return_counts=True)
        # #     self.base_nums[meta] = dict(zip(base, counts))
        # # print self.base_nums

        # Evolutionary parameters (valid for all but the common ancestor)
        self.t = 0
        self.rho = 0
        self.num_rearr = 0
        self.num_inv = 0

    def add_child(self, node):
        self.children.append(node)

    def pickle_me(self, some_file):
        # THIS WILL ERASE ALL GENOME DATA. 
        # DO NOT CALL THIS BEFORE YOU ARE SURE YOU HAVE THIS DATA WRITTEN TO FILE!
        self.delete_genomes()
        pickle.dump(self, some_file)

    # THIS WILL IRRECOVERABLY ERASE ALL GENOME DATA. ONLY CALL FROM PICKLER!
    def delete_genomes(self):
        self.genome = ""
        for tn in self.children:
            tn.delete_genomes()

    def depickle(some_file):
        return pickle.load(some_file)

    def write_all_leaves_to_files(self):
        if self.children:
            for child in self.children:
                child.write_all_leaves_to_files()
        else:
            # handle genome deletion? for the pickler?
            if not os.path.exists(OUT_DIRECTORY):
                os.makedirs(OUT_DIRECTORY)
            s = "O" if self.state == HOMOTHALLIC else "E"
            export_genome(self.genome, OUT_DIRECTORY + self.branch_name + s + ".fasta")

def main():
    # This script stores fasta files as triplists using the metadata as the first
    # value and the sequences as the second value (in numpy array form). All evolution-related 
    # operations are then performed on these triplists. 

    anc_genome = fasta_to_dict(IN_FILE)

    # Otherwise, we will run the evolution simulator.
    print "Common ancestor is " + IN_FILE.split('/')[-1]
    print "Reading common ancestor genome . . ."
    np_genome = []

    genome_len = sum([len(seq) for metadata, seq in anc_genome.items()])
    print "Genome has length", genome_len
    abs_pos = 0
    i_abs = 0
    sum_bin_lens = 0
    prev_remaining = 0
    prev_p = 0
    for i, (metadata, seq) in enumerate(anc_genome.items()): 
        sl_bins = [] if i == 0 else [Slyce(i, chrom_len, 0, prev_remaining, \
             1, prev_p)]
        i_chrom = 0 if i == 0 else 1
        pos_in_chrom = 0 if i == 0 else sl_bins[0].len()
        chrom_len = len(seq)
        while pos_in_chrom+BIN_INDICES[i_abs][1]-BIN_INDICES[i_abs][0] < chrom_len \
            and abs_pos+BIN_INDICES[i_abs][1]-BIN_INDICES[i_abs][0] < genome_len:
            this_bin = Slyce(i, chrom_len, pos_in_chrom, pos_in_chrom + BIN_INDICES[i_abs][1] - BIN_INDICES[i_abs][0], \
                1, BIN_INDICES[i_abs][2])
            bin_len = this_bin.len()

            sl_bins += [this_bin]
            i_chrom += 1
            i_abs += 1
            pos_in_chrom += bin_len
            abs_pos += bin_len
        final_bin = Slyce(i, chrom_len, pos_in_chrom, chrom_len, \
             1, BIN_INDICES[i_abs][2])
        prev_remaining = BIN_INDICES[i_abs][1] - BIN_INDICES[i_abs][0] - final_bin.len()
        prev_p = BIN_INDICES[i_abs][2]
        sl_bins += [final_bin]
        pos_in_chrom += final_bin.len()
        pos_in_chrom += final_bin.len()

        # print SlyceList(sl_bins).len(), len(seq)
        sum_bin_lens += SlyceList(sl_bins).len()
        # print SlyceList(sl_bins)
        # if i == 1:
        #     sys.exit(0)
        #np_genome += [[metadata, np.array(list(seq)), SlyceList([Slyce(i, len(seq))])]]
        np_genome += [[metadata, np.array(list(seq)), SlyceList(sl_bins)]]
    
    # print genome_len, sum_bin_lens

    # sys.exit()
    # for i, b in enumerate(bin_slyces):
    #     print "F", i, b
    # sys.exit(0)


    # Fork the genome recursively to generate a tree DEPTH deep, 
    # evolving it at each step...
    tree = fork_genome(TreeNode(np_genome, "1.", HETEROTHALLIC), "C.", DEPTH)
    # write .fasta files that represent the genomes of all the _leaves_. 
    # of the tree. We are not interested in the evolutionary intermediates. 
    tree.write_all_leaves_to_files()
    print "Pickling metadata to " + OUT_DIRECTORY + META_LOG_NAME + " . . ."
    tree.pickle_me(open(OUT_DIRECTORY + META_LOG_NAME, 'w'))
    print "Mission accomplished."

# fork the genome to a depth of DEPTH until we get 2^DEPTH terminal descendants
def fork_genome(node, branch_name, depth):
    if depth == DEPTH: # Makes sure we branch immediately (don't evolve the parent)
        node.add_child(fork_genome(node, branch_name + "1.", depth - 1))
        node.add_child(fork_genome(node, branch_name + "2.", depth - 1))
        return node
    # Goal: evolve node.genome into new_genome!
    n_a = node.size
    new_genome = []
    print "Evolving " + branch_name
    # Determine branch length t (by drawing from exponential
    # distribution with rate lambda). 
    t = random.expovariate(TIME_LAMBDA)
    print "t leading to node " + branch_name + " = " + str(t)
    # Am I homo/hetero-thallic?
    thallic = homo_or_hetero_thallic(node.state)
    if thallic == HOMOTHALLIC:
        print "Node " + branch_name + " is homothallic."
    else:
        print "Node " + branch_name + " is heterothallic."

    # Determine rearrangement rate rho (the freqs are the parameters we vary)
    rho = 0.0
    mean_rearr_size = 0.0
    sd_rearr_size = 0.0
    if thallic == HOMOTHALLIC: 
        rho = FREQ_REARR_HOMO
        mean_rearr_size = MEAN_INV_SIZE_HOMO
        sd_rearr_size = SD_INV_SIZE_HOMO
    else: # thallic == HETEROTHALLIC:
        rho = FREQ_REARR_HETERO
        mean_rearr_size = MEAN_INV_SIZE_HETERO
        sd_rearr_size = SD_INV_SIZE_HETERO

    # Find number of rearrangements
    # print rho 
    # print t 
    # print n_a
    rearrangements = int(m.floor(rho * t * n_a))
    print "This genome will undergo " + str(rearrangements) + " rearrangements. "

    # Point mutate the genome. lambda = mu * t
    if NUM_PARTITIONS == 1: # Use sub_mat and the old, 1-partition model for substitution
        new_genome = point_mut_single_partition(node.genome, MEAN_SUB_RATE * t, \
            MEAN_INS_RATE * t, MEAN_DEL_RATE * t)
    else: # Use the multipartition model
        new_genome = point_mut_multipartition(node.genome, t)
    """ Rearranger being debugged... """
    #new_genome = node.genome #TOBEREMOVED!!!


    #Rearrange the genome.
    rearr_genome = rearranger(new_genome, rearrangements, \
        mean_rearr_size, sd_rearr_size)
    #sys.exit(0)

    # After evolution finishes, create a new node for the new organism...
    new_node = TreeNode(rearr_genome, branch_name, thallic)
    new_node.t = t
    new_node.rho = rho
    new_node.num_rearr = rearrangements
    new_node.num_inv = 0 #FIXME

    if depth <= 0: # If we achieve the desired depth, STOP EVOLVING.
        return new_node
    # otherwise, evolve two children
    new_node.add_child(fork_genome(new_node, branch_name + "1.", depth - 1))
    new_node.add_child(fork_genome(new_node, branch_name + "2.", depth - 1))
    return new_node

# HOMO- OR HETERO- THALLIC DETERMINER ---------------------

def homo_or_hetero_thallic(current_state):
    p = random.uniform(0.0, 1.0)
    m = OUTPUT_HETERO
    if current_state == HOMOTHALLIC: 
        m = OUTPUT_HOMO
    if p < m.item(0): 
        return HOMOTHALLIC
    return HETEROTHALLIC

# POINT MUTATION GENERATOR --------------------------------

def single_partition_mutator(scaffold_triple):
    metadata, sequence, bins, lambda_sub, lambda_ins, lambda_del = scaffold_triple

    scaf_len = len(sequence)
    base_nums = cl.Counter(sequence)
    print "    Mutating scaffold " + metadata.split(' ')[0]
    # print bins.len(), len(sequence)
    poissons = [pmf_poisson(lambda_sub, i) for i in range(1,TOP_POISSON+1)]
    num_muts = {}
    coord_set_dict = {}
    for b in acgt:
        coord_set_dict[b] = set()
        num_muts[b] = []
        for ps in poissons:
            num_muts[b] += [np.random.binomial(base_nums[b], ps)]

    for b in acgt: # I am performing the 0th sequence of mutation
        b_to_x = set()          
        index = 0
        while index < num_muts[b][0]:
            rand_index = np.random.randint(0, scaf_len)
            if sequence[rand_index] == b:
                b_to_x.add(rand_index)
                index += 1
        probs = np.random.rand(num_muts[b][0])
        for prob_index, b_index in enumerate(b_to_x):
            coord_set_dict[BASES[np.searchsorted(CUM_OUT_MAP[b], \
                probs[prob_index])]].add(b_index)

    for p in range(1, TOP_POISSON): # I am performing the pth consecutive mutation
        for b in acgt: 
            if num_muts[b][p] == 0: # don't waste time/memory if I have no mutating to do
                continue
            b_to_x = set()
            for i in range(0, num_muts[b][p]):
                if not coord_set_dict[b]: # in the off chance I can't perform a mutation on 'b'; SHOULD BE VERY RARE
                    print 'Cannot mutate', b
                    continue
                    # print coord_set_dict
                    # print b, coord_set_dict[b], tuple(coord_set_dict[b])
                    # print "---"
                rand_element = random.choice(tuple(coord_set_dict[b]))
                b_to_x.add(rand_element)
                coord_set_dict[b].remove(rand_element)

            probs = np.random.rand(num_muts[b][p])
            for prob_index, b_index in enumerate(b_to_x):
                coord_set_dict[BASES[np.searchsorted(CUM_OUT_MAP[b], \
                    probs[prob_index])]].add(b_index)

    # deletions
    del_indexes = set()
    for b in acgt:
        del_lambda = lambda_del * OUTPUT_DEL[BASE_MAP[b]] * NUM_BASES
        del_num = np.random.binomial(base_nums[b], pmf_poisson(del_lambda, 1))
        index = 0
        while index < del_num:
            rand_index = np.random.randint(0, scaf_len)
            if sequence[rand_index] == b:
                del_indexes.add((rand_index, None))
                index += 1

    # insertions
    ins_indexes = set()
    for b in acgt:
        ins_lambda = lambda_ins * OUTPUT_INS[BASE_MAP[b]] * NUM_BASES
        ins_num = np.random.binomial(base_nums[b], pmf_poisson(ins_lambda, 1))
        for p in range(ins_num):
            ins_indexes.add((np.random.randint(0, scaf_len), b))

    # modify sequence according to coord_set_dict (point substitutions) here...
    for b in acgt:
        for index in coord_set_dict[b]:
            sequence[index] = b

    # for i in del_indexes: # delete stuff here...
    #     sequence[i] = ""
    #     si, tot_bases = bins.slyce_index_at_abs_pos(i)
    #     
    # for i in ins_indexes: # handle insertion tuples here...
    #     print "SI", sequence[i[0]]
    #     print i[1]
    #     print sequence[i[0]] + i[1]
    #     sequence[i[0]] = sequence[i[0]] + i[1]
    #     print "SQ", sequence[i[0]]
    #     sys.exit(0)
    #     si, tot_bases = bins.slyce_index_at_abs_pos(i)
    #     bins.sl[si].add_one()
    indels = sorted(ins_indexes | del_indexes)
    ##ins_indexes = list(sorted(ins_indexes, key=lambda inst: inst[0]))
    ##del_indexes = list(sorted(del_indexes))

    # # The following line is a huge runtime sink:
    # sequence = np.array(list("".join(sequence))) # recreate array

    ##ii = 0 # insertion list index
    ##di = 0 # deletion list index
    oi = 0 # index in new master sequence
    ni = 0 # index in old master sequence
    is_ins = False
    new_sequence = np.empty(shape=(scaf_len + len(ins_indexes) - len(del_indexes)), dtype=str)
    for s in bins.sl:
        s.length = len(new_sequence)


    for indel in indels:
        if indel[1]:
            # insertion
            pass
        else:
            # deletion
            pass
    # print len(ins_indexes), len(del_indexes)
    while di < len(del_indexes) and ii < len(ins_indexes):
        if di >= len(del_indexes) or ins_indexes[ii][0] < del_indexes[di]: # handle insertion
            #print "i.", ii, ins_indexes[ii]
            bins.sl[bins.slyce_index_at_abs_pos(ii)[0]].add_one()
            #print oi, ins_indexes[ii][0] , ni
            #print oi + ins_indexes[ii][0] - ni
            new_sequence[ni:ins_indexes[ii][0]] = sequence[oi:min(oi + ins_indexes[ii][0] - ni, len(new_sequence))]
            new_sequence[ins_indexes[ii][0]] = ins_indexes[ii][1]
            oi = oi + ins_indexes[ii][0] - ni
            ni = ins_indexes[ii][0] + 1
            ii += 1
            #print "i", ii, "d", di, "o", oi, "n", ni
        else: # handle deletion
            #print "d.", di, del_indexes[di]
            bins.sl[bins.slyce_index_at_abs_pos(di)[0]].zap_one()
            #print oi, del_indexes[di] , ni
            #print oi + del_indexes[di] - ni
            new_sequence[ni:del_indexes[di] - 1] = sequence[oi:oi + del_indexes[di] - 1 - ni]
            oi = oi + del_indexes[di] - 1 - ni
            ni = del_indexes[di]
            di += 1
            #print "i", ii, "d", di, "o", oi, "n", ni
    ti = max(oi, ni) # top index
    new_sequence[ti:len(new_sequence)] = sequence[oi:oi+len(new_sequence)-ti]

    # print new_sequence, len(new_sequence)

    # if len(new_sequence) != bins.len():
    #     print "ERROR: Misaligned bins.", len(new_sequence), "!=", bins.len()
    #     sys.exit(0)

    # the following are temporary hacks to fix the bin alignment bug in the above code:
    new_sequence = np.array([random.choice(BASES) if x != '' else x for x in new_sequence])
    
    if len(new_sequence) > bins.len():
        new_sequence = new_sequence[0:bins.len()]
    elif len(new_sequence) < bins.len():
        new_sequence = np.array(list(new_sequence) + [random.choice(BASES) for i in range(bins.len()-len(new_sequence))])

    # print bins.len(), len(new_sequence)
    return metadata, new_sequence, bins

# def multipartition_mutator(scaffold_triple):
#     metadata, sequence, lambda_sub, lambda_ins, lambda_del = scaffold_triple
#     return metadata, sequence

def pmf_poisson(lambda_m, k):
    return m.exp(-lambda_m) * lambda_m**k / m.factorial(k)

# simulates point mutations (base changes, insertions, deletions) in a genome triplist
def point_mut_single_partition(genome_triplist, lambda_sub, lambda_ins, lambda_del): 
    print "Generating point mutations using single partition model."
    # Use mutator(scaff) on all scaff in scaffolds
    pool = Pool(NUM_PROCESSORS)
    # change following to pool.map()
    mapped_vals = map(single_partition_mutator, \
        [[p[0], p[1], p[2], float(lambda_sub), float(lambda_ins), float(lambda_del)] \
        for p in genome_triplist])
    #   sys.exit(0)
    return mapped_vals

# # same as above, but uses the multipartition model
# def point_mut_multipartition(genome_triplist, evolutionary_time):
#     # multipartition MapReduce
#     print "Generating point mutations using " + str(NUM_PARTITIONS) + "-partition model."
#     # pool = Pool(NUM_PROCESSORS)


#     mapped_vals = map(multipartition_mutator, \
#         [[p[0], p[1], float(lambda_sub), float(lambda_ins), float(lambda_del)] \
#         for p in genome_triplist])
#     return mapped_vals

# CHROMOSOME REARRANGER -----------------------------------
def rearranger(genome_triplist, num_rearr, mean_rearr_size, sd_rearr_size):
    #TODO need to handle (rare) negative values (reselect from distribution)
    # these values can technically be > scaffold size (WE ASSUME THIS IS NOT THE CASE)
    # ^^ this is ludicrously unlikely; due to performance overhead, don't test for this
    rearr_sizes = np.round(np.random.normal(mean_rearr_size, sd_rearr_size, num_rearr)) \
        .astype(int)

    # print genome_triplist
    metadata = [m for m, s, b in genome_triplist]
    seq = [s for m, s, b in genome_triplist]
    bins = [b for m, s, b in genome_triplist]

    # Nx3 array of random numbers in (0, 1] used in the rearrangement model
    io_determiners = np.random.random_sample((num_rearr,3))
    # lengths of each chromosome
    lens = [len(s) for s in seq]
    sum_lens = sum(lens)
    num_scaff = len(lens)

    # possible start locations
    # sum_poss_start_locs = np.add(np.multiply(rearr_sizes, -num_scaff), num_scaff+sum_lens)

    # print rearr_sizes

    # chrom_slyce_lists = [] # List of slyce lists for each chromosome
    # chrom_partition_list = []
    # for i, (metadata, seq, bins) in enumerate(genome_triplist): 
    #     chrom_slyce_lists += [SlyceList([Slyce(i, lens[i])])]

    chrom_slyce_lists = bins

    #print lens
    x = sum_lens

    for i in range(num_rearr): # perform the ith rearrangement
        print "Performing rearrangement #" + str(i), "--------------------------------------------"
        size_rearr = rearr_sizes[i] # grab the ith rearrangement size
        while size_rearr < 0: # make sure it's not negative! (possible under conditions of normal distribution)
            size_rearr = int(np.random.normal(mean_rearr_size, sd_rearr_size, num_rearr))
        sum_poss_start_loc = sum_lens - num_scaff * size_rearr

        # # Single partition uniform uses the following locator: 
        # abs_pos = np.random.randint(0, sum_poss_start_loc) # absolute position of start loc of rearr

        # N-bin partition rearranger uses the following locator: 
        print np.searchsorted(np.cumsum(PARTITION_PROBS), io_determiners[i][2], side='left')
        partition = np.searchsorted(np.cumsum(PARTITION_PROBS), io_determiners[i][2])
        partition_len = sum([chrom_slyce_lists[k].len_partition(partition) for k in range(num_scaff)])

        print io_determiners[i][2], PARTITION_PROBS, partition
        print num_scaff
        print [chrom_slyce_lists[i].len_partition(partition) for i in range(num_scaff)]
        print chrom_slyce_lists[0].len_partition(1)
        print chrom_slyce_lists[0].len()

        print "partition_len", partition_len
        partition_pos = np.random.randint(0, partition_len) # absolute position of start loc of rearr WITHIN THE PARTITION
        my_chrom = 0
        while chrom_slyce_lists[my_chrom].len_partition(partition) - size_rearr + 1 < partition_pos:
            partition_pos -= (chrom_slyce_lists[my_chrom].len_partition(partition) - size_rearr + 1)
            my_chrom += 1

        abs_pos = chrom_slyce_lists[my_chrom].partition_pos_to_abs_pos(partition, partition_pos)

        #print "Chrom", my_chrom, "looks like", chrom_slyce_lists[my_chrom]
        print "From chrom", str(my_chrom) + ", length=(" + str(lens[my_chrom]) + ") at position", abs_pos, "grab length", size_rearr

        insert_sl = None
        if io_determiners[i][0] < LIKELIHOOD_INV: # is it an in-place inversion?
            print "Performing in-place inversion."
            new_parent, insert_sl = chrom_slyce_lists[my_chrom].excise(abs_pos, size_rearr)
            chrom_slyce_lists[my_chrom] = new_parent
            insert_sl = insert_sl.invert()
            chrom_slyce_lists[my_chrom] = chrom_slyce_lists[my_chrom].insert(abs_pos, insert_sl)
            #print "Chrom", my_chrom, "now looks like", chrom_slyce_lists[my_chrom]

        elif io_determiners[i][0] < LIKELIHOOD_INV*2: # is it an excision of some sort?
            new_parent, insert_sl = chrom_slyce_lists[my_chrom].excise(abs_pos, size_rearr)
            chrom_slyce_lists[my_chrom] = new_parent
            lens[my_chrom] -= size_rearr
            sum_lens  -= size_rearr
            print "Excising segment:", insert_sl
            #print "Parent SlyceList is now:", new_parent
            if io_determiners[i][0] < LIKELIHOOD_INV*2 - LIKELIHOOD_IO_REARR: # is it a translocation excision?
                print "This is a translocation."
                if random.getrandbits(1): 
                    insert_sl = insert_sl.invert()
                    print "Inverting segment into", insert_sl

                ins_pos = np.random.randint(0, sum_lens) # absolute position of start loc of rearr
                ins_chrom = 0 # chromosome where we are inserting
                # print abs_pos
                while lens[ins_chrom] <= ins_pos:
                    ins_pos -= lens[ins_chrom]
                    ins_chrom += 1

                # print type(insert_sl), insert_sl
                print "Translocating segment into chrom", ins_chrom, "at pos", ins_pos
                chrom_slyce_lists[ins_chrom] = chrom_slyce_lists[ins_chrom].insert(ins_pos, insert_sl)
                #print "Chrom", ins_chrom, "now looks like", chrom_slyce_lists[ins_chrom]

                sum_lens += size_rearr
                lens[ins_chrom] += size_rearr
            else: # is it a deletion?
                print "This is a deletion. This segment will now be junked."

        elif io_determiners[i][0] < LIKELIHOOD_INV*2 + LIKELIHOOD_IO_REARR: # is it a duplication (extraction & pasting) event?
            insert_sl = chrom_slyce_lists[my_chrom].extract(abs_pos, size_rearr)
            print "Extracting segment:", insert_sl
            if io_determiners[i][1] < LIKELIHOOD_TANDEM: # is it a tandem (postfix) duplication
                print "This is a tandem duplication."
                insert_sl = chrom_slyce_lists[my_chrom].extract(abs_pos, size_rearr)

                chrom_slyce_lists[my_chrom] = chrom_slyce_lists[my_chrom].insert(abs_pos, insert_sl)
                #print "Chrom", my_chrom, "now looks like", chrom_slyce_lists[my_chrom]
            else: # is it a non-tandem, randomly inserted duplication?
                print "This is a non-tandem duplication."
                if random.getrandbits(1): 
                    insert_sl = insert_sl.invert()
                    print "Inverting segment into", insert_sl

                ins_pos = np.random.randint(0, sum_lens) # absolute position of start loc of rearr
                ins_chrom = 0 # chromosome where we are inserting
                # print abs_pos
                while lens[ins_chrom] <= ins_pos:
                    ins_pos -= lens[ins_chrom]
                    ins_chrom += 1

                # print type(insert_sl), insert_sl
                print "Duplicating segment into chrom", ins_chrom, "at pos", ins_pos
                chrom_slyce_lists[ins_chrom] = chrom_slyce_lists[ins_chrom].insert(ins_pos, insert_sl)
                #print "Chrom", ins_chrom, "now looks like", chrom_slyce_lists[ins_chrom]
                lens[ins_chrom] += size_rearr
                sum_lens += size_rearr
        else: # is it an assisted translocation (transposon, etc...)?
            print "This is an assisted translocation."
            new_parent, insert_sl = chrom_slyce_lists[my_chrom].excise(abs_pos, size_rearr)
            chrom_slyce_lists[my_chrom] = new_parent
            lens[my_chrom] -= size_rearr

            print "Excising segment:", insert_sl
            #print "Parent SlyceList is now:", new_parent
            if random.getrandbits(1): 
                insert_sl = insert_sl.invert()
                print "Inverting segment into", insert_sl

            ins_pos = np.random.randint(0, sum_lens) # absolute position of start loc of rearr
            ins_chrom = 0 # chromosome where we are inserting
            # print abs_pos
            while lens[ins_chrom] <= ins_pos:
                ins_pos -= lens[ins_chrom]
                ins_chrom += 1

            # print type(insert_sl), insert_sl
            print "Translocating segment into chrom", ins_chrom, "at pos", ins_pos
            chrom_slyce_lists[ins_chrom] = chrom_slyce_lists[ins_chrom].insert(ins_pos, insert_sl)
            #print "Chrom", ins_chrom, "now looks like", chrom_slyce_lists[ins_chrom]
            lens[ins_chrom] += size_rearr

    # print "x:", x, sum_lens

    # f = open('slycelist'+str(num_rearr) + '.txt', 'w')
    # for i, csl in enumerate(chrom_slyce_lists):
    #     #print "L", i, csl.len()
    #     f.write(str(i))
    #     f.write(str(csl))
    # f.close()

    #Reconstructing original chromosomes according to SlyceLists
    new_genome_triplist = []
    for i, (metadata, seq, old_bins) in enumerate(genome_triplist): 
        new_seq = np.empty(chrom_slyce_lists[i].len(), dtype='string')         #TODO bug due to faulty len tracking (need to find!!) squashed???
        # print 'csl', chrom_slyce_lists[i].len(), 

        loc = 0
        for s in chrom_slyce_lists[i].sl:
            # print new_seq[loc:loc+s.len()]
            # print genome_triplist[s.l][1][s.i1:s.i2:s.i3]

            # print loc, s.len(), s.l, s.i1, s.i2, s.i3, len(genome_triplist[s.l][1]), len(new_seq)
            if len(genome_triplist[s.l][1][s.i1:s.i2:s.i3]) < len(new_seq[loc:loc+s.len()]):
                new_seq[loc:loc+len(genome_triplist[s.l][1][s.i1:s.i2:s.i3])] = genome_triplist[s.l][1][s.i1:s.i2:s.i3]
            elif len(genome_triplist[s.l][1][s.i1:s.i2:s.i3]) > len(new_seq[loc:loc+s.len()]):
                new_seq[loc:loc+s.len()] = genome_triplist[s.l][1][s.i1:s.i2:s.i3][0:s.len()]
            else:
                new_seq[loc:loc+s.len()] = genome_triplist[s.l][1][s.i1:s.i2:s.i3]
            loc += s.len()
        new_genome_triplist += [(metadata, new_seq, chrom_slyce_lists[i])]

    # print list(enumerate(genome_triplist))
    # for i, (old_m, old_seq) in enumerate(genome_triplist): 
    #     new_seq = new_genome_triplist[i][1]
    #     print i
    #     for j, b in enumerate(list(new_seq)):
    #         if b != new_seq[i]:
    #             print j, ": ", b, '!=', new_seq
        
    # sys.exit()
    return new_genome_triplist

if __name__ == '__main__':
    main()
