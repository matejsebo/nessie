#!/usr/bin/env python
import subprocess as sp

# def seq_dump(fastafile, begin, end): 
#     """Opens and reads a fasta formatted genome and returns the sequence corresponding to the inputted Murasaki coords."""
#     from Bio import SeqIO
#     records = SeqIO.parse(open(fastafile, 'r'), 'fasta')
#     full_genome = list(records)
#     ''.join(full_genome)

def get_absolute_positions(fastafile):
    """Opens and reads a fasta formated genome and returns a dictionary suitable for use in parse_mummer_coords()."""
    from Bio import SeqIO
    records = SeqIO.parse(open(fastafile, 'r'), 'fasta')
    d = {}
    genomelength = 0
    for rec in records:
        d[rec.name] = genomelength
        genomelength += len(rec.seq)
    d['end'] = genomelength
    return d

def get_scaffold_lengths(fastafile):
    """Opens and reads a fasta formated genome and returns a dictionary with sequence names as keys and the length of the sequence as values."""
    from Bio import SeqIO
    records = SeqIO.parse(open(fastafile, 'r'), 'fasta')
    d = SeqIO.to_dict(records)
    d = {key: len(rec.seq) for key, rec in d.items()}
    return d

def seq_length_list(fastafile):
    """Returns a list of the lengths of the sequences in a fasta file in order."""
    data = open(fastafile, 'r').read()

    seqs = data[1:].split('>')
    pairs = [seq.split('\n', 1) for seq in seqs]
    lengths = [len(pair[1].replace('\n', '')) for pair in pairs]

    return lengths

def total_length(fastafile):
    grep = sp.Popen(['grep', '-v', '">"', fastafile], stdout=sp.PIPE)
    wc = sp.Popen(['wc'], stdin=grep.stdout, stdout=sp.PIPE)
    out, err = wc.communicate()
    lines, words, characters = map(int, out.split())
    return characters - lines

def fasta_to_dict(fastafile):
    data = open(fastafile, 'r').read()

    genes = data[1:].split('>')
    genes = [gene.split('\n', 1) for gene in genes]
    try:
        genes = {gene[0] : gene[1].replace('\n', '') for gene in genes}
    except IndexError:
        print(fastafile)
        raise

    return genes
