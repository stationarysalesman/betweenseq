"""betweenseq: A tool to find repeat sequences that can mediate deletions in DNA sequences"""

from Bio import SeqIO
import os
import argparse


def remove_tandems(index_list):
    i = 0
    # Find tandem repeats by removing adjacent repeat frames

    # TODO: implement this algorithm or check before storing repeats
    while (i < len(index_list) - 1):
        while (index_list[i]  == index_list[i+1]):
            pass

def process_repeats(repeats, frame_size):
    """Remove tandem repeats and simple sequence repeats"""

    for repeat_seq in repeats:
        index_list = repeats[repeat_seq]
        new_list = remove_tandems(index_list)
        
def get_frame_repeats(seq, frame_size):
    """Store all indices of repeated sequences in a dictionary,
    with the sequence fragment as the key and a list of indices as the value."""

    repeats = dict() # Dictionary to store our repeats
    x = 0 # index into array of nucleotides in sequence

    while (x < len(seq) - frame_size):
        # Hash the current frame and check for collisions
        seq_frag = seq[x:x+frame_size]
        if seq_frag in repeats: # repeated sequence element
            repeats[seq_frag].append(x)
        else:
            repeats[seq_frag] = list()
            repeats[seq_frag].append(x)

    # Remove tandem repeats and simple sequence repeats
    # note: we can change this later if we want all of this data
    repeats = process_repeats(repeats, frame_size)
    
    return repeats

def get_repeats(seq_file, frame_min, frame_max):
    """Find repeat sequences for all frame sizes between frame_min and frame_max
    and store them in a dictionary."""

    # Get sequence file metadata
    seq_id = seq_file.id
    seq_seq = seq.seq

    # Dictionary that will store repeats
    repeats = dict()
    for x in range(frame_min, frame_max+1): # +1 to include max frame size
        repeats = get_frame_repeats(seq.seq, x)
        
        

    
def main():
    """Controller for analysis work flow. """

    # Define error strings
    ERR_NO_DIR = "ERROR: directory not accessible:"

    # Miscellaneous program varibles
    FRAME_MIN = 5
    FRAME_MAX = 10
    
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=str, help="directory containing sequences to analyze")
    args = parser.parse_args()
    seq_dir = args.directory
    if not(os.access(seq_dir, os.EX_OK)):
        print ERR_NO_DIR, seq_dir
        return -1

    seq_list = list()
    for dirName, subdirList, fileList in os.walk(seq_dir):
        for f in fileList:
            seq_list.append(seq_dir+f)
        
    
    for seq in seq_list:
        seq_file = SeqIO.read(seq, "genbank")
        get_repeats(seq_file, FRAME_MIN, FRAME_MAX)
        
    return

main()
