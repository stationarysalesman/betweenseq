"""betweenseq: A tool to find repeat sequences that can mediate deletions in DNA sequences"""

from Bio import SeqIO
import os
import argparse


def get_runs(indices, distances, max_dist):
    """Recursively process indices list to return a list
    of frames in which a repeat sequence could potentially 
    mediate a deletion.
    
    precondition: indices is a sorted list"""
    
    start = 0
    x = 0 # index
    dist_len = len(distances)

    # base cases
    if dist_len == 1:
        if distances[0] > max_dist:
            return None
        else:
            return indices
    if not (indices or distances):
        return

    # recursive case
    while (x < dist_len and distances[x] <= max_dist):
        x += 1

    # done
    if (x == dist_len):
        return indices
    
    if start == x: # No good runs
        return get_runs(indices[1:], distances[1:], max_dist)
    else:
        if x == 1: # only one good, need to check return before building list
            item = indices[0:2]
            recurse = get_runs(indices[x:], distances[x:], max_dist)
            if not recurse:
                return
            return [item, recurse]

        else:
            item = indices[0:x]
            recurse = get_runs(indices[x:], distances[x:], max_dist)
            if not recurse:
                return item
            return [item, recurse]
    
def process_index_list(index_list):
    """Remove repeats from index_list that are too far apart.

    Do this by building a new list that stores the differences
    between each index in the list. Example:

    index_list = [1, 7, 10, 24, 456]

    new_list = [6, 3, 14, 432]

    If an element in new_list is greater than our maximum 
    allowable distance, kick out the previous index (so that 
    if the current one is within the distance of the next one,
    we don't kick it out prematurely."""

    max_dist = 100 # maximum nucleotide length allowed between repeat sequences
    new_list = list() # build a new list with distances between indices
    x = 1 # start index at 2nd element and compare to previous
    
    while (x < len(index_list)):
        new_list.append(index_list[x] - index_list[x-1])
        x += 1

    # helper function
    return get_runs(index_list, new_list, max_dist)
            
def verify_indices(repeats):
    """Verify that any repeat sequences are not too far apart."""
    for k in repeats.keys():
        index_list = repeats[k]
        # remove indices that are too far apart
        repeats[k] = process_index_list(index_list)
    return repeats

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
        seq_frag = str(seq[x:x+frame_size])

        if seq_frag in repeats: # repeated sequence element
            repeats[seq_frag].append(x)
        else:
            repeats[seq_frag] = list()
            repeats[seq_frag].append(x)
        x += 1
    # Remove tandem repeats and simple sequence repeats
    # note: we can change this later if we want all of this data
    # repeats = process_repeats(repeats, frame_size)

    # verification
    repeats = verify_indices(repeats)
    return repeats

def get_repeats(seq_file, frame_min, frame_max):
    """Find repeat sequences for all frame sizes between frame_min and frame_max
    and store them in a dictionary."""

    # Get sequence file metadata
    seq_id = seq_file.id
    seq_seq = seq_file.seq

    # Dictionary that will store repeats
    repeats = dict()
    for x in range(frame_min, frame_max+1): # +1 to include max frame size
        repeats = get_frame_repeats(seq_seq, x)

        # TODO: format output
        if not repeats:
            return
        print "Repeats in sequence", seq_id + " (frame_size = "+str(x) + "):"
        for k in repeats.keys():
            lst = repeats[k]
            if (lst and len(lst) > 1):
                print str(k) + ": " + str(repeats[k])
        
    return
    
def main():
    """Controller for analysis work flow. """

    # Define error strings
    ERR_NO_DIR = "ERROR: directory not accessible:"

    # Miscellaneous program varibles
    FRAME_MIN = 5
    FRAME_MAX = 10
    
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str, help="file containing sequences to analyze")
    #parser.add_argument('--group', type=str, help="analyze a group of sequences")
    #parser.add_argument('--format', type=str, help="specify sequence file format")
    args = parser.parse_args()
    seq_dir = args.file
    #group = args.group
    #file_format = args.format
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
