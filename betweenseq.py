"""betweenseq: A tool to find repeat sequences that can mediate deletions in DNA sequences"""

from Bio import SeqIO
import os
import argparse
from repeatsequences import RepeatSequences


def is_tandem(seq_str):
    """Return true if seq_str is a string of repeated nucleotides."""
    x = 0
    while (x<len(seq_str)-1):
        if (seq_str[x] != (seq_str[x+1])):
            return False
        x += 1
    return True
        
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

    max_dist =  60 # maximum nucleotide length allowed between repeat sequences
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
        if (is_tandem(k)):
            del(repeats[k])
            continue
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

    return


def list_intersection(a_, b_):
    """Return the set intersection of two sets obtained from two lists."""
    set_a = set(a_)
    set_b = set(b_)
    return set_a & set_b
    
def get_subseq_by_frame_size(subseq_frame_size, seq_frag, subsequences, index_list):
    """Helper for get_frame_repeats. Determine if subsequences of 
    size subseq_frame_size have been found at index x in dictionary
    subsequences."""

    ret = False
    
    for x in range(len(seq_frag)):
        subseq = seq_frag[x:x+subseq_frame_size]
        # Check if the subsequence has been seen, and occurred at index
        if subseq in subsequences:
            subseq_indices = subsequences[subseq]
            if list_intersection(subseq_indices, index_list):
                ret = True
            else:
                subsequences[subseq].append(x)
        else:
            subsequences[subseq] = list()
            subsequences[subseq].append(x)
        x += 1

    return ret


def find_sub_repeats(seq_frag, frame_min, subsequences, index_list):
    """Determine if any subsequences in seq_frag have been encountered
    around this index.

    index: """

    x = 0
    ret = False
    start_size = len(seq_frag)-1
    for subseq_frame_size in range(start_size, frame_min, -1):
        ret = get_subseq_by_frame_size(subseq_frame_size, seq_frag, subsequences, index_list) or ret
        
    return ret
            


def check_bitmap(bitmap, seq_frag, index):
    """Return whether the sequence fragment at location index
    is covered by another repeat sequence."""

    
    y = len(seq_frag)
    index_range_list = range(index, len(seq_frag)+1) # inlude last index
    
    
    
    if (index+y > bitmap.size()):
        return False
    
    for x in range(index, index+y+1):
        

    
def get_frame_repeats(BetweenSeqMaps, seq, frame_size):
    """Store all indices of repeated sequences in a dictionary,
    with the sequence fragment as the key and a list of indices as the value."""

    repeats = BetweenSeqMaps.get_repeat_dict() # Dictionary to store our repeats
    subsequences = BetweenSeqMaps.get_subsequences_dict() # Dictionary to store sub-repeats
    frame_min = BetweenSeqMaps.get_frame_min()
    frame_max = BetweenSeqMaps.get_frame_max()
    bitmap = BetweenSeqMaps.get_bitmap()
    
    x = 0 # index into array of nucleotides in sequence

    while (x < len(seq) - frame_size):
        # Hash the current frame and check for collisions
        seq_frag = str(seq[x:x+frame_size])
        check = check_bitmap(bitmap, seq_frag, x)
        if seq_frag in repeats: # repeated sequence element
            repeats[seq_frag].append(x)
            update_bitmap(
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

    BetweenSeqMaps = RepeatSequences(len(seq_seq), frame_min, frame_max)
    
    for x in range(frame_max+1, frame_min, -1): # +1 to include max frame size
        get_frame_repeats(BetweenSeqMaps, seq_seq, x)

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
    ERR_NO_FILE = "ERROR: file not accessible:"

    # Miscellaneous program varibles
    FRAME_MIN = 6
    FRAME_MAX = 8
    file_format = "genbank" # default
    
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str, help="file containing sequences to analyze")
    # parser.add_argument('--group', type=str, help="analyze a group of sequences")
    parser.add_argument('--format', type=str, help="specify sequence file format")
    args = parser.parse_args()
    seq_file = args.file
    # group = args.group
    if args.format:
        file_format = args.format
    if not(os.access(seq_file, os.R_OK)):
        print ERR_NO_FILE, seq_file
        return -1
    
    seq = SeqIO.read(seq_file, file_format)
    get_repeats(seq, FRAME_MIN, FRAME_MAX)

    return 0

main()
