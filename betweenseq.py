"""betweenseq: A tool to find repeat sequences that can mediate deletions in DNA sequences"""

from Bio import SeqIO
import os
import argparse
from repeatsequences import RepeatSequences


class RepeatParameters:

    def __init__(self, seq_path, file_format, maxdist, group, maxframe,
                 minframe):
        self.seq_path = seq_path
        self.file_format = file_format
        self.maxdist = maxdist
        self.group = group
        self.maxframe = maxframe
        self.minframe = minframe
        return
    
        
def is_ssr(seq_str):
    """Return true if seq_str is a string of repeated nucleotides."""
    x = 0
    while (x<len(seq_str)-1):
        if (seq_str[x] != (seq_str[x+1])):
            return False
        x += 1
    return True
        
def get_runs(indices, distances, maxdist):
    """Recursively process indices list to return a list
    of frames in which a repeat sequence could potentially 
    mediate a deletion.
    
    precondition: indices is a sorted list"""
    
    start = 0
    x = 0 # index
    dist_len = len(distances)

    # base cases
    if dist_len == 1:
        if distances[0] > maxdist:
            return None
        else:
            return indices
    if not (indices or distances):
        return

    # recursive case
    while (x < dist_len and distances[x] <= maxdist):
        x += 1

    # done
    if (x == dist_len):
        return indices
    
    if start == x: # No good runs
        return get_runs(indices[1:], distances[1:], maxdist)
    else:
        if x == 1: # only one good, need to check return before building list
            item = indices[0:2]
            recurse = get_runs(indices[x:], distances[x:], maxdist)
            if not recurse:
                return
            return [item, recurse]

        else:
            item = indices[0:x]
            recurse = get_runs(indices[x:], distances[x:], maxdist)
            if not recurse:
                return item
            return [item, recurse]
    
def process_index_list(index_list, maxdist):
    """Remove repeats from index_list that are too far apart.

    Do this by building a new list that stores the differences
    between each index in the list. Example:

    index_list = [1, 7, 10, 24, 456]

    new_list = [6, 3, 14, 432]

    If an element in new_list is greater than our maximum 
    allowable distance, kick out the previous index (so that 
    if the current one is within the distance of the next one,
    we don't kick it out prematurely."""

    new_list = list() # build a new list with distances between indices
    x = 1 # start index at 2nd element and compare to previous
    
    while (x < len(index_list)):
        new_list.append(index_list[x] - index_list[x-1])
        x += 1

    # helper function
    return get_runs(index_list, new_list, maxdist)
            
def verify_indices(repeats, params):
    """Verify that any repeat sequences are not too far apart."""
    for k in repeats.keys():
        if (is_ssr(k)):
            del(repeats[k])
            continue
        index_list = repeats[k]
        if not index_list:
            continue
        # remove indices that are too far apart
        repeats[k] = process_index_list(index_list, params.maxdist)
    return repeats


def check_bitmap(bitmap, seq_frag, index):
    """Return whether the sequence fragment at location index
    is covered by another repeat sequence."""

    
    y = len(seq_frag)
    index_range_list = range(index, index+y+1) # include last index
    if (index+y > bitmap.size()):
        return False

    # Check a range of values in the bitmap
    return bitmap.test_range(index_range_list)
        

def update_bitmap(bitmap, lst, frame_size):
    """Set a range of values equal to frame_size 
    in a bitmap according to indices in lst and 
    frame_size."""
    for index in lst:
        bitmap.set_range(range(index, index+frame_size+1)) # include last index

def add_subseq(BetweenSeqMaps, seq_frag, index, minframe):
    """Add all subsequences of seq_frag to a subsequences dictionary."""

    subseq = BetweenSeqMaps.get_subseq_dict()
    maxframe = len(seq_frag)
    for frame_size in range(maxframe-1, minframe-1, -1): # include min size
        idx = 0
        while (idx+frame_size < len(seq_frag)+1):
            seq = seq_frag[idx:idx+frame_size]
            if not seq in subseq:
                subseq[seq] = list()
            if not (index+idx) in subseq[seq]:
                subseq[seq].append(index+idx)
            idx += 1
    return


def get_frame_repeats(BetweenSeqMaps, params, seq, frame_size):
    """Store all indices of repeated sequences in a dictionary,
    with the sequence fragment as the key and a list of indices as the value."""

    repeats = BetweenSeqMaps.get_repeat_dict() # Dictionary to store our repeats
    frame_min = BetweenSeqMaps.get_frame_min()
    frame_max = BetweenSeqMaps.get_frame_max()
    bitmap = BetweenSeqMaps.get_bitmap()
    
    x = 0 # index into array of nucleotides in sequence

    while (x < len(seq) - frame_size):
        # Hash the current frame and check for collisions
        seq_frag = str(seq[x:x+frame_size])
        if check_bitmap(bitmap, seq_frag, x): # already something there
            add_subseq(BetweenSeqMaps, seq_frag, x, params.minframe)
            x += 1
            continue
        
        if seq_frag in repeats: # repeated sequence element
            repeats[seq_frag].append(x)
            # Update bitmap
            update_bitmap(bitmap, repeats[seq_frag], len(seq_frag))
            
        else:
            repeats[seq_frag] = list()
            repeats[seq_frag].append(x)
            
        x += 1

    # verification
    if params.maxdist != -1: # Default value means no maximum distance
        repeats = verify_indices(repeats, params)
    return repeats

def get_repeats(seq_file, params):
    """Find repeat sequences for all frame sizes between frame_min and frame_max
    and store them in a dictionary."""

    # Get sequence file metadata
    seq_id = seq_file.id
    seq_seq = seq_file.seq
    # Unpack arguments
    maxframe = params.maxframe
    minframe = params.minframe

    BetweenSeqMaps = RepeatSequences(len(seq_seq), minframe, maxframe)
    
    for x in range(maxframe, minframe-1, -1): # +1 to include max frame size
        get_frame_repeats(BetweenSeqMaps, params, seq_seq, x)
        repeats = BetweenSeqMaps.get_repeat_dict()
        # TODO: format output
        if not repeats:
            return

    # Sort the subsequence dictionary
    subseqs = BetweenSeqMaps.get_subseq_dict()
    for k in subseqs:
        subseqs[k] = sorted(subseqs[k])
        
    print "Repeats in sequence", seq_id + ":"
        
    for k in repeats.keys():
        lst = repeats[k]
        if (lst and len(lst) > 1):
            print str(k) + ": " + str(repeats[k])

    print "----------------------------"
    print "Subsequences:"
    for k in subseqs.keys():
        lst = subseqs[k]
        if (lst and len(lst) > 1):
            print str(k) + ": " + str(subseqs[k])

        
    return


def process_files(params):
    """Process a group of sequence files."""

    # Unpack parameters
    seq_path = params.seq_path
    file_format = params.file_format
    maxframe = params.maxframe
    minframe = params.minframe
    group = params.group

    if group:
        seq_lst = list()
        for dirName, subdirList, fileList in os.walk(seq_path):
            seq_lst = fileList
        for seq in seq_lst:
            seq_file = SeqIO.read(seq_path+seq, file_format)
            get_repeats(seq_file, maxframe, minframe)

    else:
        seq_file = SeqIO.read(seq_path, file_format)
        get_repeats(seq_file, params)
        
def main():
    """Controller for analysis work flow. """

    # Define error strings
    ERR_NO_FILE = "ERROR: file not accessible:"

    
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str,
                        help="file containing sequences to analyze")
    parser.add_argument('--format', type=str, default="genbank",
                        help="specify sequence file format")
    parser.add_argument('--maxdist', type=int, default=-1,
                        help="maximum allowable distance between repeats")
    parser.add_argument('-group', help="analyze a group of sequences",
                        action='store_true')
    parser.add_argument('--maxframe', type=int, default=10,
                        help="maximum size of repeats")
    parser.add_argument('--minframe', type=int, default=5,
                        help="minimum size of repeats")
    args = parser.parse_args()

    # Encapsulate parameters in object
    params = RepeatParameters(args.file, args.format, args.maxdist,
                              args.group, args.maxframe, args.minframe)
    # Error checking
    if not(os.access(params.seq_path, os.R_OK)):
        print ERR_NO_FILE, params.seq_path
        return -1

    #Run analysis
    process_files(params)

    return 0

main()
