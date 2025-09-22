from Executor import Executor
from Similarity import *

# example for one pair of sequences
def execute_for_one_pair():
    k_list = [4, 5]  # k-mer list
    r_list = [0, 2] # markove order list
    windows = 800 # windows size
    shift = 800 # shift size
    function_list = [ChSimilarity(), EuSimilarity(), D2Similarity(), MaSimilarity(), D2SSimilarity(), D2StarSimilarity()] # similarity function list
    seq_1 = r"D:\test\seqY\seq2.txt"  # sequence 1 file path
    seq_2 = r"D:\test\seqY\seq3.txt"  # sequence 2 file path
    save_dir = "." # dir to save result
    Executor(k_list, r_list, windows, shift, function_list, seq_1, seq_2, save_dir).execute()


if __name__ == '__main__':
    execute_for_one_pair()