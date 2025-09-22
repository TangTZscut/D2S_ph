from typing import List, Tuple
from unittest.mock import Base
import numpy as np

from .BaseSimilarity import BaseSimilarity
from Data import Sequence

class D2SSimilarity(BaseSimilarity):
    
    def calculate(self, seq_1: Sequence, seq_2: Sequence, k: int, r: int = -1) -> float:
        if k <= r:
            return np.nan
        x_kmer = seq_1.get_kmer(k)
        y_kmer = seq_2.get_kmer(k)
        x_markov = seq_1.get_markov(k, r)
        y_markov = seq_2.get_markov(k, r)
        x = x_kmer.data
        y = y_kmer.data
        px = x_markov.probablity
        py = y_markov.probablity
        total_x = x_kmer.total
        total_y = y_kmer.total
        tmp_x = px * total_x
        tmp_y = py * total_y
    
        x_bar = x - tmp_x
        y_bar = y - tmp_y
        tmp_1 = x_bar * x_bar
        tmp_2 = y_bar * y_bar
        
        tmp_D2S = np.sqrt(tmp_1 + tmp_2)
        tmp_D2S = np.where(tmp_D2S == 0, 1, tmp_D2S)
        res = x_bar * y_bar / tmp_D2S
        res = np.sum(res)
        tmp_x = np.sqrt(np.sum(tmp_1 / tmp_D2S))
        tmp_y = np.sqrt(np.sum(tmp_2 / tmp_D2S))
        res = (1 - res / (tmp_x * tmp_y)) * 0.5
        return res
    
    def get_name(self) -> str:
        return "D2S"