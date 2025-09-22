from typing import List, Tuple
from unittest.mock import Base
import numpy as np

from .BaseSimilarity import BaseSimilarity
from Data import Sequence

class D2StarSimilarity(BaseSimilarity):
    
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
        
        tmp_3 = np.sqrt(tmp_x * tmp_y)
        tmp_x = np.where(tmp_x == 0, 1, tmp_x)
        tmp_y = np.where(tmp_y == 0, 1, tmp_y)
        tmp_3 = np.where(tmp_3 == 0, 1, tmp_3)
        
        res = np.sum(x_bar * y_bar / tmp_3)
        tmp_x = np.sqrt(np.sum(tmp_1 / tmp_x))
        tmp_y = np.sqrt(np.sum(tmp_2 / tmp_y))
        res = 0.5 * (1 - res / ( tmp_x * tmp_y))
        return res
    
    def get_name(self) -> str:
        return "D2Star"