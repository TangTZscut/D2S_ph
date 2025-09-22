from typing import List, Tuple
from unittest.mock import Base
import numpy as np

from .BaseSimilarity import BaseSimilarity
from Data import Sequence

class D2Similarity(BaseSimilarity):
    
    def calculate(self, seq_1: Sequence, seq_2: Sequence, k: int, r: int = -1) -> float:
        x_kmer = seq_1.get_kmer(k)
        y_kmer = seq_2.get_kmer(k)
        x = x_kmer.data
        y = y_kmer.data
        tmp_x = np.power(x, 2)
        tmp_x = np.sum(tmp_x)
        tmp_y = np.power(y, 2)
        tmp_y = np.sum(tmp_y)
        tmp = np.sum(x * y)
        res = tmp / (np.sqrt(tmp_x) * np.sqrt(tmp_y))
        res = 0.5 * (1 - res)
        return res
    
    def get_name(self) -> str:
        return "D2"