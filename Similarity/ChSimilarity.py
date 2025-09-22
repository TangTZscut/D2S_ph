from typing import List, Tuple
from unittest.mock import Base
import numpy as np

from .BaseSimilarity import BaseSimilarity
from Data import Sequence

class ChSimilarity(BaseSimilarity):
    
    def calculate(self, seq_1: Sequence, seq_2: Sequence, k: int, r: int = -1) -> float:
        x_kmer = seq_1.get_kmer(k)
        y_kmer = seq_2.get_kmer(k)
        x = x_kmer.data / x_kmer.total
        y = y_kmer.data / y_kmer.total
        z = np.absolute(x - y)
        res = np.max(z)
        return res
    
    def get_name(self) -> str:
        return "Ch"