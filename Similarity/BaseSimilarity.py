from Data import Sequence
from abc import ABC, abstractmethod
from typing import List, Tuple

class BaseSimilarity(ABC):
    
    @abstractmethod
    def calculate(self, seq_1: Sequence, seq_2: Sequence, k: int, r: int = -1) -> float:
        pass
    
    def calculate_one_to_N(self, seq_1: Sequence, seq_list: List[Sequence], k: int, r: int = -1) -> List[float]:
        return [self.calculate(seq_1, seq_2, k, r) for seq_2 in seq_list]
    
    @abstractmethod
    def get_name(self) -> str:
        pass