from Data import FullSequence
from Similarity import BaseSimilarity
import numpy as np
from typing import List
import os
import json
from Utils.PathUtils import PathUtils


class Executor:
    """
    对两条序列按多个 (k, r, 相似度函数) 组合进行滑窗比较，
    并把任意大小的窗口相似度矩阵压缩为**固定 11 维特征向量**：
      [max, mean, std, p25, p50, p75, top1-top5]
    """

    def __init__(
        self,
        k_list: List[int],
        r_list: List[int],
        winodws: int,
        shift: int,
        function_list: List[BaseSimilarity],
        seq_1: str,
        seq_2: str,
        save_dir: str,
    ):
        """
        Args:
            k_list:      要遍历的 k-mer 列表
            r_list:      D2S / D2Star 的马尔可夫阶 r 列表
            winodws:     滑动窗口长度
            shift:       窗口步长
            function_list: 相似度函数对象列表
            seq_1:       第一条序列文件路径
            seq_2:       第二条序列文件路径
            save_dir:    结果输出目录
        """
        self._k_list = k_list
        self._r_list = r_list
        self._windows = winodws
        self._shift = shift
        self._function_list = function_list
        self._seq_1 = FullSequence(seq_1)
        self._seq_2 = FullSequence(seq_2)
        self._save_dir = save_dir

    # --------------------------------------------------------------------- #
    # 内部工具函数
    # --------------------------------------------------------------------- #
    @staticmethod
    def _aggregate(matrix: np.ndarray, top_k: int = 5) -> np.ndarray:
        """
        把任意尺寸的相似度矩阵压缩成固定长度向量:
        [max, mean, std, p25, p50, p75, top1 … top_k]

        Returns
        -------
        np.ndarray
            长度 6 + top_k (默认 11) 的 1-D 特征向量
        """
        flat = matrix.ravel()
        if flat.size == 0:
            return np.zeros(6 + top_k, dtype=np.float64)

        stats = np.array(
            [
                flat.max(),
                flat.mean(),
                flat.std(ddof=0),
                np.percentile(flat, 25),
                np.percentile(flat, 50),
                np.percentile(flat, 75),
            ],
            dtype=np.float64,
        )

        top_vals = np.sort(flat)[-top_k:][::-1]  # 从大到小
        if top_vals.size < top_k:  # 短序列补零
            top_vals = np.pad(top_vals, (0, top_k - top_vals.size))

        return np.concatenate([stats, top_vals])

    def _write_result(self, result: dict) -> None:
        """保存 JSON 结果文件"""
        if result["r"] != -1:
            name = (
                f"{result['seqX']}_to_{result['seqY']}_k{result['k']}"
                f"_r{result['r']}_{result['function']}.json"
            )
        else:
            name = (
                f"{result['seqX']}_to_{result['seqY']}_k{result['k']}"
                f"_{result['function']}.json"
            )
        os.makedirs(self._save_dir, exist_ok=True)
        with open(os.path.join(self._save_dir, name), "w") as f:
            json.dump(result, f, indent=2)

    # --------------------------------------------------------------------- #
    # 相似度计算主流程
    # --------------------------------------------------------------------- #
    def _compare_sequences(self, k: int, function: BaseSimilarity, r: int = -1):
        """计算两条序列所有窗口对的相似度，并聚合成固定维特征"""
        rows = len(self._seq_1.sequence_parts)
        cols = len(self._seq_2.sequence_parts)
        sim_mat = np.zeros((rows, cols), dtype=np.float64)

        # 逐窗口计算
        for i, part_1 in enumerate(self._seq_1.sequence_parts):
            for j, part_2 in enumerate(self._seq_2.sequence_parts):
                sim_mat[i, j] = function.calculate(part_1, part_2, k, r)
                part_2.remove_kmer(k)
                part_2.remove_markov(k, r)
            part_1.remove_kmer(k)
            part_1.remove_markov(k, r)

        # ① seq1→seq2   ② seq2→seq1 (矩阵转置)
        features_1to2 = self._aggregate(sim_mat)
        # features_2to1 = self._aggregate(sim_mat.T)

        meta = {
            "k": k,
            "r": r,
            "windows": self._windows,
            "shift": self._shift,
            "function": function.get_name(),
        }

        self._write_result(
            {
                **meta,
                "seqX": self._seq_1.file_name,
                "seqY": self._seq_2.file_name,
                "features": features_1to2.tolist(),
            }
        )
        # self._write_result(
        #     {
        #         **meta,
        #         "seqX": self._seq_2.file_name,
        #         "seqY": self._seq_1.file_name,
        #         "features": features_2to1.tolist(),
        #     }
        # )

    # --------------------------------------------------------------------- #
    # 外部调用入口
    # --------------------------------------------------------------------- #
    def execute(self) -> None:
        """切分窗口 → 遍历所有 (k,r,函数) → 输出固定维特征 JSON"""
        self._seq_1.divide(self._windows, self._shift)
        self._seq_2.divide(self._windows, self._shift)

        for k in self._k_list:
            for func in self._function_list:
                if func.get_name() in {"D2S", "D2Star"}:
                    for r in self._r_list:
                        self._compare_sequences(k, func, r)
                else:
                    self._compare_sequences(k, func)

        PathUtils.del_tmp_dir()
