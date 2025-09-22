import os
from unittest import result
import pandas as pd
import json
import glob
import traceback
from multiprocessing import Pool, cpu_count, pool
from Executor import Executor
from Similarity import *
from typing import Dict, List, Tuple
from Data import FullSequence
import time

# 全局参数
k_list = [5]
r_list = [1]
windows = 400
shift = 400
function_list = [D2StarSimilarity()]
save_dir = "./run/Pn_1k"
# save_dir=r"D:\Tmp\Biology\py"
seq_pairs_file = "./run/data_pn/seq_pairs.json"
THREAD_COUNT = 8  # 进程数，根据机器状态修改

# 创建结果目录
os.makedirs(save_dir, exist_ok=True)

def analyze_pair(dic: Dict[str, str]) -> str:
    seq_1 = dic['seq_1']
    seq_2 = dic['seq_2']
    try:
        start = time.perf_counter()
        executor = Executor(k_list, r_list, windows, shift, function_list, seq_1, seq_2, save_dir)
        executor.execute2()
        end = time.perf_counter()

        res = f"OK:{seq_1} vs {seq_2}"
        print(f"执行时间: {end - start:.6f} 秒")
        print(res)
    except Exception as e:
        res = f"Error:{seq_1} vs {seq_2}"
        print(res)
        traceback.print_exc()


def show_batch_results(save_dir):
    print("\n=== 批量分析结果统计 ===")

    result_files = glob.glob(f"{save_dir}/*.json")
    if result_files:
        print(f"生成了 {len(result_files)} 个结果文件:")
        for i, file in enumerate(result_files):
            try:
                with open(file, 'r') as f:
                    result = json.load(f)
                    print(f"  {i+1}. {result['seqX']} → {result['seqY']}: {len(result['features'])} 个特征值")
            except Exception as e:
                print(f"  {i+1}. {file}: 读取失败 - {e}")

        # 示例展示第一个结果
        try:
            with open(result_files[0], 'r') as f:
                result = json.load(f)
                print(f"\n示例结果文件内容 ({os.path.basename(result_files[0])}):")
                print(f"  序列X: {result['seqX']}")
                print(f"  序列Y: {result['seqY']}")
                print(f"  k值: {result['k']}")
                print(f"  r值: {result['r']}")
                print(f"  算法: {result['function']}")
                print(f"  窗口大小: {result['windows']}")
                print(f"  步长: {result['shift']}")
                print(f"  特征值数量: {len(result['features'])}")
                print(f"  前5个特征值: {result['features'][:5]}")
        except Exception as e:
            print(f"  读取文件出错: {e}")
    else:
        print("未找到结果文件。")


def prepare_sequences(seq: FullSequence):
    """
    计算每个序列切分的子序列的kmer和markov，保存npy文件
    """
    current_part_name = ""
    try:
        seq.load_parts(windows, shift)
        for part in seq.sequence_parts:
            current_part_name = part.file_name
            for k in k_list:
                part.get_kmer(k)
                for r in r_list:
                    part.get_markov(k, r)
                    part.remove_markov(k, r)
                part.remove_kmer(k)
        seq.clear()
        print(f"OK:{seq.file_name}")
    except Exception as e:
        print(f"Error:{seq.file_name}, current_part_name: {current_part_name}")
        traceback.print_exc()


def get_full_sequence_list() -> List[FullSequence]:
    with open(seq_pairs_file, "r", encoding="utf-8") as f:
        seq_pairs: List[Dict[str, str]] = json.load(f) 
        dic = dict()
        for pair in seq_pairs:
            seq_1 = pair['seq_1']
            seq_2 = pair['seq_2']
            dic[seq_1] = 1
            dic[seq_2] = 1             
        
        full_seq_list = list()
        for seq in dic.keys():
            full_seq = FullSequence(seq)
            full_seq_list.append(full_seq)
        return full_seq_list


def prepare() -> None:
    full_seq_list = get_full_sequence_list()  # type: List[FullSequence]
    with Pool(processes=THREAD_COUNT) as pool:
        pool.map(prepare_sequences, full_seq_list)
    print("所有序列对准备完成。")


def divide_seq(seq: FullSequence) -> None:
    """
    切分序列
    """
    seq.divide(windows, shift)


def divide() -> None: 
    full_seq_list = get_full_sequence_list()  # type: List[FullSequence]
    with Pool(processes=THREAD_COUNT) as pool:
        pool.map(divide_seq, full_seq_list)
    print("切分完成")
         
                 

def batch_analysis():
    with open(seq_pairs_file, "r", encoding="utf-8") as f:
        seq_pairs: List[Dict[str, str]] = json.load(f)
        
        current = len(os.listdir(save_dir)) / 2
        print(f"准备处理 {len(seq_pairs)} 个序列对。已处理：{current}")

 
        with Pool(processes=THREAD_COUNT) as pool:
            pool.map(analyze_pair, seq_pairs)

        print("所有序列对分析完成。")



if __name__ == "__main__":

    # 首先执行 get_data_pn.py，这个是为了从里面拆出每个基因序列

    # 然后下面的分步执行，执行完一个再进行下一个，中途可以中断，里面有检测，处理过的数据不会二次处理

    # 第一步切分序列
    divide()

    # 第二步，计算kmer和markov
    prepare()

    # 第三步，计算结果
    batch_analysis()


    