#!/usr/bin/env python3
"""
序列相似性分析示例脚本
演示如何使用Executor类进行DNA序列相似性分析
"""

from Executor import Executor
from Similarity import *
import os

def custom_example():
    """自定义参数示例"""
    print("\n=== 自定义参数分析示例 ===")
    
    # 1. 自定义参数设置
    k_list = [6]             # 只分析6-mer
    r_list = [2]             # 只使用2阶马尔可夫模型
    windows = 1000            # 大窗口500bp
    shift = 500              # 无重叠切分
    
    # 2. 只使用特定算法
    function_list = [
        # D2Similarity(),      # D2统计量
        # D2SSimilarity(),     # D2S统计量
        D2StarSimilarity()   # D2Star统计量
    ]
    
    # 3. 设置输入输出路径
    seq_1 = "seq1-1.fasta"
    seq_2 = "seq1-2.fasta"
    save_dir = "results_Phage_Host7_2"
    
    # 创建结果目录
    os.makedirs(save_dir, exist_ok=True)
    
    # 4. 执行分析
    print(f"分析序列: {seq_1} vs {seq_2}")
    print(f"窗口大小: {windows}bp, 步长: {shift}bp")
    print(f"k-mer长度: {k_list}")
    print(f"马尔可夫阶数: {r_list}")
    print(f"算法: {[f.get_name() for f in function_list]}")
    
    executor = Executor(k_list, r_list, windows, shift, function_list, seq_1, seq_2, save_dir)
    executor.execute()
    
    print(f"分析完成！结果保存在: {save_dir}/")

def show_results():
    """展示结果文件"""
    print("\n=== 结果文件示例 ===")
    
    import json
    import glob
    
    # 查找结果文件
    result_files = glob.glob("results_Phage_Host2/*.json")
    
    if result_files:
        print(f"找到 {len(result_files)} 个结果文件:")
        for file in result_files[:3]:  # 只显示前3个
            print(f"  - {file}")
        
        # 显示第一个结果文件的内容
        if result_files:
            print(f"\n示例结果文件内容 ({result_files[0]}):")
            try:
                with open(result_files[0], 'r') as f:
                    result = json.load(f)
                    print(f"  序列X: {result['seqX']}")
                    print(f"  序列Y: {result['seqY']}")
                    print(f"  k值: {result['k']}")
                    print(f"  r值: {result['r']}")
                    print(f"  算法: {result['function']}")
                    print(f"  窗口大小: {result['windows']}")
                    print(f"  步长: {result['shift']}")
                    print(f"  相似性得分数量: {len(result['scores'])}")
                    print(f"  前5个得分: {result['scores'][:5]}")
            except Exception as e:
                print(f"  读取文件出错: {e}")
    else:
        print("未找到结果文件，请先运行分析。")

if __name__ == "__main__":
    # # 检查测试文件是否存在
    # if not os.path.exists("test_seq1.txt") or not os.path.exists("test_seq2.txt"):
    #     print("错误: 测试序列文件不存在！")
    #     print("请确保 test_seq1.txt 和 test_seq2.txt 文件存在。")
    #     exit(1)
    
    try:
        custom_example()
        # show_results()
        
        print("\n=== 所有示例运行完成 ===")
        print("你可以查看各个results_*目录中的JSON文件来分析结果。")
        
    except Exception as e:
        print(f"运行出错: {e}")
        print("请检查依赖模块是否正确安装。")