import csv
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

# 读取文件并提取数据
with open('yourfilepath') as f:    
    reader = csv.reader(f)
    next(reader)  # 跳过表头
    Position, Identity, Nucleotide = zip(
        *[(row[0], float(row[1]), row[2]) for row in reader])  # 提取Position、Identity和核苷酸信息

# 将Identity转换为NumPy数组，便于后续向量化操作
Identity = np.array(Identity)
Position = list(map(int, Position))  # 将Position转换为int  

# 用户输入
bp1 = int(input('扩增子最小值：'))
bp2 = int(input('扩增子最大值：'))
length_f = int(input('上游引物的长度：'))
length_r = int(input('下游引物的长度：'))
modification = input("上游引物添加修饰（输入碱基序列或输入 'NA' 跳过修饰）：")

# 定义互补碱基映射
complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# 检查序列中的GC含量是否在30%-70%之间
def calculate_gc_content(sequence):
    gc_count = sum(1 for base in sequence if base in "GC")
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0

# 检查是否有连续4个及以上互补碱基
def check_self_complementarity(sequence):
    count = 0
    for i in range(len(sequence) - 1):
        if sequence[i] == complement_map.get(sequence[i + 1], None):
            count += 1
            if count >= 3:
                return 0  # 存在连续4个互补碱基
        else:
            count = 0
    return 1  # 不存在连续4个互补碱基

# 检查上游引物和下游引物转换物之间是否有连续3个以上的碱基互补配对
def check_interaction_score(forward_seq, reverse_complement_seq):
    count = 0
    min_length = min(len(forward_seq), len(reverse_complement_seq))
    for i in range(min_length):
        if forward_seq[i] == complement_map.get(reverse_complement_seq[i], None):
            count += 1
            if count >= 3:
                return 0  # 存在连续3个以上的碱基互补配对
        else:
            count = 0
    return 1  # 不存在连续3个以上的碱基互补配对

# 将下游引物转换为互补序列
def get_complement_sequence(sequence):
    return ''.join(complement_map.get(base, base) for base in sequence)

# 处理单个扩增子长度的计算逻辑
def process_bp(bp):
    bottom_of_forward = Position[-1] - bp + 1
    top_of_reverse = bp - length_r + 1
    bottom_of_reverse = Position[-1] - length_r + 1

    # 初始化变量 v
    v = 0

    result_data = {
        '扩增子长度': [],
        '起始位点范围': [],
        '最大得分': [],
        '扩增子得分': [],
        '扩增子序列': [],
        '上游引物序列': [],
        '下游引物转换物序列': [],
        '上游引物自身配对得分': [],
        '下游引物转换配对得分': [],
        '上游引物与下游引物转换物互作得分': [],
        '修饰后上游引物自身互补配对得分': [],
        '修饰后上游引物与下游引物互作得分': [],
        '上游引物GC含量': [],  # 新增列
        '下游引物GC含量': []  # 新增列
    }

    # 计算上游引物得分
    forward_scores = []
    forward_self_scores = []
    forward_sequences = []
    modified_forward_self_scores = []
    modified_interaction_scores = []
    for i in range(0, bottom_of_forward):
        forward_seq = Nucleotide[i:i + length_f]
        forward_seq_str = "".join(forward_seq)  # 将上游引物序列转换为字符串
        forward_sequences.append(forward_seq_str)  # 记录上游引物序列  

        if Nucleotide[i + length_f - 1] == 'A':
            score = 0
        else:
            score = Identity[i:i + length_f].sum()
            if 30 <= calculate_gc_content(forward_seq_str) <= 70:
                score += 1
        forward_scores.append(score)
        forward_self_scores.append(check_self_complementarity(forward_seq_str))

        # 修饰上游引物序列
        if modification != "NA":
            modified_forward_seq = modification + forward_seq_str
            modified_forward_self_scores.append(check_self_complementarity(modified_forward_seq))  # 修饰后上游引物自身互补配对得分
        else:
            modified_forward_self_scores.append("NA")  # 如果不修饰，则得分为 NA

    forward_scores = np.array(forward_scores)

    # 计算下游引物得分
    reverse_scores = []
    reverse_conversion_scores = []
    reverse_converted_sequences = []
    interaction_scores = []
    for k in range(top_of_reverse - 1, bottom_of_reverse):
        reverse_seq = Nucleotide[k:k + length_r]
        reverse_seq_str = "".join(reverse_seq)  # 转换为字符串
        if Nucleotide[k] == 'T':
            score = 0
        else:
            score = Identity[k:k + length_r].sum()
            if 30 <= calculate_gc_content(reverse_seq_str) <= 70:
                score += 1
        reverse_scores.append(score)

        # 计算下游引物转换配对得分
        reverse_complement_seq = get_complement_sequence(reverse_seq_str)
        reverse_converted_sequences.append(reverse_complement_seq)  # 记录下游引物转换物序列
        reverse_conversion_scores.append(check_self_complementarity(reverse_complement_seq))

        # 计算上游引物与下游引物转换物的互作得分
        interaction_scores.append(check_interaction_score(forward_seq_str, reverse_complement_seq))

        # 计算修饰后上游引物与下游引物的互作得分
        if modification != "NA":
            modified_forward_seq = modification + forward_seq_str
            modified_interaction_scores.append(check_interaction_score(modified_forward_seq, reverse_complement_seq))
        else:
            modified_interaction_scores.append("NA")

    reverse_scores = np.array(reverse_scores)

    # 计算扩增子得分
    amplicon_scores = np.array([Identity[j:j + bp].sum() for j in range(0, Position[-1] - bp + 1)])

    # 计算上游和下游得分之和
    f_and_r = forward_scores + reverse_scores[:len(forward_scores)]
    BPmax = f_and_r.max()

    # 记录具有最大得分的扩增子信息
    if v <= BPmax:
        v = BPmax
        max_indices = np.where(f_and_r == BPmax)[0]

        # 找出连续的起始位点并合并
        combined_start_positions = []
        for index in max_indices:
            start_position = index + 1
            if not combined_start_positions or start_position != combined_start_positions[-1][-1] + 1:
                combined_start_positions.append([start_position])
            else:
                combined_start_positions[-1].append(start_position)

        # 合并范围内的序列和得分信息
        for position_range in combined_start_positions:
            min_position = position_range[0]
            max_position = position_range[-1]
            combined_sequence = "".join(Nucleotide[min_position - 1:max_position + bp - 1])

            # 计算上游引物和下游引物的GC含量
            forward_gc_content = calculate_gc_content(forward_sequences[min_position - 1])  
            reverse_gc_content = calculate_gc_content(reverse_converted_sequences[min_position - 1])

            result_data['扩增子长度'].append(bp)
            result_data['起始位点范围'].append(f"{min_position}-{max_position}")
            result_data['最大得分'].append(BPmax)
            result_data['扩增子得分'].append(amplicon_scores[min_position - 1])
            result_data['扩增子序列'].append(combined_sequence)
            result_data['上游引物序列'].append(forward_sequences[min_position - 1])  
            result_data['下游引物转换物序列'].append(reverse_converted_sequences[min_position - 1])
            result_data['上游引物自身配对得分'].append(forward_self_scores[min_position - 1])
            result_data['下游引物转换配对得分'].append(reverse_conversion_scores[min_position - 1])
            result_data['上游引物与下游引物转换物互作得分'].append(interaction_scores[min_position - 1])
            result_data['修饰后上游引物自身互补配对得分'].append(modified_forward_self_scores[min_position - 1])
            result_data['修饰后上游引物与下游引物互作得分'].append(modified_interaction_scores[min_position - 1])
            result_data['上游引物GC含量'].append(forward_gc_content)  # 添加上游引物GC含量
            result_data['下游引物GC含量'].append(reverse_gc_content)  # 添加下游引物GC含量

    return result_data

# 创建并行执行器和进度条
results = []
with ThreadPoolExecutor() as executor:
    futures = [executor.submit(process_bp, bp) for bp in range(bp1, bp2 + 1)]
    for future in tqdm(futures, desc="计算进度", unit="任务"):
        results.append(future.result())

# 合并结果
df = pd.concat([pd.DataFrame(res) for res in results], ignore_index=True)

# 保存结果到Excel
with pd.ExcelWriter('output.xlsx') as writer:  
    for bp_length, group in df.groupby('扩增子长度'):
        group.to_excel(writer, sheet_name=f'扩增子长度_{bp_length}', index=False)

print("结果已成功输出到 'output.xlsx' 文件中，每个扩增子长度的数据在单独的sheet中。")
