import csv
import numpy as np
import pandas as pd

# 读取文件并提取数据
with open('F:\Desktop\IBR\BoHV_Genome\BoHV_1_Identity.csv') as f:
    reader = csv.reader(f)
    next(reader)  # 跳过表头
    Position, Identity, Nucleotide = zip(*[(row[0], float(row[1]), row[2]) for row in reader])  # 提取Position、Identity和核苷酸信息

# 将Identity转换为NumPy数组，便于后续向量化操作
Identity = np.array(Identity)
Position = list(map(int, Position))  # 将Position转换为int

# 用户输入
bp1 = int(input('扩增子最小值：'))
bp2 = int(input('扩增子最大值：'))
length_f = int(input('上游引物的长度：'))
length_r = int(input('下游引物的长度：'))

# 初始化变量
v = 0
group_of_position = []
group_of_bp = []
group_of_max = []
group_of_amplicon_score = []
group_of_amplicon_sequence = []
group_of_forward_primer = []  # 新增列表用于存储上游引物序列
group_of_reverse_primer = []  # 新增列表用于存储下游引物序列

# 遍历扩增子长度范围
for bp in range(bp1, bp2 + 1):
    bottom_of_forward = Position[-1] - bp + 1
    top_of_reverse = bp - length_r + 1
    bottom_of_reverse = Position[-1] - length_r + 1

    # 计算上游引物得分
    forward_scores = []
    for i in range(0, bottom_of_forward):
        if Nucleotide[i] == 'T':  # 如果上游引物的第一个碱基为 T，得分为 0
            forward_scores.append(0)
        else:
            forward_scores.append(Identity[i:i + length_f].sum())
    forward_scores = np.array(forward_scores)

    # 计算下游引物得分
    reverse_scores = np.array([Identity[k:k + length_r].sum() for k in range(top_of_reverse - 1, bottom_of_reverse)])

    # 计算扩增子得分（从每个可能的起始位置，长度为 bp）
    amplicon_scores = np.array([Identity[j:j + bp].sum() for j in range(0, Position[-1] - bp + 1)])

    # 计算上游和下游得分之和
    f_and_r = forward_scores + reverse_scores[:len(forward_scores)]
    BPmax = f_and_r.max()

    # 记录具有最大得分的扩增子信息
    if v <= BPmax:
        v = BPmax
        max_indices = np.where(f_and_r == BPmax)[0]
        for index in max_indices:
            group_of_position.append(index + 1)
            group_of_bp.append(bp)
            group_of_max.append(BPmax)
            group_of_amplicon_score.append(amplicon_scores[index])
            group_of_amplicon_sequence.append("".join(Nucleotide[index:index + bp]))
            group_of_forward_primer.append("".join(Nucleotide[index:index + length_f]))  # 上游引物序列
            group_of_reverse_primer.append("".join(Nucleotide[index + bp - length_r:index + bp]))  # 下游引物序列

# 创建 DataFrame，并添加行号
df = pd.DataFrame({
    '行号': range(1, len(group_of_bp) + 1),  # 新增列：行号，从1开始递增
    '扩增子长度': group_of_bp,
    '起始位点': group_of_position,
    '最大得分': group_of_max,
    '扩增子得分': group_of_amplicon_score,
    '扩增子序列': group_of_amplicon_sequence,
    '上游引物': group_of_forward_primer,  # 新增列：上游引物
    '下游引物': group_of_reverse_primer  # 新增列：下游引物
})

# 创建Excel文件，并根据扩增子长度分类写入到不同的sheet
with pd.ExcelWriter('BoHV_1_result.xlsx') as writer:
    for bp_length, group in df.groupby('扩增子长度'):
        group.to_excel(writer, sheet_name=f'扩增子长度_{bp_length}', index=False)

print("结果已成功输出到 'BoHV_1_result.xlsx' 文件中，每个扩增子长度的数据在单独的sheet中。")
import csv
import numpy as np
import pandas as pd

# 读取文件并提取数据
with open('F:\Desktop\IBR\BoHV_Genome\BoHV_1_Identity.csv') as f:
    reader = csv.reader(f)
    next(reader)  # 跳过表头
    Position, Identity, Nucleotide = zip(*[(row[0], float(row[1]), row[2]) for row in reader])  # 提取Position、Identity和核苷酸信息

# 将Identity转换为NumPy数组，便于后续向量化操作
Identity = np.array(Identity)
Position = list(map(int, Position))  # 将Position转换为int  

# 用户输入
bp1 = int(input('扩增子最小值：'))
bp2 = int(input('扩增子最大值：'))
length_f = int(input('上游引物的长度：'))
length_r = int(input('下游引物的长度：'))

# 初始化变量
v = 0
group_of_position = []
group_of_bp = []
group_of_max = []
group_of_amplicon_score = []
group_of_amplicon_sequence = []
group_of_forward_primer = []  # 新增列表用于存储上游引物序列
group_of_reverse_primer = []  # 新增列表用于存储下游引物序列

# 遍历扩增子长度范围
for bp in range(bp1, bp2 + 1):
    bottom_of_forward = Position[-1] - bp + 1
    top_of_reverse = bp - length_r + 1
    bottom_of_reverse = Position[-1] - length_r + 1

    # 计算上游引物得分
    forward_scores = []
    for i in range(0, bottom_of_forward):
        if Nucleotide[i] == 'T':  # 如果上游引物的第一个碱基为 T，得分为 0
            forward_scores.append(0)
        else:
            forward_scores.append(Identity[i:i + length_f].sum())  
    forward_scores = np.array(forward_scores)

    # 计算下游引物得分
    reverse_scores = np.array([Identity[k:k + length_r].sum() for k in range(top_of_reverse - 1, bottom_of_reverse)])

    # 计算扩增子得分（从每个可能的起始位置，长度为 bp）
    amplicon_scores = np.array([Identity[j:j + bp].sum() for j in range(0, Position[-1] - bp + 1)])

    # 计算上游和下游得分之和
    f_and_r = forward_scores + reverse_scores[:len(forward_scores)]
    BPmax = f_and_r.max()

    # 记录具有最大得分的扩增子信息
    if v <= BPmax:
        v = BPmax
        max_indices = np.where(f_and_r == BPmax)[0]
        for index in max_indices:
            group_of_position.append(index + 1)
            group_of_bp.append(bp)
            group_of_max.append(BPmax)
            group_of_amplicon_score.append(amplicon_scores[index])
            group_of_amplicon_sequence.append("".join(Nucleotide[index:index + bp]))
            group_of_forward_primer.append("".join(Nucleotide[index:index + length_f]))  # 上游引物序列
            group_of_reverse_primer.append("".join(Nucleotide[index + bp - length_r:index + bp]))  # 下游引物序列

# 创建 DataFrame，并添加行号
df = pd.DataFrame({
    '行号': range(1, len(group_of_bp) + 1),  # 新增列：行号，从1开始递增
    '扩增子长度': group_of_bp,
    '起始位点': group_of_position,
    '最大得分': group_of_max,
    '扩增子得分': group_of_amplicon_score,
    '扩增子序列': group_of_amplicon_sequence,
    '上游引物': group_of_forward_primer,  # 新增列：上游引物
    '下游引物': group_of_reverse_primer  # 新增列：下游引物
})

# 创建Excel文件，并根据扩增子长度分类写入到不同的sheet
with pd.ExcelWriter('BoHV_1_result.xlsx') as writer:
    for bp_length, group in df.groupby('扩增子长度'):
        group.to_excel(writer, sheet_name=f'扩增子长度_{bp_length}', index=False)

print("结果已成功输出到 'BoHV_1_result.xlsx' 文件中，每个扩增子长度的数据在单独的sheet中。")
