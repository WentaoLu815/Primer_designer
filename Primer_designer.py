import csv
import numpy as np

# 读取文件并提取数据
with open('your_file_path') as f:
    reader = csv.reader(f)
    next(reader)  # 跳过表头
    Position, Identity = zip(*[(row[0], float(row[1])) for row in reader])  # 提取Position和Identity列并转换为float

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
group_of_amplicon_score = []  # 新增列表以存储扩增子的得分

# 遍历扩增子长度范围
for bp in range(bp1, bp2 + 1):
    bottom_of_forward = Position[-1] - bp + 1
    top_of_reverse = bp - length_r + 1
    bottom_of_reverse = Position[-1] - length_r + 1

    # 计算上游引物得分
    forward_scores = np.array([Identity[i:i + length_f].sum() for i in range(0, bottom_of_forward)])

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
            group_of_amplicon_score.append(amplicon_scores[index])  # 记录当前扩增子得分

# 输出到CSV文件
with open('result.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['扩增子长度', '起始位点', '最大得分', '扩增子得分'])  # 写入表头
    for n in range(len(group_of_position)):
        writer.writerow([group_of_bp[n], group_of_position[n], group_of_max[n], group_of_amplicon_score[n]])  # 写入每一行数据

print("结果已成功输出到 'result.csv' 文件中。")
