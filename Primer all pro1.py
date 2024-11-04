import pandas as pd

# 定义碱基互补配对规则
complementary_base = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}


def convert_to_complementary(sequence):
    """根据碱基互补规则将序列转换为互补序列"""
    return ''.join(complementary_base[base] for base in sequence)


def find_complementary_matches(upstream_sequence, downstream_complementary_sequence):
    """检测上游引物和下游引物互补序列是否存在大于等于4个碱基的连续互补配对"""
    max_match_length = 0  # 记录最长互补配对的碱基个数
    position_pair = None  # 记录符合条件的位点对

    # 迭代上游引物和下游引物转换后的互补序列，找出最大互补配对
    for i in range(len(upstream_sequence) - 3):  # 从每个位置检查4个及以上碱基
        for j in range(len(downstream_complementary_sequence) - 3):
            match_length = 0  # 当前连续互补配对的碱基个数

            # 逐个比较碱基，检查互补配对
            while (i + match_length < len(upstream_sequence) and
                   j + match_length < len(downstream_complementary_sequence) and
                   upstream_sequence[i + match_length] == downstream_complementary_sequence[j + match_length]):
                match_length += 1

            # 如果找到大于或等于4个的连续互补配对
            if match_length >= 4 and match_length > max_match_length:
                max_match_length = match_length
                position_pair = (i + 1, j + 1)  # 记录1基准的起始位点

    return position_pair, max_match_length


# 读取 Excel 文件
df = pd.read_excel('BoHV_1_result_final_complementary.xlsx', sheet_name=None)

# 遍历每个sheet并处理上游和下游引物互补配对
for sheet_name, sheet_df in df.items():
    results = []
    for idx, row in sheet_df.iterrows():
        # 获取上游和下游引物序列
        upstream_sequence = row['上游引物']
        downstream_sequence = row['下游引物']

        # 将下游引物转换为互补序列
        downstream_complementary_sequence = convert_to_complementary(downstream_sequence[::-1])

        # 检查上游引物和下游引物的互补配对
        position_pair, match_length = find_complementary_matches(upstream_sequence, downstream_complementary_sequence)

        # 将结果添加到列表中
        results.append({
            '上游引物互补位点': position_pair[0] if position_pair else None,
            '下游引物互补位点': position_pair[1] if position_pair else None,
            '互补碱基个数': match_length
        })

    # 将结果添加到 DataFrame 中
    result_df = pd.DataFrame(results)
    sheet_df['上游引物互补位点'] = result_df['上游引物互补位点']
    sheet_df['下游引物互补位点'] = result_df['下游引物互补位点']
    sheet_df['互补碱基个数'] = result_df['互补碱基个数']

# 将最终结果写入一个新的 Excel 文件
with pd.ExcelWriter('BoHV_1_final_complementary_match.xlsx') as writer:
    for sheet_name, sheet_df in df.items():
        sheet_df.to_excel(writer, sheet_name=sheet_name, index=False)

print("处理完成，结果已保存至 'BoHV_1_final_complementary_match.xlsx' 文件中。")
