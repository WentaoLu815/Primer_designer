import pandas as pd

# 定义互补配对规则
complementary_base = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}


def is_complementary(seq1, seq2):
    """检查 seq1 和 seq2 是否为互补序列"""
    return all(complementary_base[base1] == base2 for base1, base2 in zip(seq1, seq2))


def find_complementary_positions(primer_sequence):
    """在给定的下游引物序列中查找大于四个碱基的连续互补配对子序列的位置"""
    length = len(primer_sequence)
    max_match_length = 0  # 记录最长互补碱基个数
    position_pair = None  # 记录符合条件的位点对

    for match_length in range(4, length + 1):  # 从4个碱基开始递增检查更长的互补序列
        for i in range(length - match_length + 1):  # 从每个位置检查长度为 match_length 的子序列
            forward_seq = primer_sequence[i:i + match_length]
            for j in range(i + match_length, length - match_length + 1):
                reverse_seq = primer_sequence[j:j + match_length]

                # 检查正向互补或反向互补
                if is_complementary(forward_seq, reverse_seq) or is_complementary(forward_seq, reverse_seq[::-1]):
                    max_match_length = match_length
                    position_pair = f"{i + 1},{j + 1}"
                    break  # 找到更长的匹配，跳出内层循环

            if max_match_length == match_length:
                break  # 找到更长的匹配，跳出外层循环

    return position_pair, max_match_length if position_pair else (None, 0)


# 读取上游引物处理后的结果文件
df = pd.read_excel('BoHV_1_result_upstream_complementary.xlsx', sheet_name=None)

# 处理下游引物数据
for sheet_name, sheet_df in df.items():
    downstream_results = sheet_df['下游引物'].apply(find_complementary_positions)
    sheet_df['下游引物互补配对位点'] = downstream_results.apply(lambda x: x[0])
    sheet_df['下游引物互补碱基个数'] = downstream_results.apply(lambda x: x[1])

# 保存最终结果到新的 Excel 文件
with pd.ExcelWriter('BoHV_1_result_final_complementary.xlsx') as writer:
    for sheet_name, sheet_df in df.items():
        sheet_df.to_excel(writer, sheet_name=sheet_name, index=False)

print("下游引物处理完成，最终结果已保存至 'BoHV_1_result_final_complementary.xlsx' 文件中。")
