import csv


# 读取CSV文件的第二列数据
def read_second_column(file_path):
    with open(file_path, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # 跳过表头
        return [int(row[1]) for row in reader]


# 查找并记录连续数值的起始和终止值
def find_continuous_ranges(numbers):
    if not numbers:
        return []

    ranges = []
    start = numbers[0]
    end = numbers[0]

    for n in numbers[1:]:
        if n == end + 1:  # 检查数字是否是连续的
            end = n
        else:
            # 只有在 start 和 end 不同的时候才记录区间
            if start != end:
                ranges.append((start, end))
            else:
                ranges.append((start, start))  # 单独的数字区间
            start = n  # 更新新的连续区间的起始值
            end = n  # 重置end为新的起始值

    # 添加最后一个连续区间
    if start != end:
        ranges.append((start, end))
    else:
        ranges.append((start, start))

    return ranges


# 写入合并结果到新的CSV文件
def write_ranges_to_csv(ranges, output_file_path):
    with open(output_file_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Start', 'End'])  # 写入表头
        for start, end in ranges:
            writer.writerow([start, end])


# 主程序
input_file_path = 'input_file_path'  # 替换为实际文件路径
output_file_path = 'output_file_path'  # 输出文件路径

# 读取第二列数据
numbers = read_second_column(input_file_path)
# 对数据进行排序，确保连续数值的顺序
numbers.sort()
# 查找连续区间
ranges = find_continuous_ranges(numbers)
# 将结果写入新的CSV文件
write_ranges_to_csv(ranges, output_file_path)

print("连续区间已写入到 'output_ranges.csv' 文件中。")
