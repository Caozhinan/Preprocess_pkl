import re

with open('test_time.txt', encoding='utf-8') as f:
    log = f.read()

def sum_times(pattern):
    times = re.findall(pattern, log)
    return sum(float(x) for x in times)

src_sum = sum_times(r'find_closest_atom_by_coord\(src\): ([0-9.]+)s')
tgt_sum = sum_times(r'find_closest_atom_by_coord\(tgt\): ([0-9.]+)s')
analyze_sum = sum_times(r'analyze_atom_pair: ([0-9.]+)s')
single_sum = sum_times(r'single edge total: ([0-9.]+)s')

print(f'find_closest_atom_by_coord(src) 总时间: {src_sum:.6f} s')
print(f'find_closest_atom_by_coord(tgt) 总时间: {tgt_sum:.6f} s')
print(f'analyze_atom_pair 总时间: {analyze_sum:.6f} s')
print(f'single edge total 总时间: {single_sum:.6f} s')
print(f'所有src+tgt+analyze+single edge合计: {src_sum + tgt_sum + analyze_sum + single_sum:.6f} s')