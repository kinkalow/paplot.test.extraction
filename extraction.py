import pandas as pd

# set
node_size = 6503523

# utils
def p(countent, msg='', newline=True):
    print('# ' + msg)
    print(countent)
    if newline:
        print()

def get_node(bp):
    return str(bp // node_size).zfill(4)

# check
#p(get_node(6503522), '6503522')    # 0000
#p(get_node(6503523), '6503523')    # 0001
#p(get_node(13007045), '13007045')  # 0001
#p(get_node(13007046), '13007046')  # 0002

# 確認用: ノードの切れ目
sections = [node_size * i - 1 for i in range(40)]
node2sec = {get_node(i + 1): [i + 1, j] for i, j in zip(sections, sections[1:])}
p(sections, 'sections')
p(node2sec, 'node2sec')

# hg19をロード
# Node2(end node)を調整するのに必要
filepath = '../paplot/scripts/paplot/templates/genome_size_hg19.csv'
hg19_df = pd.read_csv(filepath, header=None, names=['Chr', 'Size'])[:24]  # Chr=1からYまで
hg19_df = hg19_df.astype({'Chr': 'object', 'Size': int})
hg19_df['Node'] = [get_node(bp) for bp in hg19_df['Size']]
p(hg19_df, 'hg19_df')
p(hg19_df.dtypes, 'hg19_df.dtypes')
chr_list = hg19_df['Chr'].tolist()

# ncan
filepath = '../data/ncan/data_ca.csv'
ncan_df = pd.read_csv(filepath, header=0)
df = ncan_df.copy()
df = df[['Sample', 'Chr1', 'Break1', 'Chr2', 'Break2']]
p(df.dtypes, 'df.dtypes')
# 行を削除
p(len(df), 'len(df) ... オリジナルのデータ数')
df = df[df['Chr1'].isin(chr_list)]; p(len(df), 'len(df) ... Chr1のカラム要素に対象でない要素が含まれていたらその行を削除する')
df = df[df['Chr2'].isin(chr_list)]; p(len(df), 'len(df) ... Chr2のカラム要素に対象でない要素が含まれていたらその行を削除する')
df['Node1'] = [get_node(bp) for bp in df['Break1']]
df['Node2'] = [get_node(bp) for bp in df['Break2']]
#df = df.reset_index(drop=True)
# 同じノードを抽出
same_node_df = df[(df['Node1'] == df['Node2']) & (df['Chr1'] == df['Chr2'])].copy()
p(len(same_node_df), 'len(same_node_df)')
p(same_node_df.head(), 'same_node_df.head()')
# Node2を調整
node2fix = []
for _, row in same_node_df.iterrows():
    chr = row['Chr2']
    node2 = row['Node2']
    node = hg19_df[chr == hg19_df['Chr']]['Node'].values[0]
    if node2 == node:
        node2fix.append(str(int(node2) - 1).zfill(4))
    else:
        node2fix.append(str(int(node2) + 1).zfill(4))
same_node_df['Node2Fix'] = node2fix
p(same_node_df[['Node2', 'Node2Fix']], "same_node_df[['Node2', 'Node2Fix']]")
# Node2FixをオリジナルデータのNode2に代入する
df['Node2Fix'] = df['Node2']
for idx, row in same_node_df.iterrows():
    df.at[idx, 'Node2Fix'] = row['Node2Fix']
p(len(df['Node2'] != df['Node2Fix']), "len(df['Node2'] != df['Node2Fix'])")
df['Node2'] = df['Node2Fix']
df = df.drop('Node2Fix', axis=1)

# break pointのstartとendの位置をチェックする
# 位置が間違っている場合、Break1などの並び順を変える
break1 = df['Break1'].tolist()
break2 = df['Break2'].tolist()
chr1 = df['Chr1'].tolist()
chr2 = df['Chr2'].tolist()
node1 = df['Node1'].tolist()
node2 = df['Node2'].tolist()
b1_start = [True if c1 < c2 and (c1 == c2 and b1 < b2) else False for b1, b2, c1, c2 in zip(break1, break2, chr1, chr2)]
org_df = df.copy()
df['Break1'] = [b1 if b1_start else b2 for b1, b2 in zip(break1, break2)]
df['Break2'] = [b2 if b1_start else b1 for b1, b2 in zip(break1, break2)]
df['Chr1'] = [c1 if b1_start else c2 for c1, c2 in zip(chr1, chr2)]
df['Chr2'] = [c2 if b1_start else c1 for c1, c2 in zip(chr1, chr2)]
df['Node1'] = [n1 if b1_start else n2 for n1, n2 in zip(node1, node2)]
df['Node2'] = [n2 if b1_start else n1 for n1, n2 in zip(node1, node2)]

# 最終データを出力
# 不要な行を削除してNode2Fixを追加した
# startとendを考慮し並び順を変えた
p(df, '最終データ')
p(df.dtypes, 'df.dtypes')

# =============================================================================
# =============================================================================
# =============================================================================

print('#========================== 例 ==========================\n')
pd.set_option('display.max_rows', None)

# 例1: circos1:
ch_s = 1          # start
bp_s = 185183787  # start 0028
ch_e = 1          # end
bp_e = 209060188  # end   0032
# 例2: circos1:
ch_s = 4          # start
bp_s = 26068433   # start 0004
ch_e = 12         # end
bp_e = 66232349   # end   0010
# 例3: circosLast:
ch_s = 7          # start
bp_s = 26241365   # start 0004
ch_e = 15         # end
bp_e = 40854180   # end   0006

def change(ch, bp):
    ch = str(ch)
    node = str(int(bp / node_size)).zfill(4)
    return ch, node
[ch_s, node_s] = change(ch_s, bp_s)
[ch_e, node_e] = change(ch_e, bp_e)

# paplotと同じ書式で出力
def change_data_format(idx):
    sample = org_df.at[idx, 'Sample']
    chr1 = org_df.at[idx, 'Chr1']
    bp1 = org_df.at[idx, 'Break1']
    chr2 = org_df.at[idx, 'Chr2']
    bp2 = org_df.at[idx, 'Break2']
    if chr1 == 'X':
        chr1 = '23'
    elif chr1 == 'Y':
        chr1 = '24'
    if chr2 == 'X':
        chr2 = '23'
    elif chr2 == 'Y':
        chr2 = '24'
    return f'{sample}__{chr1}__{bp1}__{chr2}__{bp2}\n'

# startとendをつなぐbp
print(f'# startとendをつなぐ線: ch_s(={ch_s}),node_s(={node_s}),bp_s(={bp_s}),ch_e(={ch_e}),node_e(={node_e}),bp_e(={bp_e})')
res = df[((ch_s == df['Chr1']) & (node_s == df['Node1']) & (ch_e == df['Chr2']) & (node_e == df['Node2'])) |
         ((ch_s == df['Chr2']) & (node_s == df['Node2']) & (ch_e == df['Chr1']) & (node_e == df['Node1']))]
p(res, '結果')
p(len(res), '個数')
with open(f'result/pd_b_{bp_s}', 'w') as f:
    for idx, _ in res.iterrows():
        f.write(change_data_format(idx))
print('#==============================')

# start/end から開始するbp
for i in range(2):
    if i == 0:
        [ch, node, bp, f] = [ch_s, node_s, bp_s, 's']
    elif i == 1:
        [ch, node, bp, f] = [ch_e, node_e, bp_e, 'e']
    print('# ch, node:', ch, node)

    print(f'# start/endからスタートする線: ch(={ch}),node(={node}),bp(={bp})')
    res = df[(ch == df['Chr1']) & (node == df['Node1']) |
             (ch == df['Chr2']) & (node == df['Node2'])]
    p(res, '結果')
    p(len(res), '個数')
    with open(f'result/pd_{f}_{bp_s}', 'w') as f:
        for idx, _ in res.iterrows():
            f.write(change_data_format(idx))
    print('#==============================')

# =============================================================================
# =============================================================================
# =============================================================================

print('')
print('# 同じbpをもつデータを確認する')
break1 = df['Break1'].tolist()
break2 = df['Break2'].tolist()
add_idx = 2  # ヘッダ行と1から行数を開始するために+2
with open('result/same_bp', 'w') as f:
    for idx1 in range(len(break1) - 1):
        bp_s = break1[idx1]
        bp_e = break2[idx1]
        for idx2 in range(idx1 + 1, len(break1)):
            bp2_s = break1[idx2]
            bp2_e = break2[idx2]
            if((bp_s == bp2_s and bp_e == bp2_e) or (bp_s == bp2_e and bp_e == bp2_s)):
                i1 = df.index[idx1]
                i2 = df.index[idx2]
                data = ' '.join([str(i1 + add_idx), str(i2 + add_idx),
                                 ','.join(map(str, ncan_df.iloc[i1].values)),
                                 ','.join(map(str, ncan_df.iloc[i2].values))])
                print(data)
                f.write(data + '\n')
