# 宇宙天気（日ごとの Ap と F10.7）を GFZ から取り、アプリが読む小さな CSV にする。
# OP's LAB Maps の軌道履歴グラフ（下段）用。以前は CelesTrak の SW-All.csv を読んでいたが、
# 利用者の回線から CelesTrak に届かないため、元データの GFZ を GitHub 経由で配る。
# GFZ はブラウザからの読み込みを許可するヘッダー（CORS）を返さないので、直接は読めない。
#
# 出典：GFZ Helmholtz Centre for Geosciences（CC BY 4.0）
#   Matzka et al. (2021) https://doi.org/10.1029/2020SW002641
#   F10.7 は Dominion Radio Astrophysical Observatory / Natural Resources Canada
# 黒点数（SN）は CC BY-NC 4.0 なので出力に含めない。
import sys
import urllib.request

SRC = 'https://kp.gfz.de/app/files/Kp_ap_Ap_SN_F107_since_1932.txt'
OUT = 'data/spaceweather.csv'

req = urllib.request.Request(SRC, headers={'User-Agent': 'Get_TLE (github.com/iqpslover-byte/Get_TLE)'})
with urllib.request.urlopen(req, timeout=120) as r:
    txt = r.read().decode('utf-8', 'replace')

rows = []
for line in txt.splitlines():
    if not line or line.startswith('#'):
        continue
    c = line.split()
    if len(c) < 28:
        continue
    y, m, d = int(c[0]), int(c[1]), int(c[2])
    if y < 1957:                      # CelesTrak の SW-All.csv と同じく 1957 年から
        continue
    ap = int(c[23])                   # 日平均 Ap（欠測は -1）
    f107 = float(c[25])               # F10.7 観測値（欠測は -1.0）
    rows.append('%04d-%02d-%02d,%s,%s' % (
        y, m, d,
        str(ap) if ap >= 0 else '',
        ('%.1f' % f107) if f107 >= 0 else ''))

# 取得が壊れていたら前回のファイルを残す（空や途中切れで上書きしない）
if len(rows) < 20000:
    print(f'✗ 宇宙天気の行数が少なすぎます（{len(rows)} 行）。更新しません')
    sys.exit(1)

# 列名はアプリの読み取り（CelesTrak の SW-All.csv と同じ名前）に合わせる
with open(OUT, 'w', encoding='utf-8', newline='\n') as f:
    f.write('DATE,AP_AVG,F10.7_OBS\n')
    f.write('\n'.join(rows) + '\n')
print(f'✓ 宇宙天気 {len(rows)} 日ぶんを保存しました（{rows[0][:10]}〜{rows[-1][:10]}）')
