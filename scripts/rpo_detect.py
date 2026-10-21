# キラー衛星（近接運用 RPO）の検知。OP's LAB Maps の知らせ用に data/rpo.json を出す。
# ・全物体どうしの総当たりはしない。「同じ軌道面・同じ高さ」の組だけに絞ってから2段で判定する。
#   ① 警報（張り付き）＝この先3日のうち 50km 以内に12時間以上いる
#   ② 注意（前触れ）  ＝面 0.3°以内・平均高度の差 1km 以内・1周平均の距離 1,000km 以内・
#                       1日で 10km 以上縮む・どちらかが14日以内に大きく噴射（前のTLEで進めた位置と 20km 超ずれる）
# ・答え合わせは OP's LAB Maps の tools/rpo（Cosmos 2612×ICEYE-X36・2026-10-18）。
# ・噴射を見るため、候補の機の TLE を data/rpo_state.json に控える（新しく候補に入った機だけ履歴を取る）。
# ・Space-Track の目安は全カタログの問い合わせが1時間に1回まで＝6時間ごとのワークフローで1回だけ呼ぶ。
# ・静止軌道は面では絞れない（経度の近さの別判定が要る）ので、今は低軌道（平均高度 2,000km 未満）だけ。
import bisect
import json
import os
import re
import sys
import time
from datetime import datetime, timedelta, timezone

import numpy as np
from sgp4.api import Satrec, SatrecArray, jday

BASE = 'https://www.space-track.org'
OUT = 'data/rpo.json'
STATE = 'data/rpo_state.json'
RE = 6378.135

PLANE = 0.3          # 軌道面の角度差（°）
MARGIN = 10.0        # 近地点〜遠地点の範囲の重なりに足す余裕（km）
LEO = 2000.0         # これより低い平均高度だけ見る（km）
# ①
DWELL_KM, DWELL_H, DWELL_DAYS, STEP = 50.0, 12.0, 3, 30
# ②
POOL_DALT = 5.0      # 控えを持つ候補＝平均高度の差がこれ未満の組（km）
DALT, SEP, CLOSE, MAN, MAN_DAYS = 1.0, 1000.0, 10.0, 20.0, 14
KEEP_TLE_DAYS = 3    # 1日前の状態を出すのに要る長さ＋余裕
SINCE_DAYS = 30      # 開始の日時をさかのぼる上限（これより前から続いていれば「30日以上前から」）
NEW_DAYS = 30        # 打上げからこの日数以内の機は外す（分離直後の機どうしが寄って見えるだけ）
# ステーションと結合・ドッキングする機（名前の頭）。結合中は張り付いて当然なので外す
STN = ('ISS', 'CSS', 'TIANGONG', 'SZ-', 'TIANZHOU', 'PROGRESS', 'SOYUZ', 'DRAGON', 'CYGNUS', 'HTV',
       'POISK', 'STARLINER')
# いつも並んで飛んでいると分かっている組。一覧の下へ回し、知らせは鳴らさない（消さないのは、崩れた時に気づけるように）
KNOWN = {('31698', '36605'): 'TerraSAR-X と TanDEM-X の編隊飛行'}
HIST = 'data/rpo_hist.json'   # 鳴った組の履歴（消さずに貯める・アプリの ARCHIVE）
NEW_MD = '_new_rpo.md'   # 新しい組が鳴った時だけ書く＝ワークフローが Issue を立てる（コミットしない）
# 軌道面が近い衛星（アプリ v3.63）用に、衛星（PAYLOAD）全機の軌道面と TLE を書き出す。約4MB を6時間ごとに main へ積むと
# 履歴が膨らむので、ワークフローが履歴を持たない planes ブランチへ上書きで置く（main にはコミットしない）
PLANES = '_planes.json'
FIELDS ='NORAD_CAT_ID,OBJECT_NAME,OBJECT_ID,OBJECT_TYPE,COUNTRY_CODE,LAUNCH_DATE,EPOCH,SEMIMAJOR_AXIS,' \
         'PERIAPSIS,APOAPSIS,TLE_LINE1,TLE_LINE2,INCLINATION,RA_OF_ASC_NODE,ECCENTRICITY'


def utc(s):
    return datetime.fromisoformat(s.replace('Z', '')).replace(tzinfo=timezone.utc)


def iso(t):
    return t.strftime('%Y-%m-%dT%H:%M:%SZ')


class ST:
    def __init__(self):
        import requests
        self.s = requests.Session()
        self.s.verify = os.environ.get('ST_NOVERIFY') != '1'   # 手元の試験でSSLを傍受される回線用
        r = self.s.post(BASE + '/ajaxauth/login',
                        data={'identity': os.environ['SPACETRACK_USERNAME'],
                              'password': os.environ['SPACETRACK_PASSWORD']}, timeout=60)
        r.raise_for_status()

    def q(self, path):
        r = self.s.get(BASE + '/basicspacedata/query' + path, timeout=600)
        r.raise_for_status()
        return r.json()


# ── 候補の組（同じ面・高さが重なる・別の打上げ・名前の系列が違う） ──
def fam(n):
    m = re.match(r'[A-Z]+', n.upper())
    return m.group(0) if m else n.upper()


def is_stn(n):
    n = n.upper()
    return any(n.startswith(s) for s in STN)


def candidates(objs, t0):
    sats = [Satrec.twoline2rv(x['TLE_LINE1'], x['TLE_LINE2']) for x in objs]
    j, f = jday(t0.year, t0.month, t0.day, t0.hour, t0.minute, t0.second)
    e, r, v = SatrecArray(sats).sgp4(np.array([j]), np.array([f]))
    r, v, ok = r[:, 0], v[:, 0], e[:, 0] == 0
    h = np.cross(r, v)
    h /= np.linalg.norm(h, axis=1)[:, None]
    per = np.array([float(x['PERIAPSIS']) for x in objs])
    apo = np.array([float(x['APOAPSIS']) for x in objs])
    fams = [fam(x['OBJECT_NAME']) for x in objs]
    lns = [(x.get('OBJECT_ID') or '')[:8] for x in objs]
    cos_lim = np.cos(np.radians(PLANE))
    idx = [i for i in np.argsort(per) if ok[i]]
    pairs = []
    for p, i in enumerate(idx):
        for k in idx[p + 1:]:
            if per[k] > apo[i] + MARGIN:
                break
            if per[i] > apo[k] + MARGIN or lns[i] == lns[k] or fams[i] == fams[k]:
                continue
            if h[i] @ h[k] >= cos_lim:
                pairs.append((i, k))
    return pairs


# ── ① 張り付き ──
def dwell(objs, pairs, t0):
    ids = sorted({i for p in pairs for i in p})
    pos = {k: n for n, k in enumerate(ids)}
    arr = SatrecArray([Satrec.twoline2rv(objs[k]['TLE_LINE1'], objs[k]['TLE_LINE2']) for k in ids])
    d = np.empty((len(pairs), DWELL_DAYS * 86400 // STEP), dtype=np.float32)
    a_i = np.array([pos[a] for a, b in pairs])
    b_i = np.array([pos[b] for a, b in pairs])
    n = 0
    for c in range(DWELL_DAYS * 4):           # 6時間ずつ（全部を一度に持つと数GBになる）
        ts = [t0 + timedelta(seconds=s) for s in range(c * 21600, (c + 1) * 21600, STEP)]
        jj, ff = zip(*[jday(t.year, t.month, t.day, t.hour, t.minute, t.second) for t in ts])
        _, r, _ = arr.sgp4(np.array(jj), np.array(ff))
        d[:, n:n + len(ts)] = np.linalg.norm(r[a_i] - r[b_i], axis=2)
        n += len(ts)
    out = []
    for p, (a, b) in enumerate(pairs):
        hrs = float((d[p] < DWELL_KM).sum()) * STEP / 3600
        if hrs < DWELL_H:
            continue
        m = int(np.argmin(d[p]))
        out.append({'a': a, 'b': b, 'hours': round(hrs, 1), 'nowKm': round(float(d[p, 0]), 1),
                    'minKm': round(float(d[p, m]), 2), 'minAt': iso(t0 + timedelta(seconds=m * STEP)),
                    'medKm': round(float(np.median(d[p])), 1)})
    return out


# ── 候補の機の TLE の控え（噴射の検出用） ──
def rv(sat, ts):
    j, f = zip(*[jday(t.year, t.month, t.day, t.hour, t.minute, t.second + t.microsecond / 1e6) for t in ts])
    _, r, v = sat.sgp4_array(np.array(j), np.array(f))
    return r, v


def jump(prev, cur):
    """前のTLEで次の元期へ進めた位置と次のTLEの差（km）"""
    t1 = utc(cur[0])
    if (t1 - utc(prev[0])).total_seconds() < 3600:
        return 0.0
    r0, _ = rv(Satrec.twoline2rv(prev[1], prev[2]), [t1])
    r1, _ = rv(Satrec.twoline2rv(cur[1], cur[2]), [t1])
    return float(np.linalg.norm(r0[0] - r1[0]))


def add_tles(rec, rows, now):
    """rows＝[(元期, 1行目, 2行目, 平均高度)] を元期順に足し、噴射を拾い、古いものを捨てる"""
    for row in sorted(rows):
        if rec['tles'] and row[0] <= rec['tles'][-1][0]:
            continue
        if rec['tles'] and jump(rec['tles'][-1], row) > MAN:
            rec['man'].append(row[0])
        rec['tles'].append(list(row))
    cut = iso(now - timedelta(days=KEEP_TLE_DAYS))
    keep = [x for x in rec['tles'] if x[0] >= cut]
    rec['tles'] = keep if keep else rec['tles'][-1:]
    cut = iso(now - timedelta(days=MAN_DAYS))
    rec['man'] = [m for m in rec['man'] if m >= cut]


def row_of(x):
    return (x['EPOCH'][:19] + 'Z', x['TLE_LINE1'], x['TLE_LINE2'], round(float(x['SEMIMAJOR_AXIS']) - RE, 3))


def latest(rec, t):
    es = [x[0] for x in rec['tles']]
    i = bisect.bisect_right(es, iso(t)) - 1
    return rec['tles'][i] if i >= 0 else None


ORB = [timedelta(minutes=m) for m in range(0, 96, 8)]


def pstate(A, B, t):
    if not A or not B:
        return None
    ts = [t + o for o in ORB]
    ra, va = rv(Satrec.twoline2rv(A[1], A[2]), ts)
    rb, vb = rv(Satrec.twoline2rv(B[1], B[2]), ts)
    ha, hb = np.cross(ra[0], va[0]), np.cross(rb[0], vb[0])
    pl = float(np.degrees(np.arccos(np.clip(ha @ hb / np.linalg.norm(ha) / np.linalg.norm(hb), -1, 1))))
    return pl, abs(A[3] - B[3]), float(np.linalg.norm(ra - rb, axis=1).mean())


def approach(pool, sats, t):
    """② 注意。pool＝[(番号A, 番号B)]・sats＝控え"""
    out = []
    for a, b in pool:
        if a not in sats or b not in sats:
            continue
        s = pstate(latest(sats[a], t), latest(sats[b], t), t)
        if not s or s[0] > PLANE or s[1] > DALT or s[2] > SEP:
            continue
        t1 = t - timedelta(days=1)
        s1 = pstate(latest(sats[a], t1), latest(sats[b], t1), t1)
        if not s1 or s1[2] - s[2] < CLOSE:
            continue
        cut = iso(t - timedelta(days=MAN_DAYS))
        ma = [m for m in sats[a]['man'] if cut <= m <= iso(t)]
        mb = [m for m in sats[b]['man'] if cut <= m <= iso(t)]
        if not (ma or mb):
            continue
        out.append({'a': a, 'b': b, 'planeDeg': round(s[0], 3), 'dAltKm': round(s[1], 2),
                    'sepKm': round(s[2], 1), 'closeKmPerDay': round(s1[2] - s[2], 1), 'manA': ma, 'manB': mb})
    return out


def info(x):
    return {'id': x['NORAD_CAT_ID'], 'name': x['OBJECT_NAME'], 'intl': x.get('OBJECT_ID') or '',
            'country': x.get('COUNTRY_CODE') or '', 'altKm': round(float(x['SEMIMAJOR_AXIS']) - RE, 1),
            # 地図に描くため（大半はアプリの tle_recent.json に居ない）
            'epoch': x['EPOCH'][:19] + 'Z', 'tle1': x['TLE_LINE1'], 'tle2': x['TLE_LINE2']}


def jst(t):
    """'2026-10-16T14:14:11Z' → '10-16 23:14'（JST）。噴射の時刻は見つけた TLE の元期＝実際はその前の TLE との間のどこか"""
    return (utc(t) + timedelta(hours=9)).strftime('%m-%d %H:%M')


def notify(new, now):
    """新しく鳴った組の知らせ文。1行目が Issue の題名"""
    kinds = {'dwell': '① 警報（張り付き）', 'approach': '② 注意（前触れ）'}
    L = ['キラー衛星の検知：' + '、'.join(f"{o['a']['name']} × {o['b']['name']}" for o in new), '']
    for o in new:
        L.append(f"### {kinds[o['kind']]}　{o['a']['name']}（{o['a']['country']}）× {o['b']['name']}（{o['b']['country']}）")
        if o['kind'] == 'dwell':
            L.append(f"- この先3日のうち {DWELL_KM:.0f}km 以内に {o['hours']} 時間・最小 {o['minKm']}km・今 {o['nowKm']}km")
        else:
            L.append(f"- 1周平均の距離 {o['sepKm']}km・1日で {o['closeKmPerDay']}km 縮む・面の差 {o['planeDeg']}°・高度差 {o['dAltKm']}km")
        for side in ('A', 'B'):
            if o.get('man' + side):
                L.append(f"- 噴射 {o[side.lower()]['name']}　" + '、'.join(jst(m) for m in o['man' + side]) + ' JST')
        L.append(f"- 番号 {o['a']['id']} / {o['b']['id']}・平均高度 {o['a']['altKm']}km")
        L.append('')
    L.append(f"判定 {(now + timedelta(hours=9)).strftime('%Y-%m-%d %H:%M')} JST")
    with open(NEW_MD, 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')


# ── 開始の日時（その組が ① か ② の条件に入り始めた時） ──
def dwell_at(A, B, t):
    """その時点で最新の TLE で ① を判定し直す（本番より粗い60秒刻み）"""
    if not A or not B:
        return False
    ts = [t + timedelta(seconds=x) for x in range(0, DWELL_DAYS * 86400, 60)]
    ra, _ = rv(Satrec.twoline2rv(A[1], A[2]), ts)
    rb, _ = rv(Satrec.twoline2rv(B[1], B[2]), ts)
    return float((np.linalg.norm(ra - rb, axis=1) < DWELL_KM).sum()) * 60 / 3600 >= DWELL_H


def since_of(pairs, st, now):
    """pairs＝[(番号A, 番号B)]。6時間ずつさかのぼり、近くに居続けていた最初の時刻を返す
    （その時点で最新だった TLE で判定し直す）。{(A, B): (時刻, 上限で打ち切ったか)}
    「近くに居る」＝① か、同じ面（PLANE）・高度差 POOL_DALT 以内・1周平均の距離 SEP 以内。
    ★② そのもの（高度差1km）で続きを見ると、最後の詰めの噴射で高度差が数km開いた間が抜け落ちる
      （Cosmos 2612：10/15 18Z〜10/16 12Z の24時間＝高度差2.4kmで距離180→11km）。"""
    ids = sorted({n for p in pairs for n in p}, key=int)
    frm = (now - timedelta(days=SINCE_DAYS + MAN_DAYS + 1)).strftime('%Y-%m-%d')
    rows = st.q('/class/gp_history/NORAD_CAT_ID/' + ','.join(ids) + '/EPOCH/%3E' + frm +
                '/predicates/NORAD_CAT_ID,EPOCH,SEMIMAJOR_AXIS,TLE_LINE1,TLE_LINE2/orderby/EPOCH%20asc/format/json')
    rec = {k: {'tles': [], 'man': []} for k in ids}
    got = {}
    for x in rows:
        got.setdefault(x['NORAD_CAT_ID'], set()).add(row_of(x))
    for k in ids:
        for row in sorted(got.get(k, [])):   # 噴射を拾いながら全部並べる（古いものは捨てない）
            r = rec[k]
            if r['tles'] and jump(r['tles'][-1], row) > MAN:
                r['man'].append(row[0])
            r['tles'].append(list(row))
    out = {}
    for a, b in pairs:
        t, ok = now.replace(minute=0, second=0), None
        for _ in range(SINCE_DAYS * 4 + 1):
            A, B = latest(rec[a], t), latest(rec[b], t)
            p = pstate(A, B, t)
            near = p and p[0] <= PLANE and p[1] <= POOL_DALT and p[2] <= SEP
            if not (near or dwell_at(A, B, t)):
                break
            ok = t
            t -= timedelta(hours=6)
        else:
            out[(a, b)] = (iso(ok), True)
            continue
        out[(a, b)] = (iso(ok or now), False)
    return out


def main():
    now = datetime.now(timezone.utc).replace(microsecond=0)
    st = ST()
    G = st.q('/class/gp/decay_date/null-val/epoch/%3Enow-30/OBJECT_TYPE/PAYLOAD,UNKNOWN'
             '/predicates/' + FIELDS + '/orderby/norad_cat_id/format/json')
    print(f'全カタログ {len(G)} 件')
    write_planes(G, now)
    run(G, now, st)


def write_planes(G, now):
    """1行 = [番号, 名前, 軌道傾斜角°, RAAN°, 平均高度km, 離心率, 元期(UNIX秒), TLE1, TLE2]"""
    out = []
    for x in G:
        if x.get('OBJECT_TYPE') != 'PAYLOAD' or not x.get('TLE_LINE1') or not x.get('TLE_LINE2'):
            continue
        try:
            out.append([int(x['NORAD_CAT_ID']), x['OBJECT_NAME'] or '', round(float(x['INCLINATION']), 4),
                        round(float(x['RA_OF_ASC_NODE']), 4), round(float(x['SEMIMAJOR_AXIS']) - RE, 1),
                        round(float(x['ECCENTRICITY']), 6), int(utc(x['EPOCH']).timestamp()), x['TLE_LINE1'], x['TLE_LINE2']])
        except Exception:
            pass
    with open(PLANES, 'w', encoding='utf-8') as f:
        json.dump({'updated_at': int(now.timestamp()), 'sats': out}, f, separators=(',', ':'), ensure_ascii=False)
    print(f'✓ {PLANES} {len(out)} 機・{os.path.getsize(PLANES) // 1024} KB')


def run(G, now, st=None):
    t0 = now.replace(minute=0, second=0)
    new_cut = (now - timedelta(days=NEW_DAYS)).strftime('%Y-%m-%d')
    objs = [x for x in G if x.get('TLE_LINE1') and x.get('TLE_LINE2')
            and (float(x['PERIAPSIS']) + float(x['APOAPSIS'])) / 2 < LEO
            and not is_stn(x['OBJECT_NAME']) and not x['OBJECT_NAME'].upper().startswith('TBA')
            and (x.get('LAUNCH_DATE') or '0000') < new_cut]
    by = {x['NORAD_CAT_ID']: x for x in objs}
    pairs = candidates(objs, t0)
    print(f'低軌道 {len(objs)} 機・同じ面の組 {len(pairs)}')

    alerts = dwell(objs, pairs, t0)
    for o in alerts:
        o['a'], o['b'] = objs[o['a']]['NORAD_CAT_ID'], objs[o['b']]['NORAD_CAT_ID']
    print(f'① 警報 {len(alerts)} 組')

    # ② 候補と控え
    pool = []
    for i, k in pairs:
        A, B = objs[i], objs[k]
        if abs(float(A['SEMIMAJOR_AXIS']) - float(B['SEMIMAJOR_AXIS'])) < POOL_DALT:
            pool.append((A['NORAD_CAT_ID'], B['NORAD_CAT_ID']))
    try:
        state = json.load(open(STATE, encoding='utf-8'))
    except Exception:
        state = {'sats': {}, 'seen': {}}
    sats = state['sats']
    # ① の組の2機も控える（噴射の時刻を載せるため）
    need = sorted({n for p in pool for n in p} | {n for o in alerts for n in (o['a'], o['b'])}, key=int)
    fresh = [n for n in need if n not in sats]
    print(f'② 候補 {len(pool)} 組・{len(need)} 機（新しく控える {len(fresh)} 機）')
    if fresh and st:
        since = (now - timedelta(days=MAN_DAYS + 1)).strftime('%Y-%m-%d')
        for n in range(0, len(fresh), 100):
            rows = st.q('/class/gp_history/NORAD_CAT_ID/' + ','.join(fresh[n:n + 100]) +
                        '/EPOCH/%3E' + since + '/predicates/NORAD_CAT_ID,EPOCH,SEMIMAJOR_AXIS,TLE_LINE1,TLE_LINE2'
                        '/orderby/EPOCH%20asc/format/json')
            got = {}
            for x in rows:
                got.setdefault(x['NORAD_CAT_ID'], []).append(row_of(x))
            for k in fresh[n:n + 100]:
                sats[k] = {'tles': [], 'man': []}
                add_tles(sats[k], got.get(k, []), now)
            time.sleep(3)
    for k in need:
        sats.setdefault(k, {'tles': [], 'man': []})
        add_tles(sats[k], [row_of(by[k])], now)
    for k in list(sats):
        if k not in by or k not in need:   # 候補から外れた機は捨てる（戻ってきたら取り直す）
            del sats[k]

    for o in alerts:
        o['manA'], o['manB'] = list(sats[o['a']]['man']), list(sats[o['b']]['man'])
    warns = approach(pool, sats, now)
    print(f'② 注意 {len(warns)} 組')

    # 最初に鳴った時刻を覚える（同じ組が鳴り続けても知らせは1回で済むように）
    seen = state['seen']
    items, new = [], []
    for kind, arr in (('dwell', alerts), ('approach', warns)):
        for o in arr:
            ab = tuple(sorted((o['a'], o['b']), key=int))
            key = kind + ':' + '-'.join(ab)
            first = seen.get(key, {}).get('first', iso(now))
            seen[key] = {**seen.get(key, {}), 'first': first, 'last': iso(now)}
            a, b = o.pop('a'), o.pop('b')
            it = {'kind': kind, 'first': first, 'a': info(by[a]), 'b': info(by[b]), **o, '_key': key, '_ab': ab}
            if ab in KNOWN:
                it['known'] = KNOWN[ab]
            elif first == iso(now):
                new.append(it)
            items.append(it)
    # 開始の日時は、まだ出していない組だけ履歴から求める（既知の組は求めない＝ずっと続いている）
    todo = [it for it in items if 'known' not in it and 'since' not in seen[it['_key']]]
    if todo and st:
        try:
            got = since_of(list({it['_ab'] for it in todo}), st, now)
            for it in todo:
                seen[it['_key']]['since'], seen[it['_key']]['sinceCapped'] = got[it['_ab']]
        except Exception as e:
            print(f'! 開始の日時を出せませんでした: {e}')
    hist_update(items, now)
    for it in items:
        k = it.pop('_key'); it.pop('_ab')
        if 'since' in seen[k]:
            it['since'] = seen[k]['since']
            if seen[k].get('sinceCapped'):
                it['sinceCapped'] = True
    cut = iso(now - timedelta(days=60))
    state['seen'] = {k: v for k, v in seen.items() if v['last'] >= cut}
    items.sort(key=lambda o: ('known' in o, o['kind'] != 'dwell', o.get('minKm', o.get('sepKm', 0))))
    if new:
        notify(new, now)
        print(f'新しく鳴った組 {len(new)}（{NEW_MD}）')

    os.makedirs('data', exist_ok=True)
    with open(OUT, 'w', encoding='utf-8') as f:
        json.dump({'updated_at': iso(now), 'count': len(items), 'items': items},
                  f, separators=(',', ':'), ensure_ascii=False)
    with open(STATE, 'w', encoding='utf-8') as f:
        json.dump(state, f, separators=(',', ':'))
    print(f'✓ {OUT} {len(items)} 件・{STATE} {os.path.getsize(STATE) // 1024} KB')


def hist_add(H, it, t):
    """1回の検知の1組を履歴へ。H＝{key: 記録}・t＝その回の時刻(ISO)"""
    key = it['_key']
    d = it.get('minKm', it.get('sepKm'))
    h = H.get(key)
    if not h:
        h = H[key] = {'key': key, 'kind': it['kind'], 'first': it.get('first', t), 'last': t, 'minKm': d,
                      'a': {k: it['a'][k] for k in ('id', 'name', 'intl', 'country') if k in it['a']},
                      'b': {k: it['b'][k] for k in ('id', 'name', 'intl', 'country') if k in it['b']}}
    h['last'] = max(h['last'], t)
    h['first'] = min(h['first'], it.get('first', t))
    if d is not None and (h.get('minKm') is None or d < h['minKm']):
        h['minKm'] = d
    for k in ('since', 'known'):
        if it.get(k):
            h[k] = it[k]
    for k in ('a', 'b'):   # 名前は最新へ（TBA などから付け直されることがある）
        h[k].update({x: it[k][x] for x in ('name', 'intl', 'country') if x in it[k]})


def hist_update(items, now):
    """鳴った組の履歴 data/rpo_hist.json を更新する。いま鳴っている組は active"""
    try:
        H = {h['key']: h for h in json.load(open(HIST, encoding='utf-8'))['items']}
    except Exception:
        H = {}
    t = iso(now)
    for it in items:
        hist_add(H, it, t)
    keys = {it['_key'] for it in items}
    for h in H.values():
        h['active'] = h['key'] in keys
    out = sorted(H.values(), key=lambda h: h['last'], reverse=True)
    with open(HIST, 'w', encoding='utf-8') as f:
        json.dump({'updated_at': t, 'count': len(out), 'items': out}, f, separators=(',', ':'), ensure_ascii=False)


if __name__ == '__main__':
    sys.exit(main())
