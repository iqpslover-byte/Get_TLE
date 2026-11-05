# 「狙っている機」の検知（rpo_detect.py から呼ぶ）。距離は見ない。
# 噴射した機について、噴射の前と後で「ある1機への近づき方」がどう変わったかで判定する。
#   plane  ＝噴射で傾斜角か高度を変え、この先ある1機の軌道面に重なるようになった（何か月先でも）
#            例：Cosmos 2611・2613（2026-05 に傾斜角 +0.87° → ICEYE-X37 の軌道面へ約半年で重なる）
#   phase  ＝同じ軌道面で高度を変え、ある1機との位相（軌道上の前後）が詰まる向きになった（先回り・待ち伏せ）
#            例：Cosmos 2612（2026-09-30 に 590km へ上げ、10-12 に 543km へ下げて ICEYE-X36 に寄せた）
#   keep   ＝ある1機と同じ軌道面・同じ高度にいて、噴射で保っている（2回以上）
#   launch ＝打上げ後に初めて出てきた時点で、ある1機の軌道面・高度に合っている
# 噴射＝控え（6時間ごと）との比較で、傾斜角・RAAN（歳差を除いた残り）・高度が自然な変化より大きく変わったこと。
#   弱い推力で少しずつ変える機のため、約4週間前の控えとも比べる（slow）。
# ふつうの軌道保持（太陽同期の傾斜角の直しなど）で鳴りすぎないよう、狙いが数機に絞れる時だけ候補にする。
# 答え合わせ：OP's LAB Maps の tools/rpo（2026-11-05・Cosmos 2611／2612／2613 × ICEYE-X36／X37）。
import math
import numpy as np

RE = 6378.135
MU = 398600.4418
J2 = 1.08262668e-3
DAY = 86400.0

# 噴射の見分け（隣の控えとの差）
MAN_DI = 0.02        # 傾斜角（°）＝TLE のばらつきは 0.003° 程度・太陽同期の自然な変化は1週間で 0.001° 程度
MAN_DW = 0.03        # RAAN の歳差を除いた残り（°）
MAN_UP = 0.5         # 高度が上がった（km）＝大気抵抗では上がらない
MAN_DOWN = 2.0       # 3日以内に下がった（km）＝大気抵抗より速い
SLOW_DI = 0.05       # 約4週間で傾斜角がこれ以上変わった（°）＝弱い推力
SLOW_DAYS = 24
# 狙いの判定
HORIZON = 365        # 軌道面が重なる日を探す先（日）
PL_HIT = 0.1         # 重なったとみなす軌道面の差（°）
PL_DI = 0.1          # 傾斜角の差（°）
PL_DALT = 5.0        # 平均高度の差（km）
PH_PLANE = 0.3       # 位相を見る組の軌道面の差（°）
PH_DALT = 60.0       # 先回り・待ち伏せの高度差の上限（km）
PH_DAYS = 45         # 位相がそろうまでの日数の上限
KEEP_PLANE, KEEP_DALT, KEEP_HITS = 0.15, 5.0, 2
LAUNCH_DAYS = 45     # 打上げからこの日数以内に初めて出てきた機を launch で見る
LAUNCH_DI, LAUNCH_DALT = 0.05, 3.0
MAX_TARGETS = 3      # 狙いの候補がこれより多い噴射はふつうの軌道保持とみなす
STN = ('ISS', 'CSS', 'TIANGONG', 'SZ-', 'TIANZHOU', 'PROGRESS', 'SOYUZ', 'DRAGON', 'CYGNUS', 'HTV',
       'POISK', 'STARLINER')


def rates(sma, ecc, inc_deg):
    """J2 による RAAN と緯度引数（ω+M）の進み（°/日）"""
    n = np.sqrt(MU / sma ** 3)                     # rad/s
    i = np.radians(inc_deg)
    p = (RE / (sma * (1 - ecc ** 2))) ** 2
    dW = -1.5 * n * J2 * p * np.cos(i)
    dw = 0.75 * n * J2 * p * (5 * np.cos(i) ** 2 - 1)
    dM = n * (1 + 0.75 * J2 * p * np.sqrt(1 - ecc ** 2) * (3 * np.cos(i) ** 2 - 1))
    return np.degrees(dW) * DAY, np.degrees(dw + dM) * DAY


def el_of(x):
    """Space-Track の1行 → 控えの1件 [元期(秒), 傾斜角, RAAN, 平均高度, 離心率, 緯度引数]"""
    from datetime import datetime, timezone
    t = datetime.fromisoformat(x['EPOCH'].replace('Z', '')).replace(tzinfo=timezone.utc).timestamp()
    return [int(t), round(float(x['INCLINATION']), 4), round(float(x['RA_OF_ASC_NODE']), 4),
            round(float(x['SEMIMAJOR_AXIS']) - RE, 3), round(float(x['ECCENTRICITY']), 6),
            round((float(x['ARG_OF_PERICENTER']) + float(x['MEAN_ANOMALY'])) % 360, 3)]


def at(el, t):
    """控え el を時刻 t（秒）へ進めた RAAN と緯度引数（°）"""
    dW, du = rates(el[3] + RE, el[4], el[1])
    d = (t - el[0]) / DAY
    return (el[2] + dW * d) % 360, (el[5] + du * d) % 360


def maneuver(prev, cur):
    """隣の控えどうし。噴射なら {'di','dw','dalt'}（変わった量）、無ければ None"""
    if cur[0] - prev[0] < 3 * 3600:
        return None
    w, _ = at(prev, cur[0])
    di = cur[1] - prev[1]
    dw = (cur[2] - w + 180) % 360 - 180
    dalt = cur[3] - prev[3]
    days = (cur[0] - prev[0]) / DAY
    if abs(di) >= MAN_DI or abs(dw) >= MAN_DW or dalt >= MAN_UP or (dalt <= -MAN_DOWN and days <= 3):
        return {'di': round(di, 4), 'dw': round(dw, 4), 'dalt': round(dalt, 2)}
    return None


class Cat:
    """狙われる側の全カタログ（衛星だけ）。配列で持って1回の噴射ごとに全機と比べる"""
    def __init__(self, objs, now):
        ok = [x for x in objs if x.get('OBJECT_TYPE') == 'PAYLOAD' and not any(x['OBJECT_NAME'].upper().startswith(s) for s in STN)]
        self.x = ok
        self.id = [x['NORAD_CAT_ID'] for x in ok]
        E = np.array([el_of(x) for x in ok], dtype=float)
        self.inc, self.alt, self.ecc = E[:, 1], E[:, 3], E[:, 4]
        self.dW, self.du = rates(self.alt + RE, self.ecc, self.inc)
        d = (now - E[:, 0]) / DAY
        self.W0 = (E[:, 2] + self.dW * d) % 360       # now での RAAN
        self.u0 = (E[:, 5] + self.du * d) % 360
        self.now = now
        self.fam = [fam(x['OBJECT_NAME']) for x in ok]
        self.ln = [(x.get('OBJECT_ID') or '')[:8] for x in ok]


def fam(n):
    import re
    m = re.match(r'[A-Z]+', n.upper())
    return m.group(0) if m else n.upper()


def plane_deg(i1, W1, i2, W2):
    i1, i2 = np.radians(i1), np.radians(i2)
    c = np.cos(i1) * np.cos(i2) + np.sin(i1) * np.sin(i2) * np.cos(np.radians(W1 - W2))
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def forecast(el, C, mask):
    """el（狙う側）と C[mask] の、今から HORIZON 日の軌道面の差の最小と、その日（今からの日数）"""
    W, _ = at(el, C.now)
    dW, _ = rates(el[3] + RE, el[4], el[1])
    days = np.arange(0, HORIZON + 1, 1.0)
    Wx = W + dW * days                                    # (D,)
    Wy = C.W0[mask][:, None] + C.dW[mask][:, None] * days  # (N,D)
    p = plane_deg(el[1], Wx[None, :], C.inc[mask][:, None], Wy)
    k = np.argmin(p, axis=1)
    return p[np.arange(len(k)), k], k, p[:, 0]


def phase_days(el, C, idx):
    """位相がそろうまでの日数（詰まる向きでなければ inf）と今の位相の差（°・相手が前なら＋）"""
    _, ux = at(el, C.now)
    _, dux = rates(el[3] + RE, el[4], el[1])
    dphi = (C.u0[idx] - ux + 540) % 360 - 180            # 相手が何度前か
    rel = dux - C.du[idx]                                # 自分が1日に何度追いつくか
    with np.errstate(divide='ignore', invalid='ignore'):
        d = np.where(rel * dphi > 0, np.abs(dphi / rel), np.inf)
    return d, dphi


def scan(chaser, pre, post, C, kind_hint=None):
    """kind_hint＝'man'（跳び）／'slow'（弱い推力・位相は見ない＝高度の自然な下がりを噴射と取り違える）／'launch'"""
    """1回の噴射（pre→post）で狙われた候補。chaser＝Space-Track の1行（名前・打上げ）。[(相手の番号, 中身)]"""
    me = chaser['NORAD_CAT_ID']
    ln, fm = (chaser.get('OBJECT_ID') or '')[:8], fam(chaser['OBJECT_NAME'])
    base = np.array([(C.id[k] != me) and C.ln[k] != ln and C.fam[k] != fm for k in range(len(C.id))])
    out = []
    # plane：噴射の後、傾斜角と高度がそろい、この先軌道面が重なる。噴射の前は重ならなかった
    m = base & (np.abs(C.inc - post[1]) <= PL_DI) & (np.abs(C.alt - post[3]) <= PL_DALT)
    if m.any():
        mn, day, now_p = forecast(post, C, m)
        ok = (mn <= PL_HIT) & (now_p > mn + 0.02)          # まだ重なっていない（重なっている組は keep で見る）
        if pre is not None:
            mn0, _, _ = forecast(pre, C, m)
            ok &= mn0 >= np.maximum(0.3, 3 * mn)
        for k, a, d, p0 in zip(np.where(m)[0][ok], mn[ok], day[ok], now_p[ok]):
            out.append((C.id[k], {'mode': 'plane', 'minDeg': round(float(a), 3), 'inDays': int(d),
                                  'planeDeg': round(float(p0), 3), 'dInc': round(float(C.inc[k] - post[1]), 3),
                                  'dAltKm': round(float(C.alt[k] - post[3]), 2)}))
    # phase と keep：今すでに同じ軌道面の相手
    Wp, _ = at(post, C.now)
    pnow = plane_deg(post[1], Wp, C.inc, C.W0)
    m2 = base & (pnow <= PH_PLANE) & (np.abs(C.alt - post[3]) <= PH_DALT)
    if m2.any() and pre is not None and kind_hint == 'man':
        idx = np.where(m2)[0]
        d1, dphi = phase_days(post, C, idx)
        d0, _ = phase_days(pre, C, idx)
        for k, a, b, f in zip(idx, d1, d0, dphi):
            if a <= PH_DAYS and (not np.isfinite(b) or a < b / 2) and abs(post[3] - pre[3]) >= MAN_UP:
                out.append((C.id[k], {'mode': 'phase', 'phaseDays': round(float(a), 1), 'phaseDeg': round(float(f), 1),
                                      'planeDeg': round(float(pnow[k]), 3), 'dAltKm': round(float(C.alt[k] - post[3]), 2)}))
    m3 = base & (pnow <= KEEP_PLANE) & (np.abs(C.alt - post[3]) <= KEEP_DALT)
    if pre is not None:
        Wq, _ = at(pre, C.now)
        m3 &= plane_deg(pre[1], Wq, C.inc, C.W0) <= KEEP_PLANE
    for k in np.where(m3)[0]:
        if not any(t == C.id[k] for t, _ in out):
            out.append((C.id[k], {'mode': 'keep', 'planeDeg': round(float(pnow[k]), 3), 'dAltKm': round(float(C.alt[k] - post[3]), 2)}))
    if kind_hint == 'launch':
        out = [(t, o) for t, o in out if o['mode'] in ('plane', 'keep')
               and abs(o.get('dInc', 0)) <= LAUNCH_DI and abs(o['dAltKm']) <= LAUNCH_DALT]
        for t, o in out:
            o['mode'] = 'launch'
    # 狙いが絞れない（多くの機と同じくらい合う）噴射は、ふつうの軌道保持とみなして捨てる
    if len({t for t, _ in out}) > MAX_TARGETS:
        return []
    return out


# ── 6時間ごとの処理（rpo_detect.py から呼ぶ） ──
KEEP_EL_DAYS = 30    # 控えを残す日数
EL_GAP_DAYS = 3      # 控えは3日おき（＋いちばん新しい1件）＝大きさを抑える
ACTIVE_DAYS = 365    # 最後の噴射からこの日数を過ぎた組は外す（それまでは幾何だけで決める＝Cosmos 2611 は噴射から214日で重なる）
MAN_KEEP_DAYS = 60   # 噴射の時刻を控える日数（弱い推力の判定で最近噴射した機を外すため）
STILL_PLANE = 0.15   # まだ近づいているか＝この先の軌道面の差の最小がこれ以下（plane・launch）
STILL_SAME = 0.3     # 同じ軌道面にいるか（phase・keep）


def _iso(t):
    from datetime import datetime, timezone
    return datetime.fromtimestamp(int(t), timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def _el_store(hist, cur, now):
    """控えに cur を足す（いちばん新しい1件は入れ替え・3日おきに残す・30日より古いものは捨てる）"""
    if hist and cur[0] - hist[-1][0] < EL_GAP_DAYS * DAY and len(hist) >= 2 and hist[-1][0] - hist[-2][0] < EL_GAP_DAYS * DAY:
        hist[-1] = cur
    else:
        hist.append(cur)
    cut = now - KEEP_EL_DAYS * DAY
    while len(hist) > 2 and hist[0][0] < cut:
        hist.pop(0)


def events(objs, els, now, new_days=LAUNCH_DAYS):
    """控え els（{番号: {'h': [控え], 'm': [噴射の時刻(秒)]}}）と今の全カタログを比べて噴射を拾う。els は更新する"""
    from datetime import datetime, timezone
    out = []
    seen = set()
    for x in objs:
        k = x['NORAD_CAT_ID']; seen.add(k)
        try:
            cur = el_of(x)
        except Exception:
            continue
        rec = els.get(k)
        if rec is None:
            els[k] = {'h': [cur], 'm': []}
            ld = x.get('LAUNCH_DATE')
            if ld and els.get('_ready') and (now - datetime.fromisoformat(ld).replace(tzinfo=timezone.utc).timestamp()) / DAY <= new_days:
                out.append((x, None, cur, 'launch', None))
            continue
        h = rec['h']
        if cur[0] <= h[-1][0]:
            continue
        m = maneuver(h[-1], cur)
        if m:
            rec['m'].append(cur[0])
            out.append((x, h[-1], cur, 'man', m))
        else:
            old = [e for e in h if cur[0] - e[0] >= SLOW_DAYS * DAY]
            recent = [t for t in rec['m'] if cur[0] - t < (SLOW_DAYS + 7) * DAY]
            if old and not recent and abs(cur[1] - old[-1][1]) >= SLOW_DI:
                out.append((x, old[-1], cur, 'slow', {'di': round(cur[1] - old[-1][1], 4)}))
        _el_store(h, cur, now)
        rec['m'] = [t for t in rec['m'] if now - t < MAN_KEEP_DAYS * DAY]
    for k in [k for k in els if not k.startswith('_') and k not in seen]:
        del els[k]
    els['_ready'] = True   # 1回目（控えが空）の全機を「打上げ直後」と取り違えないため
    return out


def regeom(a_el, b_el, now):
    """今の幾何（控えの最新どうし）：軌道面の差・この先の最小とその日・高度差・傾斜角の差・位相"""
    Wa, ua = at(a_el, now); Wb, ub = at(b_el, now)
    dWa, dua = rates(a_el[3] + RE, a_el[4], a_el[1]); dWb, dub = rates(b_el[3] + RE, b_el[4], b_el[1])
    days = np.arange(0, HORIZON + 1, 1.0)
    p = plane_deg(a_el[1], Wa + dWa * days, b_el[1], Wb + dWb * days)
    k = int(np.argmin(p))
    dphi = (ub - ua + 540) % 360 - 180
    rel = dua - dub
    ph = abs(dphi / rel) if rel * dphi > 0 else None
    return {'planeDeg': round(float(p[0]), 3), 'minDeg': round(float(p[k]), 3), 'inDays': k,
            'dInc': round(b_el[1] - a_el[1], 3), 'dAltKm': round(b_el[3] - a_el[3], 2),
            'phaseDeg': round(float(dphi), 1), 'phaseDays': round(float(ph), 1) if ph is not None and ph < 3650 else None}


def update(objs, now, state):
    """state['els']・state['targets'] を更新し、一覧に出す組を返す。[(a番号, b番号, 中身)]"""
    els = state.setdefault('els', {})
    T = state.setdefault('targets', {})
    by = {x['NORAD_CAT_ID']: x for x in objs}
    ev = events(objs, els, now)
    if ev:
        C = Cat(objs, now)
        for x, pre, cur, kind, m in ev:
            for tid, o in scan(x, pre, cur, C, kind):
                key = x['NORAD_CAT_ID'] + '>' + tid
                r = T.setdefault(key, {'first': _iso(now), 'mans': [], 'modes': []})
                if o['mode'] not in r['modes']:
                    r['modes'].append(o['mode'])
                if m and not any(abs(t['t'] - cur[0]) < 3600 for t in r['mans']):
                    r['mans'].append({'t': cur[0], 'kind': kind, **{k: float(v) for k, v in m.items()}})
                r['last'] = cur[0]
    out = []
    for key in list(T):
        a, b = key.split('>')
        r = T[key]
        if a not in by or b not in by or now - r.get('last', now) > ACTIVE_DAYS * DAY:
            del T[key]; continue
        try:
            g = regeom(el_of(by[a]), el_of(by[b]), now)
        except Exception:
            continue
        still = g['minDeg'] <= STILL_PLANE or (g['planeDeg'] <= STILL_SAME and abs(g['dAltKm']) <= PH_DALT)
        if not still:
            del T[key]; continue
        if r['modes'] == ['keep'] and len(r['mans']) < KEEP_HITS:
            continue   # 同じ軌道面で保つだけの組は、噴射が2回以上で出す
        meet = _iso(now + g['inDays'] * DAY) if g['minDeg'] <= PL_HIT and g['inDays'] > 0 else None
        out.append((a, b, {'modes': r['modes'], 'mans': [_iso(t['t']) for t in r['mans']],
                           'manDetail': [{k: v for k, v in t.items() if k != 't'} for t in r['mans']],
                           'meetAt': meet, **g, '_first': r['first']}))
    return out
