# 再突入予測（Space-Track の TIP メッセージ）を取って data/tip.json にする。OP's LAB Maps の「再突入予測」レイヤー用。
# ・物体ごとに最新の1通だけを残す（TIP は再突入の4日前から段階的に出直される）。
# ・地上軌跡を描くため、その物体の最新TLEも同じファイルに入れる（大半は Starlink で tle_recent.json に居ない）。
# ・Space-Track の目安は TIP の問い合わせが1時間に1回＝毎時のワークフローで1回だけ呼ぶ。
import json
import os
import sys
from datetime import datetime, timedelta, timezone

import requests

BASE = 'https://www.space-track.org'
OUT = 'data/tip.json'
# 予測時刻がこれより前の物は外す（再突入後の確定報告をアプリで24時間見せるので、余裕を持って2日）
KEEP_PAST_DAYS = 2


def main():
    s = requests.Session()
    r = s.post(BASE + '/ajaxauth/login',
               data={'identity': os.environ['SPACETRACK_USERNAME'],
                     'password': os.environ['SPACETRACK_PASSWORD']}, timeout=60)
    r.raise_for_status()

    # 直近10日に出た通知（予測が4日先まで出るので、それより長めに取って最新の1通を選ぶ）
    r = s.get(BASE + '/basicspacedata/query/class/tip/INSERT_EPOCH/%3Enow-10'
                     '/orderby/INSERT_EPOCH%20desc/format/json', timeout=120)
    r.raise_for_status()
    msgs = r.json()

    cut = datetime.now(timezone.utc) - timedelta(days=KEEP_PAST_DAYS)
    latest = {}
    for m in msgs:   # 新しい順に並んでいるので、最初に出た1通がその物体の最新
        n = m.get('NORAD_CAT_ID')
        if not n:
            continue
        if n in latest:
            if latest[n]:
                latest[n]['nMsg'] += 1
            continue
        try:
            dec = datetime.strptime(m['DECAY_EPOCH'], '%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc)
        except (KeyError, ValueError):
            continue
        if dec < cut:
            latest[n] = None   # 古い物は印だけ付けて外す（後の通知で数え直さないように）
            continue
        latest[n] = {
            'id': n,
            'msgEpoch': m.get('MSG_EPOCH', ''),
            'decay': m.get('DECAY_EPOCH', ''),
            'window': int(float(m.get('WINDOW') or 0)),          # ±分
            'nextReport': int(float(m.get('NEXT_REPORT') or 0)),  # 時間・0＝再突入後の確定報告
            'lat': float(m.get('LAT') or 0),
            'lon': float(m.get('LON') or 0),
            'incl': float(m.get('INCL') or 0),
            'dir': m.get('DIRECTION', ''),
            'rev': m.get('REV', ''),
            'high': m.get('HIGH_INTEREST', 'N') == 'Y',
            'nMsg': 1,
        }
    items = [v for v in latest.values() if v]

    # 物体の名前・種別・最新TLE（再突入済みの物は gp に無いことがあるので gp_history の最後で補う）
    ids = [v['id'] for v in items]
    gp = {}
    if ids:
        r = s.get(BASE + '/basicspacedata/query/class/gp/NORAD_CAT_ID/' + ','.join(ids) + '/format/json',
                  timeout=120)
        r.raise_for_status()
        for g in r.json():
            gp[g['NORAD_CAT_ID']] = g
        miss = [i for i in ids if i not in gp]
        if miss:
            r = s.get(BASE + '/basicspacedata/query/class/gp_history/NORAD_CAT_ID/' + ','.join(miss) +
                      '/EPOCH/%3Enow-30/orderby/EPOCH%20desc/format/json', timeout=180)
            if r.ok:
                for g in r.json():
                    gp.setdefault(g['NORAD_CAT_ID'], g)
    for v in items:
        g = gp.get(v['id'], {})
        v['name'] = g.get('OBJECT_NAME', '')
        v['type'] = g.get('OBJECT_TYPE', '')
        v['intl'] = g.get('OBJECT_ID', '')
        v['tleEpoch'] = g.get('EPOCH', '')
        v['l1'] = g.get('TLE_LINE1', '')
        v['l2'] = g.get('TLE_LINE2', '')

    items.sort(key=lambda v: v['decay'])
    os.makedirs('data', exist_ok=True)
    with open(OUT, 'w', encoding='utf-8') as f:
        json.dump({'updated_at': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
                   'source': 'Space-Track.org TIP', 'count': len(items), 'data': items},
                  f, separators=(',', ':'), ensure_ascii=False)
    print(f'✓ 再突入予測 {len(items)} 件を保存しました（通知 {len(msgs)} 通・TLEなし '
          f'{sum(1 for v in items if not v["l1"])} 件）')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'✗ 再突入予測の取得エラー: {e}')
        sys.exit(1)
