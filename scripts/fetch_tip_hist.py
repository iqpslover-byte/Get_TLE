# 過去の再突入（Space-Track の TIP の確定報告）を取って data/tip_hist.json にする。OP's LAB Maps の再突入予測の表用。
# ・確定報告（NEXT_REPORT=0＝再突入後の報告）のうち、予測時刻が直近 KEEP_DAYS 日の物。物体ごとに最後の1通。
# ・地上軌跡を描くため、その物体の最後のTLEも同じファイルに入れる（gp は再突入済みの物も最後の1組を持っている）。
# ・確定報告はもう変わらないので1日1回で足りる＝毎時のワークフローから呼ばれても、前回から20時間たつまでは何もしない。
import json
import os
import sys
from datetime import datetime, timedelta, timezone

import requests

BASE = 'https://www.space-track.org'
OUT = 'data/tip_hist.json'
KEEP_DAYS = 365
MIN_INTERVAL_H = 20
CHUNK = 400   # gp の問い合わせ1回に入れる番号の数（URL が長くなりすぎないように）


def fresh():
    try:
        with open(OUT, encoding='utf-8') as f:
            t = datetime.strptime(json.load(f)['updated_at'], '%Y-%m-%dT%H:%M:%SZ').replace(tzinfo=timezone.utc)
        return datetime.now(timezone.utc) - t < timedelta(hours=MIN_INTERVAL_H)
    except Exception:
        return False


def main():
    if fresh() and '--force' not in sys.argv:
        print('✓ 過去の再突入は前回から20時間たっていないので取り直しません')
        return
    s = requests.Session()
    r = s.post(BASE + '/ajaxauth/login',
               data={'identity': os.environ['SPACETRACK_USERNAME'],
                     'password': os.environ['SPACETRACK_PASSWORD']}, timeout=60)
    r.raise_for_status()

    r = s.get(BASE + '/basicspacedata/query/class/tip/NEXT_REPORT/0/DECAY_EPOCH/%3Enow-' + str(KEEP_DAYS) +
                     '/orderby/INSERT_EPOCH%20desc/format/json', timeout=180)
    r.raise_for_status()
    msgs = r.json()

    latest = {}
    for m in msgs:   # 新しい順＝最初に出た1通がその物体の最後の確定報告
        n = m.get('NORAD_CAT_ID')
        if not n or n in latest or not m.get('DECAY_EPOCH'):
            continue
        latest[n] = {
            'id': n,
            'msgEpoch': m.get('MSG_EPOCH', ''),
            'decay': m.get('DECAY_EPOCH', ''),
            'window': int(float(m.get('WINDOW') or 0)),
            'nextReport': 0,
            'lat': float(m.get('LAT') or 0),
            'lon': float(m.get('LON') or 0),
            'incl': float(m.get('INCL') or 0),
            'dir': m.get('DIRECTION', ''),
            'rev': m.get('REV', ''),
            'high': m.get('HIGH_INTEREST', 'N') == 'Y',
        }
    items = list(latest.values())

    ids = [v['id'] for v in items]
    gp = {}
    for k in range(0, len(ids), CHUNK):
        part = ids[k:k + CHUNK]
        r = s.get(BASE + '/basicspacedata/query/class/gp/NORAD_CAT_ID/' + ','.join(part) + '/format/json',
                  timeout=180)
        r.raise_for_status()
        for g in r.json():
            gp[g['NORAD_CAT_ID']] = g
    for v in items:
        g = gp.get(v['id'], {})
        v['name'] = g.get('OBJECT_NAME', '')
        v['type'] = g.get('OBJECT_TYPE', '')
        v['intl'] = g.get('OBJECT_ID', '')
        v['tleEpoch'] = g.get('EPOCH', '')
        v['l1'] = g.get('TLE_LINE1', '')
        v['l2'] = g.get('TLE_LINE2', '')

    items.sort(key=lambda v: v['decay'], reverse=True)
    os.makedirs('data', exist_ok=True)
    with open(OUT, 'w', encoding='utf-8') as f:
        json.dump({'updated_at': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
                   'source': 'Space-Track.org TIP', 'days': KEEP_DAYS, 'count': len(items), 'data': items},
                  f, separators=(',', ':'), ensure_ascii=False)
    print(f'✓ 過去の再突入 {len(items)} 件を保存しました（確定報告 {len(msgs)} 通・TLEなし '
          f'{sum(1 for v in items if not v["l1"])} 件）')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(f'✗ 過去の再突入の取得エラー: {e}')
        sys.exit(1)
