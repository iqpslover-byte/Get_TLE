# Get_TLE

Space-Track から最新の衛星 TLE（Two-Line Element）データを定期取得するリポジトリです。

## 仕組み

GitHub Actions が毎時間自動実行され、Space-Track API から最新100件の衛星 TLE を取得して `data/tle.json` に保存します。

## ファイル構成

```
.github/
  └── workflows/
      └── fetch-tle.yml       # GitHub Actions ワークフロー
data/
  └── tle.json                # 取得した TLE データ（JSON形式）
README.md
```

## GitHub Secrets の設定

GitHub リポジトリの Settings → Secrets に以下を登録してください：

- `SPACETRACK_USERNAME` : Space-Track のユーザー名（メールアドレス）
- `SPACETRACK_PASSWORD` : Space-Track のパスワード

## 使用方法

### 自動実行
毎時間 UTC 0分に自動実行されます。

### 手動実行
GitHub リポジトリの Actions タブから「Fetch Space-Track TLE Data」を選択して「Run workflow」をクリックします。

### OP's LAB Maps での使用

```javascript
const response = await fetch('https://raw.githubusercontent.com/iqpslover-byte/Get_TLE/main/data/tle.json');
const tleData = await response.json();
const satellites = tleData.data;
```

## ライセンス

Space-Track API の利用規約に準拠します。
