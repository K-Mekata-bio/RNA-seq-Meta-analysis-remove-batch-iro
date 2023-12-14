## バイオインフォマティクスデータ分析スクリプトの詳細要約
- **RNA-seqのメタアナリシスの際に使用するコードです。カウントデータを統合する際に発生するバッチエフェクトをCombatかLimma-voomのremoveBatchEffectを使って除去します。発現変動遺伝子の同定は、limmaパッケージを使用しています。EnhancedVolcanoは簡単にVolcano plotを作成できるので、おすすめです。このコードはversion.1につき今後様々な変更が付け加えられる可能性があります。**
### 注意点
- **このコードはVersion.1です。適宜変更して使ってください:**
　- 気が向けば自動化しやすいコードにします。無理だ？頑張ってください(´･ω･`)
  - TMM法とLog変換による正規化の二種類のファイルを用意しています。お好きな方を使ってください。
- **バッチエフェクトを除去する際には、どのグループをリファレンスにするか決める必要があります:**
  - Set "C" group as the reference level ("H"/"C")　coldata$group <- relevel(coldata$group, ref = "C")
- **DEGの同定に使う閾値を変更してください:**
  - デフォルトは、fdr_threshold <- 0.05　logFC_threshold <- log2(1.5)　です。
- **coef=""の名前は結果に合わせて変更してください:**
  - デフォルトは、results_raw <- topTable(fit_raw, coef="groupH", number=Inf, sort.by="p")です。
- **このページはChatGPT-4に要約させたものを編集して作っています**

### ファイル
- **main-TMM.r:**
  - メインファイル　TMM method
- **main-log.r:**
  - メインファイル Log method
- **coldata.csv & countdata.csv:**
  - サンプルファイル

### インストールとライブラリの読み込み
- **CRANからのパッケージインストール:**
  - `ggplot2`, `gridExtra`, `FactoMineR`, `reshape2`, `shiny`, `ggpubr` をインストール
- **Bioconductorからのパッケージインストール:**
  - `EnhancedVolcano`, `org.Hs.eg.db`, `sva`, `edgeR`, `limma` をインストール
- **ライブラリの読み込み:**
  - 先にインストールされたすべてのパッケージを読み込み

### データの読み込みと処理
- **遺伝子発現データ(`countdata.csv`)とサンプルメタデータ(`coldata.csv`)の読み込み:**
  - ヘッダーあり、行名を最初の列に設定してCSVファイルから読み込み
- **データの前処理:**
  - 遺伝子発現データを行列に変換
  - DGEリストオブジェクトの作成とフィルタリング
  - Log変換による正規化とCSVへの保存

### バッチ効果の除去
- **バッチ効果除去のための `ComBat` 関数:**
  - Log変換されたデータに対してバッチ効果を除去
  - Log-countをCSVに保存

### 差異発現遺伝子（DEG）分析
- **バッチ効果除去前後のDEG分析:**
  - FDR閾値とfold change閾値の設定
  - 線形モデルの適用、仮説検定、結果の抽出
  - 多重検定調整後のP値の計算とDEGの選定
  - 遺伝子シンボルへのマッピングとCSVへの保存

### 可視化
- **Volcano plot（`EnhancedVolcano`を使用）:**
  - LogfoldchangeとFDRを軸にしたプロット作成
  - PDFとPNG形式での保存
- **BOX plot（`ggpubr`を使用）:**
  - サンプルごとの正規化カウントをプロット
  - PDF形式での保存
- **PCA分析（`prcomp`を使用）:**
  - バッチ効果除去前のデータにPCA分析を適用

### データの保存
- 分析結果や可視化結果をCSV、PDF、PNGファイルとして保存
