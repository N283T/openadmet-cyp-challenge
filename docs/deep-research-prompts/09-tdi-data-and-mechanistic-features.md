# TDI data and mechanistic features

CYP TDI/MBI予測に利用できる外部データと、計算可能なmechanistic featureをDeep Researchしてください。OpenADMET CYP ChallengeのCYP2D6/CYP3A4 boolean TDIとendpoint互換性があるかを最優先し、回答は日本語でお願いします。

PubChem、ChEMBL、BindingDB、FDA approval packages、文献補足資料、特許、規制資料、公開リポジトリ、代謝・site-of-metabolism・reactive metaboliteデータを探索してください。各候補について、可能な範囲で以下を確認してください。

- 化合物数、対象CYP、陽性・陰性数。
- IC50 shift、k_inact/K_I、k_obs、定性的TDI、clinical DDI、通常のCYP inhibitionのどれを収録するか。
- assay条件、構造・識別子、取得方法、形式、ライセンス、現在の入手可否。
- OpenADMET/ChEMBL/PubChemとの重複、系列リーク、ラベル互換性、陰性の信頼性。
- 直接教師、補助タスク、事前学習、機構検証用のどれに適するか。

さらに、公開実装で計算可能なTDI structural alerts、bioactivation motifs、electrophile/reactivity、CYP site-of-metabolism、metabolite enumeration、heme-access geometryなどを調べてください。商用専用機能と再現可能なオープン手法を区別してください。

最後に、今すぐ取得・実装できる候補を優先順位付けし、最小実験、必要な重複除外、成功・失敗判定を提案してください。推測で欠損情報を埋めず、重要な主張には直接URLを付けてください。
