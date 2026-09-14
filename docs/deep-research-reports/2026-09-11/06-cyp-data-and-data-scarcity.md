# CYP研究で実際に使えるデータと、データ不足を補うための実践的な戦略

## まず結論：最も価値が高い「追加情報」は何か

ご提示の「ヒトCYP1A2、2C9、2D6、3A4のdirect-inhibition pIC50が各約1,300–2,300件」「単一濃度スクリーニングあり」「2D6・3A4のTDI」という構成は、2026年8月に公開されたOpenADMET CYP Inhibition Blind Challengeのデータ構成と非常によく一致します。同チャレンジでは、direct inhibition、単一濃度一次スクリーニング、±NADPH preincubationによるIC50-shift TDI、Emaxが相互に関連付けられており、**これらを独立したラベル群としてではなく、一つの多忠実度（multi-fidelity）実験系として扱うことが、外部データを数千件増やすこと以上に重要**です。citeturn14search3turn18view0

調査の結果、ChEMBL 36と古典的PubChem qHTS以外で、特に確認価値が高いのは次の系統です。

| 候補 | 現在の入手性 | 主に得られるもの | direct pIC50への使い方 | 独立情報としての期待 |
|---|---|---|---|---|
| **EPA/Tox21 10K CYP P450-Glo系** | 公開。CompTox/invitroDB系 | CYP1A2/2C9/2C19/2D6/3A4阻害、濃度反応/AC50・hit call | **別アッセイ補助タスク**、重複化合物でassay calibration | **非常に高い**。報告上9,667化合物規模。ただしChEMBL/PubChemとの重複確認必須 |
| **EPA ToxCast NVS_ADME CYP系** | 公開 | Novascreen CYP enzyme inhibition AC50、inactiveを含む | assay-specific auxiliary head、環境化学領域への外挿 | 中～高。Veith qHTSとは別測定系 |
| **OpenADMET 2026年3月 reaction phenotyping** | 公開 | CYP3A4/CYP2J2による親化合物depletion、well-level peak area等、1,221化合物 | CYP3A4の「結合/阻害」ではなく**代謝されやすさ**の補助タスク | 高い。ただしChallenge化合物との構造重複は要計算 |
| **OpenADMET Challengeのwell-level/posterior samples** | **Challenge終了後に要申請予定** | 生シグナル、pIC50/slope posterior | label uncertainty、censoring、heteroscedastic regression | 非常に高い。同じ測定から情報量を増やす |
| **NCATS dietary supplement / natural product qHTS** | 公開、PubChem | CYP1A2等の天然物・抽出物阻害 | 化学的に定義された化合物のみ補助利用 | 中。化合物と「substance/extract」を分ける必要 |
| **BindingDB独自文献・特許抽出** | 公開 | IC50/Ki等 | 未ChEMBL文献・特許の発見 | 中。ただし**ChEMBL/PubChemとの循環重複が大きい** |
| **FDA clinical index substrate/inhibitor情報** | 公開 | clinical substrate/inhibitor/inducer分類 | direct regressionではなく外部妥当性・機構検証 | 高品質だが小規模 |
| **SABIO-RK** | 非商用は無償、商用ライセンス | Km、Ki、kcat等＋pH/温度/変異体/反応条件 | kinetics/mechanistic auxiliary | 小規模だが条件情報が濃い |
| **ATTNSOM等のSoMデータ** | 公開 | atom-level CYP site-of-metabolism | encoder事前学習・代謝機構補助 | 中～高 |
| **Certara DIDB** | 商用 | curated drug-interaction evidence | TDI/DDIの文献発見・評価系 | 特にTDIでは高い可能性 |
| **特許・regulatory文書から新規抽出** | 文書自体は公開多数 | IC50、Ki、kinact/KI、陰性、条件 | 直接ラベルになり得る | **未回収量は大きい可能性、工数も大きい** |

特に重要な新規候補は**Tox21 10KのP450阻害系**です。EPAのassay descriptionには `TOX21_p450_CYP1A2_Antagonist`、`...CYP2C9...`、`...CYP2D6...`、`...CYP3A4...` が存在し、2025年の報告はTox21 10Kライブラリ9,667化合物を5種のP450-Glo assayでスクリーニングしたとしています。これはVeith 2009の17K qHTSとは別系統です。citeturn41search1turn41search10turn39search5  
ただし、**「9,667件をそのまま9,667個の新しいdirect IC50として足せる」という意味ではありません**。EPA pipelineでのcurve fitting、hit call、QC flags、化合物標準化後のユニーク数、ChEMBL 36への既収載分を照合して初めて増分が分かります。この確認を最優先にする価値があります。

また、OpenADMET側ではChallenge終了後に**well-level assay signalsと各化合物のpIC50/slope posterior samplesをrequest可能にする予定**と明記されています。現在は2026年9月11日でChallenge終了予定は11月3日なので、これは「今すぐ」ではなく**近い将来の最重要取得候補**です。citeturn18view0turn14search9


## 直接阻害とTDIで実際に増やせるデータ

### OpenADMET自身に、まだ使い切れていない情報がかなりある

Challenge dataset cardでは、direct inhibition用ファイルは4,905化合物を行として持ち、4 CYPのpIC50、95% confidence interval、standard deviation等を収録しています。一方、TDI関連ファイルは6,145化合物規模、single-concentrationデータは**17,504 measurement / 4,376 compounds × 4 CYP**、Emaxデータは6,146化合物規模です。したがって「direct pIC50が各1,300–2,300しかない」という見方より、実際には**同じ化学空間に対して濃淡の異なる観測が数万件存在する**と見る方が適切です。citeturn14search3turn18view0

single-concentration screeningは単なる補助的なbinary hit flagではありません。公開説明によれば一次スクリーニングはおおむね50 µMで行われ、その結果からdose-response follow-up対象が選択されています。follow-upでは12-point concentration responseが実施され、さらに1A2/2C9/3A4のtop hitsから近傍化合物が追加されました。つまり現在のpIC50ラベル集合は**ランダムサンプルではなく、一次screeningによるselectionを条件とした集合**です。citeturn17view0turn18view0

これはかなり重大です。単純にpIC50既測定化合物だけで学習すると、

\[
P(y_{\mathrm{pIC50}}\mid x,\; \text{selected for DRC})
\]

を学ぶことになり、本当に欲しい

\[
P(y_{\mathrm{pIC50}}\mid x)
\]

との間にverification/selection biasが入ります。したがってsingle-concentrationデータは「教師データをちょっと増やす」以上に、**なぜその化合物にIC50測定が存在するのかというselection mechanismをモデル化するためのデータ**として使う価値があります。

さらに、この一次screeningの全体分布にはdose-responseに進まなかった大量の弱活性・不活性化合物が含まれるため、negative chemical spaceを埋める役割があります。モデルとしては、

> structure → single-concentration response  
> structure → latent potency  
> latent potency + assay condition → DRC pIC50/Emax

という階層構造にするのが、single-concentration responseを無理に「推定pIC50」へ変換するより安全です。

OpenADMETのassay条件もモデルに入れるべきです。1A2/2C9/3A4はSupersome＋Vivid fluorogenic probe系で、詳細な公開記事では3A4=DBOMF、2C9=BOMF、1A2=EOMCC、一方2D6はfluorescence系がうまく機能せず、dextromethorphan→dextrorphanをEcho-MS/QTOFで読み出す系へ変更されています。つまり**2D6と他3酵素はmeasurement modality自体が違います**。citeturn18view0

これはcross-CYP multitask modelでは重要です。CYP identityだけでなく、

`assay_platform = fluorescence / MS`  
`probe_substrate`  
`preincubation_NADPH`  
`preincubation_time`  
`readout = product formation / probe turnover`

を条件として与えないと、shared encoderが生物学的cross-CYP signalとplatform-specific signalを混同します。

### OpenADMET 2026年3月のreaction phenotypingは別タスクとして有望

Challenge以前の2026年3月OpenADMET releaseには、**1,221化合物のCYP3A4/CYP2J2 reaction phenotyping**があり、well-level peak area、parent depletion ratioなどが公開されています。同時期のCYP3A4 inhibition datasetは1,343化合物規模ですが、こちらはactive CYPと30分preincubationしているため、純粋なreversible/direct inhibitionではなく、directとmetabolism-dependent componentが混ざり得ます。citeturn15search0

したがってこの旧CYP3A4 inhibition値を現在の「direct pIC50」にそのまま追加するのは勧めません。

一方、reaction phenotyping部分はかなり面白いです。親化合物depletionは

> 「この分子がCYP3A4に結合して活性を阻害するか」

ではなく、

> 「CYP3A4によって実際に消費されるか」

に近い情報です。

これはdirect inhibitionと相関し得るが同じではなく、**substrate/inhibitor dualityを分離する補助task**として使えます。同一化合物がChallenge側にも存在すれば、むしろ重複は利点です。たとえば、

\[
\text{inhibition potency},\quad
\text{parent depletion},\quad
\text{TDI shift}
\]

を同一分子について持つと、「強いreversible binder」「rapid substrate」「metabolism-activated inhibitor」の表現を分離しやすくなります。

### Veith 2009 qHTSは「数」より、ラベルの意味を再利用する

古典的なPubChem AID 1851は17,143 samplesを各CYPに対して7–15濃度で測定したqHTSです。最高濃度は57 µM、CRCをpotency・efficacy・fit qualityから分類し、no responseも明示的にinactiveとして持っています。citeturn32search1

後続のcurationでは無機物・混合物等を除いて16,359化合物、conclusiveなCYP–compound pairは63,391件となり、isoform別にconclusive dataは1A2 12,634、2C9 12,264、2D6 13,276、3A4 12,383でした。citeturn32search0

ただし、このデータには二つの重要な制約があります。

第一に、Veith assayではpro-luciferin turnover低下を読んでおり、著者自身が**substrateとして酵素を消費・競合する化合物も「inhibitory signal」を生み得る**と注意しています。第二にcompound preincubationを行っていないため、TDI/mechanism-based inhibitionには低感度です。citeturn32search1

したがって、AID 1851のAC50をOpenADMET direct pIC50と同一の回帰目的変数へpoolするより、

\[
y_{\text{OpenADMET-direct}}
\neq
y_{\text{Veith-qHTS}}
\]

として**shared molecular encoder＋assay-specific heads**にする方が理にかなっています。

さらに、重複化合物を捨てるのではなく、同一化合物に

\[
pIC50_{\mathrm{OpenADMET}}-pAC50_{\mathrm{Veith}}
\]

が存在するpaired setを作ると、**probe/platform-dependent assay shiftそのものを学習・定量化**できます。これは新しい化合物を数千追加するより、外部データをどの程度信じてよいかを判断する上で重要です。

### Tox21 10K P450は最優先でローカル抽出を試すべき

EPA系データで見落としやすいのは、古いToxCast NVS_ADMEだけではなく、`TOX21_p450_CYP1A2_Antagonist`、`CYP2C9`、`CYP2D6`、`CYP3A4` などのTox21 P450 endpointsが存在することです。EPAのassay-description資料にもこれらのendpointが記載されています。citeturn41search1turn39search5

2025年の報告では、**Tox21 10K libraryの9,667 compoundsをCYP1A2、2C9、2C19、2D6、3A4のP450-Glo assayでscreeningした**とされています。citeturn41search10

ここで実際に行うべき抽出は「active compoundsだけを取る」ではなく、EPA invitrodb/CompToxから、

`DTXSID / chemical ID`  
`AEID / endpoint`  
`hitcall`  
`AC50`  
`top / efficacy`  
`curve-fitting/QC flags`  
`tested but inactive`

まで取ることです。CompToxはToxCast/Tox21 bioactivityとdose-response関連情報をダウンロード可能な形で提供しています。citeturn55search35

**不確実性**として、私は今回の公開ページから9,667全化合物について「各isoformで何件が有効AC50を持つか」「ChEMBL 36に何件入っているか」までは確認できませんでした。ここは推測で埋めるべきではありません。ただし規模から見て、ChEMBL 36 local databaseに対して最初にoverlap analysisする価値は非常に高いです。

### ToxCast NVS_ADMEは別の独立測定層

ToxCastには `NVS_ADME_hCYP1A2`、`hCYP2C9`、`hCYP2D6`、`hCYP3A4` 等のNovascreen系CYP enzyme inhibition endpointsがあります。公開例ではAC50がµM単位で格納され、ToxCast初期campaignは976 chemicalsを331 enzymatic/receptor signaling assaysにprofilingした研究として報告されています。citeturn39search8turn39search18

これは化学空間が医薬品中心ではなく環境化学物質にも広がるため、direct pIC50の精密な追加教師というより、

- CYP binding/inhibitionに関係する一般的chemical representationの学習
- approved/drug-like space外へのOOD robustness確認
- Veith/OpenADMETと重複した化合物のcross-platform validation

に向いています。

Tox21 P450-GloとToxCast NVS_ADMEは名前が似ていますが、**別assay familyとして維持すべき**です。

### PubChemには天然物・サプリメント系の比較的新しいCYPデータもある

NCATSのPubChem AID 1919976は2024年depositの「qHTS for inhibitors of CYP1A2 against the NCATS DSHEA and TCM libraries」です。ここでは964 tested substancesに対して、PubChem上で構造化されたtested compoundsは209、active compounds 40、inactive compounds 114と表示されています。citeturn32search4

これは「measurement rows」と「unique molecular structures」を区別すべき好例です。多数がbotanical extracts、mixtures等なので、964を分子モデルの教師件数として数えてはいけません。研究自体ではCYP1A2、2D6、3A4等を含むdrug-interaction panelが用いられています。citeturn33search0

用途としては、構造が明確な209 compoundsの範囲で補助利用するか、extract-level dataを「分子教師」ではなく天然物相互作用の外部評価に残すのが安全です。


## 基質性・代謝部位・速度・構造・変異体をどう組み合わせるか

direct inhibitionとの相関が完全でないデータほど、**別taskとして使った時に価値があります**。

### 基質性

Therapeutics Data Commonsには、Veith inhibition datasetsとは別に

`CYP2C9_Substrate_CarbonMangels`  
`CYP2D6_Substrate_CarbonMangels`  
`CYP3A4_Substrate_CarbonMangels`

が公開されています。citeturn32search5

ただし、今回の調査ではCarbon–Mangels各datasetについて、元文献から現在のTDC tableまでの「全compound単位のprovenance」とexact unique countを十分確認できませんでした。したがって、**TDCのCSVを独立測定データとして無条件に加えることは勧めません**。まずInChIKey/structureでChEMBL、PubChem、OpenADMETと重ね、reference columnや原論文を遡るべきです。

それでもタスク自体は有用です。inhibitorとsubstrateは、

- CYP active siteへ到達する
- isoform-specificなchemical recognitionを受ける

という部分を共有します。そのためsubstrate/non-substrate classificationをauxiliary headに置くと、binding-compatibleなchemical representationを学ばせられます。一方、substrate=non-inhibitorでもinhibitor=non-substrateでもないため、direct pIC50へbinary labelを変換してはいけません。Veith assay自身がsubstrateとinhibitorのsignal混同を認めていることも、この分離の重要性を示します。citeturn32search1

### site-of-metabolismとmetabolite

2026年のATTNSOMはatom-level CYP site-of-metabolism predictionについて**二つのbenchmark datasets**を用い、isoformを独立に扱うよりcross-isoform attentionを導入することで、複数CYP間の代謝パターンを利用しています。コードとdatasetも公開すると明記されています。citeturn56academia6

これはdirect pIC50の教師にはなりませんが、分子encoderのpretraining taskとして非常に自然です。

例えば、

\[
L =
L_{\mathrm{direct\ pIC50}}
+\lambda_1 L_{\mathrm{substrate}}
+\lambda_2 L_{\mathrm{SoM}}
+\lambda_3 L_{\mathrm{depletion}}
\]

とすれば、分子全体のpotencyと、原子単位のmetabolic liabilityを同時に学ばせられます。

CYP metabolismでは「原子のintrinsic reactivity」と「その原子をheme近傍に置けるbinding geometry」の両方が重要というmechanistic frameworkがあり、CYP3A4/2D6 SoM予測でもQSAR reactivityとinduced-fit dockingを組み合わせる研究があります。citeturn56academia8

したがって、SoM予測値を単に追加descriptorとして入れるより、

> molecular encoder  
> → atom reactivity head  
> → isoform-specific accessibility/SoM head  
> → molecule-level inhibition head

という構造にする方が機構的には整合的です。

ただしこれは**CYP SoMでは実証があるが、direct-inhibition pIC50改善については今回確認できた範囲では直接実証ではなく類推**です。

### kineticsはSABIO-RKが「件数は少なくても情報密度が高い」

SABIO-RKは文献からmanual curationされたreaction kinetics databaseで、substrates/products、inhibitors/activators/cofactors、wild type/mutant、Km/Ki等のkinetic parameters、rate equations、organism、tissue、cellular location、pH、temperature、buffer条件まで保持します。REST APIとtable/SBML等のexportがあり、非商用利用はfree、commercial usersはlicenseが必要です。citeturn54search2

CYP研究では大量学習データ源というより、

- 同じcompoundでIC50とKiの関係を確認する
- substrate concentration依存性を理解する
- variantによるKm/kcat/Ki変化を見る
- assay-condition featuresを定義する
- mechanism validation setを作る

用途が適しています。

特に外部文献から `IC50 = 0.3 µM` だけ拾うより、substrate concentrationとKmがあればCheng–Prusoff的な関係や競合性の違いを考えられるため、**数十～数百件でも条件付きモデルの校正には価値があります**。

### 構造は「dock scoreをそのまま特徴量」よりmechanistic sanity checkに向く

RCSB PDBにはhuman CYP ligand-bound structuresが公開されており、CYP3A4 inhibitor-bound structureなどを直接取得できます。例えば4D78はCYP3A4–inhibitor complexです。citeturn41search5

ただしCYP3A4は大きく柔軟なactive siteとsubstrate-dependent behaviorを持つため、一つのcrystal structureへのstatic docking scoreをpIC50教師に足すだけでは危険です。むしろ、

- 明らかなheme-coordinating Nを持つtype-II binder
- active-site volumeとの整合性
- CYP2D6での塩基性中心
- metabolizable atomのheme accessibility
- known ligand poseと予測attributionが整合するか

といった**機構検証・representation regularization**の方が使いやすいです。

### 変異体データはdirect potencyの増量ではなくrobustness testに使う

CYP2D6やCYP2C9には機能の異なる多数のalleleがあり、PharmVar/CPIC/PharmGKB系統でvariant-to-function knowledgeが整理されています。CYP2C9*2/*3等の機能差は臨床薬理でも利用されています。citeturn54search4

ここでの重要点は、variant activityをwild-type direct inhibition pIC50に混ぜないことです。

将来的には

\[
f(\text{compound},\text{CYP sequence/variant},\text{assay})
\]

というproteochemometric modelに拡張できます。実際、5 CYPをprotein descriptorsとcompound descriptorsで同時に扱った古典的proteochemometric studyでは、63,391 CYP–compound interactionsを単一モデルに統合し、isoform間transferが成立することを示しています。citeturn32search0

同研究でleave-one-isoform-out predictionは平均AUC約0.79まで成立した一方、CYP2D6は他isoformとactivity profileが比較的異なることも示されました。したがって「全部のCYPを完全共有する」のではなく、**shared trunk + isoform-specific adapters/heads**が自然です。citeturn32search0


## 来歴と重複：ChEMBL・PubChem・BindingDBを単純結合してはいけない

今回の調査で、最も実務上重要な見落としの一つは**データベース間の循環的な再配布**です。

BindingDBは独自のjournal literature extractionとpatent extractionに加え、selected PubChem confirmatory BioAssaysとChEMBL entriesを取り込んでいます。さらにBindingDBとChEMBLのデータはPubChem BioAssayにもdepositされます。citeturn55search4turn55search18

したがって、

> ChEMBL + PubChem + BindingDB

を3 independent sourcesとして加算するのは誤りです。

さらにChEMBL自身の公式FAQでは、**BindingDB由来のpatent dataがChEMBL内に `SRC_ID = 37` として含まれる**ことが明示されています。citeturn55search8

つまりChEMBL 36 localをすでにお持ちなら、BindingDBをdownloadする前にまず、

> ChEMBL 36のBindingDB-derived records  
> BindingDB native literature records  
> BindingDB patent records  
> PubChem-origin records

をsource-levelで分解すべきです。

推奨するprovenance unitは「compound」ではなく**measurement event**です。

概念的には以下のkeyを保持します。

```text
structure_parent_key
original_structure_key
source_database
source_record_id
source_document / DOI / patent
depositor
PubChem_AID
assay_id
target_isoform
species
enzyme_system
probe_substrate
substrate_concentration
preincubation_time
NADPH_presence
readout_platform
endpoint_type
relation (=, <, >)
value
unit
activity_comment
QC_flag
```

重複判定も、

```text
same structure + same numeric value
```

ではなく、

```text
same original experiment / same AID / same paper-table-row
```

を基準にする必要があります。

同じ化合物・同じCYPでも**異なるassayなら重複ではありません**。むしろpaired assay dataとして価値があります。

特に以下の三種類を区別すると整理しやすくなります。

| 状況 | 扱い |
|---|---|
| 同一論文・同一assayをChEMBLとBindingDBが再配布 | 一つのmeasurementとしてdeduplicate |
| PubChem qHTSをChEMBLが収載 | 一つのmeasurement。AID lineageを保持 |
| 同一compoundをVeith、Tox21、OpenADMETが独立測定 | **消さない**。assay-condition違いのrepeated measurementとして保持 |

この方針にすると、「追加情報量」は単純なunique structure数だけでなく、

\[
N_{\text{new structures}},
\quad
N_{\text{new compound-assay pairs}},
\quad
N_{\text{paired cross-assay compounds}}
\]

の三つで評価できます。

後者二つは非常に重要です。たとえば新規構造がゼロでも、500化合物がOpenADMETとTox21双方で測られていれば、assay shiftを学習する500 pairsが増えたことになります。

### pIC50・AC50・Ki・single concentrationを一つの数値列にしない

Veith qHTSのAC50はcurve behaviorとefficacyを含むHTS endpointで、OpenADMETのdirect pIC50とはassay systemが異なります。Veith側ではsubstrateもinhibition-like signalを作り得ます。citeturn32search1

さらにOpenADMET自身でもfluorescence probeにはautofluorescence/interferenceやsubstrate dependenceがあり、2D6はそのためMS-based assayへ変更されています。citeturn18view0

したがって、

```text
pIC50_openadmet
pAC50_veith
pAC50_tox21
pIC50_literature_HLM
```

を一列にして平均するのは、training setの見かけのNは増えても意味を薄める可能性があります。

より良いモデルは、

\[
z_i = \text{latent CYP interaction propensity}
\]

を共有し、

\[
y_{ij}=g_j(z_i,\text{assay metadata})+\epsilon_{ij}
\]

として各assayのobservation modelを分離するものです。

### inactiveと打ち切り値を捨てるのは特にもったいない

Veith qHTSは最高57 µMまでno responseのcompoundをexplicit inactiveとして保持します。citeturn32search1  
OpenADMETでも低potency領域、特にpXC50<4はlowest tested dose以下で不確実性が大きいと説明されています。citeturn17view0

これを、

```text
pIC50 = 4.0
```

と固定値化すると人工的なpile-upが生じます。

本来は例えば

\[
pIC50 < 4
\]

という**left-censored observation**です。

回帰lossも、

- exact value → Gaussian / Student-t likelihood
- `<4` → \(P(Y<4)\)
- `>x` → \(P(Y>x)\)
- interval → \(P(a<Y<b)\)

と分ける方が情報を失いません。

95% CIやstandard deviationまであるOpenADMETデータなら、measurementごとにnoise scaleを変えるheteroscedastic regressionも可能です。Challenge終了後にposterior samplesが得られれば、point estimate自体を教師とせずposterior distributionを教師にできます。citeturn14search3turn18view0


## データを増やさず不足を補う方法と、少数追加測定の設計

ここでは「**実測情報そのものを増やす方法**」と「**既存の情報を効率よく使う方法**」を分けます。

### 既存情報の利用効率を上げる：最優先はassay-aware multi-task

最も実証に近いのはcross-CYP multitaskです。前述のproteochemometric studyでは5 CYPを同時学習することでcross-isoform informationを利用できました。citeturn32search0  
また2026年のATTNSOMもsite-of-metabolismでcross-isoform modelingの有効性を示しています。citeturn56academia6

今回なら単なる「4 CYP pIC50 multitask」ではなく、

```text
Shared molecular encoder
 ├─ CYP1A2 direct pIC50
 ├─ CYP2C9 direct pIC50
 ├─ CYP2D6 direct pIC50
 ├─ CYP3A4 direct pIC50
 ├─ single-concentration response × CYP
 ├─ Emax × CYP
 ├─ +NADPH/TDI-condition pIC50
 ├─ ΔpIC50 / TDI probability
 ├─ substrate probability
 ├─ parent depletion
 └─ atom-level SoM
```

のような構成を優先します。

ただしloss weightsは件数比例にせず、direct pIC50をprimary endpointに置きます。大きなTox21/ToxCast datasetにencoderが支配されないよう、task-balanced samplingやgradient balancingが必要です。

薬物動態領域では、pretraining＋multitask/transfer learningにより小規模PK endpointsのgeneralizationを改善する例もありますが、これはCYP direct pIC50そのものの実証ではなく隣接領域からの支持です。citeturn56academia9

### single concentrationは弱ラベルではなく「別の観測方程式」と考える

濃度 \(C\) におけるresponseは、単純化すればHill-type relationshipでIC50、Hill slope、Emaxに依存します。

そのため、

\[
r(C)\rightarrow pIC50
\]

を直接regressするより、

\[
\{pIC50,\,E_{\max},\,h\}
\rightarrow r(C)
\]

というgenerative observation modelを入れる方が自然です。

OpenADMETにはEmaxも公開されたため、このモデルを実データで校正できる点が非常に大きいです。citeturn14search3

実装を簡単にするなら、

1. shared encoder  
2. pIC50 head  
3. single-concentration head  
4. Emax head

だけでも十分です。

より進めるなら、latent pIC50/Emax/slopeからsingle-concentration responseを再構成するphysics-informed auxiliary lossを入れられます。

### TDIはbinary classificationよりpaired-potency modelが情報量を使える

OpenADMET TDIでは同一化合物を、30分preincubationで

- NADP(H)なし
- NADP(H)あり

のpaired armsとして測定し、\(\Delta pIC50\)を計算しています。competitionでは2-fold shift thresholdがTDI labelに使われています。citeturn17view0turn18view0

したがって最初から

```text
is_TDI = 0/1
```

だけを学習するより、

\[
pIC50_{\mathrm{direct}},
\quad
pIC50_{\mathrm{+NADPH}},
\quad
\Delta pIC50
\]

をjoint predictionし、そのposteriorからTDI probabilityを出す方が情報を多く使えます。

特にboundary近傍ではhard labelが不安定です。公開記事でもcredible interval/posteriorを用いたより高度なTDI label formulationが検討されています。citeturn18view0

また同じ記事が明記するように30分preincubationだけではslow inactivatorを取り逃す可能性があり、IC50 shiftだけではreversible/covalent/metabolite-mediated mechanismを完全に区別できません。したがってTDI labelを「機構のground truth」と考えるのは危険です。citeturn18view0

### CYP3A4 TDIではsubstrate/probe依存性を意識した方がよい

CYP3A4はallosteryやmultiple binding modesの影響が特に大きく、HLM-based TDI assessmentでassay substrateによりfalse-positive/overpredictionが起こり得ることを調べた研究があります。citeturn14search12

小規模なCYP3A4 TDI SAR seriesやkinact/KI実験のsupplementary datasetsも公開されています。これらは数千件を増やす用途ではありませんが、

- modelが正しいmechanistic series rankingを出すか
- direct inhibitorとTDIを分離できるか
- metabolizable structural changeに反応するか

という**mechanistic challenge set**として有用です。citeturn14search6turn14search13

### 半教師あり学習は「pseudo-label大量投入」よりrepresentation learningを優先

unlabeled moleculesを数百万集めてpseudo-pIC50を付ける方法は、一見データを増やせますが、現在のラベル数1,300–2,300/isoformならteacher biasを増幅しやすいです。

優先順としては、

\[
\text{unlabeled molecular pretraining}
>
\text{consistency regularization}
>
\text{high-confidence pseudo-labeling}
\]

を勧めます。

特に公開vendor-like chemical space、PubChem structures、Tox21/ToxCastの構造を使ってencoderを事前学習し、CYP taskでfine-tuneする方が安全です。

pseudo-labelを使う場合は、training distribution内だけでなく、**scaffold/OOD distanceとpredictive uncertaintyの両方**で採択するべきです。

この部分はCYP-specificに強い実証を今回確認できたわけではなく、一般的small-molecule MLからの類推です。

### 文献抽出は「全文LLM抽出」より二段階にする

文献・特許・regulatory文書には未収載のIC50/Ki/TDI値が残っています。PK drug–drug interaction文献については、abstract/sentenceから実験的PK evidenceを高精度で絞り込むtext-mining研究があり、文献triage自体は自動化可能です。citeturn54academia9

CYPでは次の二段階が現実的です。

```text
Stage A: 高recall文書検索
"CYP3A4" AND ("IC50" OR "Ki" OR "kinact" OR "time-dependent")
"CYP2D6" AND ...
```

その後、

```text
Stage B: compound / CYP / value / unit / relation /
        substrate / enzyme source / preincubation /
        NADPH / reference をschema抽出
```

し、最後に人手でverificationします。

特にmedicinal chemistry特許では、

`IC50 < 1 µM`  
`1–10 µM`  
`>30 µM`

のようなbin/censored ADME tablesが多いため、数値に変換せずinterval labelsとして使う方がよいです。

BindingDBはliterature/patent extractionも行っていますが、ChEMBL/selected PubChemも取り込むため、**BindingDB検索を「発見器」として使い、原典へ戻る**のが最も安全です。citeturn55search4turn55search18

### 少数だけ追加測定できるなら、active learningの目的を「最高活性探索」にしない

追加測定が例えば各CYP 100–300化合物程度しかできない場合、predicted strongest inhibitorsを測るだけでは既知SARを確認するだけになりやすいです。

優先すべき候補は、

\[
\text{acquisition}
=
\alpha\,\text{uncertainty}
+\beta\,\text{chemical diversity}
+\gamma\,\text{assay-disagreement}
+\delta\,\text{OOD coverage}
\]

です。

特に今回の環境では「assay-disagreement」が使えます。

たとえば、

- OpenADMET modelはstrong、Tox21 modelはweak
- single concentrationはstrongだがpredicted pIC50はweak
- directはstrongだがTDI branchはuncertain
- substrate/SoM modelはstrongly metabolizableだがdirect modelはinactive

という化合物を意図的に測ると、**モデルのどの仮説が間違っているかを判別できる実験**になります。

OpenADMET自身も一次screen hitからdose-response follow-up、さらにhit近傍のchemically similar compoundsを追加する設計を採っています。citeturn17view0turn18view0  
一般のchemical-space active learningでも、uncertaintyだけでなく構造–物性仮説と探索を組み合わせる研究があります。citeturn56academia7

ここは**実測情報を本当に増やす方法**です。

direct inhibitionなら、追加100化合物をすべてランダムに選ぶより、

- 30%：high uncertainty / in-domain
- 30%：underrepresented scaffolds
- 20%：cross-assay disagreement
- 10%：near-threshold/censored領域
- 10%：known mechanistic controls

のようなbalanced designを試す価値があります。割合自体は提案であり、CYPで最適と実証済みの数字ではありません。


## 優先して取得するデータと、優先して検証する戦略

最終的な優先順位は次のように考えます。

### データ取得・確認の優先順位

| 優先度 | 対象 | 何をするか | 理由 | 主な不確実性 |
|---|---|---|---|---|
| **最優先** | **EPA Tox21 10K P450 endpoints** | invitroDB/CompToxから4 CYPの全response/AC50/hitcall/QCを抽出し、ChEMBL36・PubChem既取得分とInChIKey/AID/sourceで照合 | 9,667-compound規模の別測定campaignなら最大級の増分候補 | 実際のusable AC50件数、ChEMBL36重複率を今回確定できていない |
| **最優先** | **OpenADMET single-conf + Emax + direct + TDIの統合** | compound-level joint tableを作り、未測定ではなくmissing-by-selectionとして扱う | 新規DBを探す前に既存情報量を大幅に回収できる | model設計によるgainは要実証 |
| **最優先・将来** | **Challenge終了後のwell-level signals / pIC50-slope posterior** | 11月3日以降、公開/申請条件を再確認 | label uncertainty、censoring、curve qualityを直接使える | 実際の提供形式・時期はChallenge後まで未確定 citeturn18view0 |
| **高** | **OpenADMET March reaction phenotyping** | CYP3A4/CYP2J2 depletionを現在の4,905 structuresとcrossmatch | substrate/metabolism情報として非常に相補的 | 同一library由来の重複が大きい可能性 |
| **高** | **ToxCast NVS_ADME** | Novascreen CYP endpointsを別taskで追加 | Veith/OpenADMETとは別platform、環境化学空間も補完 | small-drug-like target domainとのshift |
| **高** | **PubChem AID1851の完全provenance化** | AC50だけでなくcurve class、inactive、quality、SID/CIDまで保持 | 既知データでも使い方を改善できる | 新規構造数はほぼ増えない可能性 |
| **中～高** | **BindingDB patent/literature差分** | ChEMBL `SRC_ID=37` とまず差分化し、原典単位で抽出 | ChEMBL未収載文献・特許の発見器として有効 | 再配布重複が非常に多い |
| **中** | **CarbonMangels substrate / ATTNSOM SoM** | 原典provenance確認後auxiliary/pretrainingへ | CYP recognitionを別角度から学べる | direct pIC50への改善量は未知 |
| **中** | **SABIO-RK** | CYP名＋human＋inhibitor/substrateでKi/Km/conditionを抽出 | 少量でもmechanistic metadataが濃い | CYP薬物代謝でのcoverageは要確認 |
| **中・評価用** | **FDA clinical index sets** | strong/moderate inhibitors、sensitive substrates等をexternal panel化 | clinical relevanceのsanity check | 定量in-vitro regression datasetではない citeturn55search9 |
| **条件付き高** | **Certara DIDB** | 予算があればTDI/DDI原典探索に利用 | curated DDI resourceとして特にTDIに魅力 | 商用、export可能範囲・ライセンス要確認 citeturn55search23 |

BindingDBについては、**「BindingDBを全部追加」ではなく「ChEMBL36にないBindingDB-native recordsだけを探す」**のが重要です。BindingDB→PubChem、ChEMBL→PubChemという流れもあるため、PubChem側でも必ずDepositor/Sourceを見てください。citeturn55search4

またOpenADMETの利用条件には一つ確認事項があります。Hugging Face dataset card側にはApache-2.0 tagが見える一方、OpenADMETの記事では**code=Apache 2.0、data=CC BY 4.0**と明記されています。再配布・派生dataset公開を行う場合は、リポジトリ内の最新LICENSEとデータ個別metadataを固定versionで確認するのが安全です。citeturn14search3turn18view0

### データ不足対策の優先順位

| 優先度 | 手法 | 実測を増やすか | CYPでの根拠 | 推奨理由 |
|---|---|---|---|---|
| **最優先** | **assay-aware multi-task / multi-fidelity** | いいえ | cross-CYP PCM、OpenADMETのpaired assays | 4 CYP、single-conf、TDI、Emaxを自然に統合できる |
| **最優先** | **censored + uncertainty-aware regression** | いいえ | OpenADMETがCI・低pXC50不確実性を明示 | `<4`等を捨てず、ラベル品質を直接利用 |
| **最優先** | **selection-biasをsingle-confで補正** | いいえ | OpenADMETのfollow-up selection design | 現在のpIC50集合のsampling biasを直接扱える |
| **最優先・TDI** | **paired direct/+NADPH potency model** | いいえ | OpenADMET IC50-shift assayそのもの | binary TDIより情報量が多い |
| **高** | **external assaysは同一labelにpoolせずassay head分離** | いいえ | Veith/OpenADMET双方がprobe/platform effectsを示す | データ数を増やしながら意味を壊しにくい |
| **高** | **substrate/depletion/SoM auxiliary learning** | いいえ | CYP cross-isoform/SoMで実証あり | binding＋metabolismに関係するrepresentationを補強 |
| **高** | **small active-learning experiment** | **はい** | OpenADMET自身がscreen→follow-up selectionを採用 | 少ない追加測定を最大限情報化できる |
| **中** | **molecular self-supervised pretraining** | いいえ | CYP直接の決定的実証は今回未確認 | low-riskだがgainはarchitecture依存 |
| **中** | **high-confidence pseudo-labeling** | いいえ | CYP-specific evidence弱い | teacher bias増幅のリスクあり |
| **中～高** | **文献・特許のschema-based extraction＋human verification** | 既存実測を回収 | DDI文献miningは実証あり | 新規TDI/Ki/IC50を得られる可能性 |
| **機構検証** | **docking/reactivity/SoM constraints** | いいえ | CYP3A4/2D6 SoMでmechanistic evidence | pure potency predictorよりsanity checkとして有効 |

私なら最初のモデル実験は、**「外部データを大量投入したモデル」ではなく、現在データだけで以下のablationを先に行います。**

\[
\begin{aligned}
M_0 &: \text{direct pIC50 only}\\
M_1 &: M_0+\text{4-CYP multitask}\\
M_2 &: M_1+\text{single-concentration head}\\
M_3 &: M_2+\text{Emax}\\
M_4 &: M_3+\text{censor/CI-aware loss}\\
M_5 &: M_4+\text{TDI paired heads}
\end{aligned}
\]

これで**「本当に足りないのは化合物数なのか、現在ある不完全ラベルの利用法なのか」**を切り分けます。

その後、

\[
M_5+\text{Veith},
\quad
M_5+\text{Tox21},
\quad
M_5+\text{ToxCast},
\quad
M_5+\text{substrate/SoM}
\]

を一つずつ追加します。

評価はrandom splitだけでは不十分です。OpenADMETではhitからchemical-similarity expansionが行われているため、近縁analogがtrain/testに跨るrandom splitは過大評価し得ます。citeturn17view0turn18view0  
最低でもscaffold splitに加え、可能なら**「元hit＋そこから追加された類似化合物」をgroupとして分離するsplit**を作るべきです。

最も優先したい研究仮説は、結局次の三つです。

**第一は、Tox21 10K P450 dataがChEMBL36に対して実際にどれだけ独立した新規compound-assay pairsを与えるか。**  
規模が大きく、今回確認できた候補の中では最大の増分可能性があります。ただしusable row数と重複率が未確認なので、不確実性も最も大きい候補です。citeturn41search1turn41search10

**第二は、single-concentration、Emax、direct pIC50、TDIを一つのmulti-fidelity modelにした時の改善量。**  
これは新規データ取得なしで検証でき、現在のpIC50 datasetが一次screenで選択された集合である以上、selection biasを補える可能性もあります。OpenADMETにはそのために必要な同一campaign内の観測が非常に揃っています。citeturn14search3turn18view0

**第三は、外部AC50/IC50を「追加pIC50」として混ぜるのではなく、assay-conditioned transferにした時に初めて性能が上がるか。**  
Veith自身がsubstrate interferenceとprobe依存性を、OpenADMET自身がfluorescence interference・substrate dependence・2D6 platform変更を報告しているため、CYPではこの問題は理論上の懸念ではなく実測上の問題です。citeturn32search1turn18view0

この三つが陽性なら、現在の1,300–2,300件/isoformというラベル数は「単純に少なすぎる」のではなく、**周囲にある数万件の異種観測をpotencyと同一視せず、来歴とassay contextを保持したまま利用できていないことが主なボトルネック**である可能性があります。逆にここで改善しなければ、次に投資すべき対象はpseudo-labelingではなく、Tox21差分・未収載文献/特許・少数のactive-learning追加測定です。