# CYPの機能・基質性・阻害・代謝を理解するための dry 研究ランドスケープ

## 全体像と結論

CYP の dry 研究を広く見渡すと、最も重要な結論は「**CYP 予測を一つの分子物性予測問題として扱わない方がよい**」ということです。CYP が関与する現象は、少なくとも **局所化学反応性 → 酵素への進入・結合 → 反応可能な配向 → 酸化反応 → 生成物分岐 → 酵素失活／再生 → assay で観測される値 → in vivo の DDI** という階層からなります。P450 の酸化化学の中心には高反応性の Compound I があり、その存在と C–H 活性化能は実験的にも直接特徴づけられていますが、実際の医薬品代謝では「反応しやすい原子がある」だけでは不十分で、その原子がヘム近傍にどの程度の確率で配置されるかが競合します。citeturn42search0turn18academia38

この観点から見ると、現在もっとも成熟している dry 問題は **site of metabolism（SoM）の候補順位付け**です。これは局所的な反応性が候補原子を強く絞るためで、SMARTCyp の「intrinsic reactivity + accessibility」、RS-Predictor の量子・トポロジカル記述子、XenoSite/FAME 系の learned local patterns、さらに近年の cross-isoform 学習へと、一貫して「局所化学」と「酵素依存性をどう足すか」が中心テーマになっています。citeturn43view0

一方、**IC50、基質／非基質、TDI、実際に観測される代謝物量**は、分子構造だけから定義される intrinsic property ではありません。例えば大規模 CYP qHTS 研究では CYP1A2、2C9、2C19、2D6、3A4 に対して 17,143 化合物の potency が測定されていますが、その値は recombinant enzyme、probe substrate、assay format という観測系の上にあります。したがって「同じ分子の CYP3A4 IC50」という表記でも、probe、基質濃度、pre-incubation、protein concentration、NADPH、microsome/recombinant enzyme などが違えば、別の観測量とみなす方が安全です。citeturn30search0turn30search26

特に **CYP3A4 は静的な単一 docking pose で説明しにくい**代表です。X 線構造から大きな ligand promiscuity と active-site plasticity が示され、速度論でも複数基質が同時に存在するモデルが必要になるケースがあります。単純な Michaelis–Menten 式を無条件に当てはめると、そこで推定された \(K_m\)、\(V_{\max}\)、阻害定数を「物理的な単一状態」と誤解する危険があります。citeturn31search6turn33search34turn33search3

逆に、**CYP2D6 と CYP2C9 は、構造計算が wet perturbation と結びついた好例**があります。CYP2D6 では計算モデルが Glu216 を含む酸性残基の役割を提案し、変異導入で basic substrate の \(K_m\) が 10–100 倍変化し、Glu216/Asp301 の電荷を同時に中和すると通常 CYP2D6 の基質ではない diclofenac に対する活性が 1000 倍超増えるという、かなり強い因果的検証が行われました。これは「attention map や docking pose がもっともらしい」という意味の説明とは質が違います。citeturn40view5

したがって、今後の研究をモデル種類ではなく「どの物理階層を解こうとしているか」で整理すると、次のようになります。

| 問題 | dry で比較的よく解ける部分 | 追加情報が特に効く部分 | 本質的な難所 |
|---|---|---|---|
| 可逆的 CYP 阻害 | chemotype、電荷、疎水性、ヘム配位 motif、類縁体 SAR | assay condition、probe、3D pose、pKa | IC50 自体が assay-dependent、3A4 の multiple binding |
| 基質性 | enzyme affinity に寄与する物性・motif、既知基質との類似性 | SoM/reactivity、productive pose、turnover data | 「結合する」と「代謝される」が別 |
| SoM | 局所反応性、atom environment | isoform、steric accessibility、ensemble docking | minor site、isoform specificity、未観測 site を negative にしがち |
| 代謝物 | SoM + reaction rule | QM reaction mechanism、逐次代謝 | branching、relative abundance、非 CYP 反応 |
| TDI/MBI | reactive-metabolite alert、親化合物の binding | QM、bioactivation pathway、\(K_I,k_{\rm inact}\)、adduct data | 多段階反応を一つの binary label に潰すこと |
| in vivo DDI | in vitro parameter の統合 | PBPK、unbound exposure、enzyme turnover | CYP だけでなく transport/induction/organ exposure が関与 |

この整理からは、既に試されている fingerprint/tree、CheMeleon、TabPFN、fine-tuning を単に別 backbone に置き換えるより、**観測条件・反応階層・cross-isoform 情報を明示的にモデル化すること**の方が、新しい情報を取り込む余地が大きいように見えます。これは以下の文献群からの総合的な推論です。citeturn30search0turn43view0turn18academia38turn33search3

## 予測対象ごとの到達点と難しさ

### 可逆的阻害強度は「分子の IC50」を予測しているわけではない

CYP inhibition QSAR/ML は実用上かなり有用ですが、機構的にはラベル定義が最初の難所です。単純な competitive inhibition なら、基質濃度と \(K_m\) が変われば IC50 と \(K_i\) の関係も変わります。さらに mixed inhibition、substrate inhibition、複数 ligand binding があれば単純な変換式自体が成立しません。CYP3A4 では multiple-substrate kinetic schemes が実験データを説明するために必要になる例があり、これは probe-dependent な CYP3A4 inhibition を単一 scalar とみなすことへの直接的な警告です。citeturn33search34turn33search3turn30search26

したがって、IC50 モデルが高精度でも説明しているものはしばしば

\[
\text{structure}
\longrightarrow
\text{assay-specific apparent inhibition}
\]

であって、

\[
\text{structure}
\longrightarrow
K_d \ \text{or microscopic binding free energy}
\]

ではありません。この区別は 1A2/2D6/2C9 の比較的特徴的な pocket では多少弱まり、3A4 では特に重要になります。CYP1A2 の結晶構造は比較的狭く planar/aromatic ligand を収容する構造的特徴を示し、CYP2D6 は酸性残基による basic substrate 認識、CYP2C9 は Arg108 を含む acidic substrate 認識が実験的に支持されています。citeturn30search1turn40view5turn31search4

ここで 3D 構造は「全面的に 2D ML より優れる」のではなく、**特定の相互作用仮説が存在するときに最も価値が高い**です。例えば「CYP2D6 の protonated amine と acidic residue の組」「CYP2C9 の carboxylate と Arg108」「ヘム鉄への type-II coordination」のような hypothesis がある場合、pKa/charge state + docking/MD + matched-pair SAR は解釈可能性を大きく上げられます。CYP2D6 の結晶構造も Asp301、Glu216、Phe483、Phe120 を substrate recognition と結びつけています。citeturn31search16turn40view5

実装上は、**複数 assay の IC50 を一つに平均しないこと**が非常に重要です。compound × CYP × probe × enzyme system × preincubation × substrate concentration × endpoint を tensor とみなし、condition を covariate/embedding として持つ方が、単純な compound→pIC50 回帰より因果構造に近くなります。これは文献から直接確立済みの「唯一の正解」ではなく、probe selection や CYP3A4 kinetics に関する実験知見から導かれる設計提案です。citeturn30search26turn33search3

### 基質性は binary classification より turnover の分解が重要

「substrate / non-substrate」は dry データセットを作りやすい一方で、かなり粗いラベルです。ある分子が active site に結合しても productive orientation を取れなければ代謝されず、逆に複数の弱い pose の一部だけが productive でも measurable turnover が生じます。Feng と Friesner の研究が intrinsic chemical reactivity と induced-fit docking を二段階に分けたのはまさにこのためで、反応しやすさと protein geometry の両方を明示的に扱っています。彼らの個別 test setting では top-2 SoM が CYP3A4 で 92.7%、CYP2D6 で 100% と報告されていますが、これは後述する XenoSite/ATTNSOM の数値とは dataset・split・endpoint が違うためランキングには使えません。citeturn18academia38

基質性については、私は次の分解の方が研究上有益だと考えます。

\[
P(\text{turnover})
\approx
P(\text{bind/access})
\times
P(\text{productive pose}\mid \text{bound})
\times
P(\text{reaction}\mid \text{pose})
\times
P(\text{product escapes})
\]

これは厳密な速度式ではなく conceptual decomposition ですが、2D 表現、3D docking、QM をどこに使うべきかが明瞭になります。SMARTCyp 系の intrinsic reactivity、protein-aware docking、QM/MM はそれぞれ異なる項を補っていると理解できます。citeturn43view0turn18academia38turn31search15

また「文献に substrate と書かれていない」を negative とみなすのは危険です。CYP substrate annotations は測られた isoform が化合物ごとに異なるため、これは典型的な **missing-not-at-random** に近い問題です。後述する ATTNSOM の「未測定 isoform を negative にせず mask する」という発想は、SoM だけでなく substrate/inhibition multitask にもそのまま一般化する価値があります。citeturn43view0

### SoM は成熟しているが、何を「正解」とするかに注意が要る

SoM は CYP dry research の中で最も体系化されています。SMARTCyp は intrinsic reactivity と steric accessibility、RS-Predictor は atom-level topological/quantum descriptors、XenoSite は atom/molecule descriptors を neural network に統合し、FAME 3 は atom-centered circular fingerprints と extremely randomized trees を用いました。近年は protein/isoform information と cross-isoform relationships を組み込む方向へ進んでいます。citeturn43view0

この歴史はかなり示唆的です。**SoM では「どの GNN が一番か」より、局所反応性という inductive bias が長く生き残っている**からです。反応可能な C–H、N、O、芳香環等を chemistry がまず限定し、その中から protein accessibility が選ぶ、という二段構造は SMARTCyp から modern hybrid methods まで繰り返し現れています。Compound I の反応化学が実験的に直接検証されていることも、この分解に物理的根拠を与えています。citeturn42search0turn43view0turn18academia38

2026 年の ATTNSOM は、この問題に対する最新の興味深い例です。XenoSite 由来データは 679 分子、2,003 molecule–isoform instances、9 CYP（1A2/2A6/2B6/2C8/2C9/2C19/2D6/2E1/3A4）からなり、1 分子当たり metabolized atom は平均約 1.4 個なので atom-level では強い class imbalance があります。モデルは shared graph encoder、molecule-conditioned representation、cross-isoform attention を組み合わせ、**実験されていない molecule–isoform pair を auxiliary loss から mask** しています。citeturn43view0

ただし ATTNSOM は leaderboard よりも **evaluation design の教材**として読む方が価値があります。論文の Zaretzki benchmark は random 10-fold CV で、ATTNSOM の平均 Top-2 は 0.871、過去の XenoSite の報告値は 0.870 ですが、XenoSite は leave-one-out なので直接比較ではありません。さらに ATTNSOM の著者ら自身が GLMCyp を「報告された evaluation set に training/test instances が混在するため fair comparison から外した」と記載しています。citeturn43view0

さらに公開コードを追うと、ATTNSOM はまず **(molecule, isoform) graph を単位に StratifiedKFold** し、同じ molecule が別 isoform の training graph として残り得ます。一方で auxiliary cross-isoform labels については held-out (molecule, CYP) pair を明示的に mask して漏洩を防いでいます。つまり「同じ test pair の label を直接見る」漏洩対策は入っていますが、**完全に未知の新規 molecule へ一般化する評価とは別物**です。新規化合物予測が目的なら canonical parent molecule 単位で全 isoform を group holdout する再評価が重要です。これはコードから確認できる点です。fileciteturn6file0L2-L2 fileciteturn7file0L2-L2

公開実装は以下です。

`https://github.com/dmis-lab/ATTNSOM` fileciteturn1file0L2-L6

### 実代謝物は SoM より一段難しい

SoM が分かっても、そこから **どの transformation が起こるか**、さらに **どの metabolite が主要になるか** は別問題です。例えば同じ atom 周辺でも hydroxylation、dealkylation、heteroatom oxidation、epoxidation といった反応候補があり、その後さらに secondary metabolism や phase-II conjugation が続き得ます。ATTNSOM もこの点では SoM prediction であって metabolite structure/yield prediction ではありません。citeturn43view0

この領域で古くから有効なのが、**reaction rules と empirical frequency の組み合わせ**です。SyGMa は expert-defined transformation rules と empirical scoring で potential metabolites を生成する方法として開発され、その後 Python reimplementation も公開されています。これは「すべてを end-to-end reaction model に学習させる」方法とは対照的ですが、rare reaction を chemical constraint で扱えること、生成経路を人間が監査できることが依然として大きな長所です。citeturn39search1

公開実装:

`https://github.com/3D-e-Chem/sygma` — Python library for drug metabolite prediction。fileciteturn2file0L2-L6

GLORYx は SoM prediction と literature-derived reaction rules を統合し、Phase I/II metabolites を生成する方向の公開ツールです。したがって実用上は、**SoM ranking → reaction-class prior → chemically valid product enumeration → downstream ranking** という modular pipeline は、現在でも非常に合理的です。citeturn39search3

特に注目すべき新しい外部データが AstraZeneca の MetID データです。2025 年に公開されたデータには 120 compounds の metabolite-identification schemes が含まれ、atom-resolved な exact SoM annotations が提供されています。ATTNSOM はこれを AZ-ExactSOM として使っていますが、重要なのはこのデータが **human hepatocyte の aggregate metabolism で CYP isoform-specific ではない**ことです。したがってこれは「CYP3A4 SoM の外部検証」ではなく、より現実的な「human cellular metabolism の metabolic soft spot」検証と解釈すべきです。citeturn42search1turn42search3turn43view0

この distinction は非常に重要です。isoform-resolved recombinant CYP data で学習したモデルを hepatocyte data に当てて性能が落ちた場合、それは必ずしも「モデルが失敗した」のではなく、**CYP isoform label → cellular metabolic outcome という domain shift**を測っている可能性があります。逆に、この shift を explicitly model できれば、薬物動態上はこちらの方が価値のある研究になります。citeturn43view0turn42search3

### TDI は「阻害分類」より反応ネットワークとして扱うべき

Time-dependent inhibition（TDI）/mechanism-based inactivation（MBI）は、通常の reversible inhibition よりはるかに multi-step です。典型的には

\[
E + I
\rightleftharpoons EI
\rightarrow EI^*
\rightarrow E_{\rm inactive}
\]

のように、binding に続いて CYP 自身が inhibitor を reactive species へ変換し、その後 enzyme inactivation が起こります。標準的な近似では

\[
k_{\rm obs}
=
\frac{k_{\rm inact}[I]}
{K_I+[I]}
\]

を複数濃度・複数 preincubation time から推定しますが、実系では inhibitor depletion、multiple metabolites、reversible binding、enzyme loss などが重なり、この単純形から外れるケースがあります。progress-curve 全体を直接 mechanistic model に fitting する研究が進んでいるのはそのためです。citeturn34search4turn34search7

したがって「TDI yes/no」を structure-only classifier で学習すること自体には screening value がありますが、説明できるのは **TDI-associated chemotype の統計的 propensity**までです。2025 年にも reversible/TDI CYP3A4 を対象とする (Q)SAR 研究が報告されていますが、mechanism の理解を目的にするなら、次の latent quantities に分ける方が情報量があります。citeturn34search1

\[
\text{parent binding}
\rightarrow
\text{bioactivation}
\rightarrow
\text{reactive intermediate}
\rightarrow
\text{enzyme adduct/heme damage}
\rightarrow
k_{\rm inact}, K_I
\]

つまり TDI では QM が特に有望です。普通の IC50 では QM barrier を計算しても assay potency に直接対応しないことがありますが、MBI では「どの intermediate が形成できるか」「どの bond activation が律速候補か」という化学機構そのものが endpoint の一部だからです。とはいえ QM だけで MBI を断定すべきではなく、NADPH dependence、time dependence、dilution/recovery、reactive-metabolite trapping、protein/heme adduct MS などと組み合わせて初めて強い機構証拠になります。progress-curve modeling の研究も、単純な \(K_I,k_{\rm inact}\) 抽出が常に十分ではないことを示しています。citeturn34search4turn34search7

Ritonavir/CYP3A4 はこの階層を考える重要な系です。CYP3A4–ritonavir complex の構造が決定され、その後 ritonavir analogues を合成して ligand–enzyme interactions を分解する研究が行われ、さらに pharmacoenhancer として cobicistat の SAR 開発へつながっています。これは「計算 pose → feature importance」で終わらず、**構造 → analog synthesis → biochemical potency/selectivity → pharmacological design** と往復した例として非常に価値があります。citeturn33search0turn31search2turn33search5turn33search2

## 構造・物性・反応機構の情報はどこで本当に効いたか

### 二次元構造と物性が強い理由

CYP の screening endpoint では、2D fingerprint や simple physicochemical descriptors が意外に強いのは不思議ではありません。lipophilicity、size、aromaticity、charge/protonation、H-bond pattern は binding/access と assay behavior の両方を proxy し、medicinal-chemistry series 内では local SAR が強いからです。Veith らの 17,143-compound qHTS のようなデータは、この種の statistical SAR を学ぶには非常に強力です。citeturn30search0

しかし、2D model の feature attribution から「この N が CYP2D6 の Asp301 と相互作用したから」という mechanistic statement へ飛ぶことはできません。同じ correlation が pKa、lipophilicity、series membership で説明できる場合があるからです。CYP2D6 の Glu216/Asp301 のケースが強いのは、計算モデルの prediction に対し site-directed mutagenesis という intervention が行われ、substrate specificity と regioselectivity が大きく変わったからです。citeturn40view5

したがって既存の descriptor/fingerprint、CheMeleon、TabPFN 等をさらに有効に使う方向としては、backbone を変えるより **mechanistically chosen auxiliary variables を出力させる**方が興味深いです。例えば pKa/charge class、candidate SoM count、local reactivity、heme-coordination propensity、isoform-specific anchor compatibility などを latent/intermediate task とし、最終 IC50 と一緒に学習させる設計です。これは既存文献から直接 benchmark 済みの処方ではなく、CYP2D6/2C9 の perturbation evidence と SoM methods の構造から導く研究仮説です。citeturn40view5turn31search4turn43view0

### 三次元構造は「pose を出す」より competing hypotheses を作るために使う

CYP1A2、2C9、2D6、3A4 の結晶構造を見ると、同じ structure-based workflow を全 isoform に適用する合理性は低いことが分かります。

**CYP1A2** は最初のヒト構造研究から polycyclic aromatic substrates に適応した比較的狭い active site が示され、平面性・aromatic shape complementarity が特に意味を持ちやすい isoform です。従って 1A2 では 3D shape、planarity、aromatic packing を 2D model の補助情報にする価値が比較的高いと考えられます。citeturn30search1turn30search3

**CYP2C9** では flurbiprofen-bound structure と Arg108 周辺の研究が重要です。Arg108 の置換は diclofenac や S-warfarin のような acidic substrates の代謝を大きく低下させ、neutral substrate では影響が異なることが示されています。最近の MD + kinetics の研究でも diclofenac carboxylate–Arg108 salt bridge や aromatic contacts を検討し、mutation experiments で計算仮説を検証しています。citeturn31search4turn31search11

**CYP2D6** では protonatable basic nitrogen と acidic active-site environment の関係が特に強く、Glu216/Asp301 を変えると substrate class 自体を変えられるほどの effect が得られました。ここでは docking electrostatics、protonation-state enumeration、pKa calculation が単なる「extra descriptors」ではなく enzyme specificity の mechanistic variables になります。citeturn40view5turn31search16

**CYP3A4** は逆です。large/flexible pocket、multiple ligand binding、cooperativity があるため「最高 docking score の一 pose」を explanation と呼ぶ危険が最大です。CYP3A4 では ensemble docking、multiple occupancy、MD、kinetic model をセットにしないと、static pose に対する過剰解釈になりやすいです。結晶構造そのものが ligand promiscuity を示し、速度論研究も multisubstrate models の必要性を示しています。citeturn31search6turn33search34turn33search3

このため structure-based features を ML に追加する場合、単一 docking score よりも

\[
P(d_{\text{SoM-heme}}<d_0),\quad
\text{pose entropy},\quad
\text{fraction of productive frames},\quad
\text{anchor occupancy},\quad
\text{tunnel accessibility}
\]

のような **ensemble quantities** の方が CYP biology に近いと考えられます。これは特に 3A4 で検証価値が高い仮説です。citeturn31search6turn42search20

### 膜を無視した構造計算にも限界がある

薬物代謝 CYP は membrane-associated protein であり、可溶化した catalytic domain の crystal structure だけが生理的状態ではありません。membrane-bound CYP2C9 の MD では membrane environment と active-site access tunnels の opening/closing の関係が調べられており、CYP2C9 と高相同性の CYP2C19 でも膜との相互作用が異なり得ることが simulations で示されています。citeturn42search20turn42search14

これは「2C9 と 2C19 は sequence が似ているので protein embedding も似るはず」という単純な multitask transfer に対する重要な反例です。**sequence similarity が高くても membrane orientation、tunnel dynamics、local electrostatics が substrate specificity を分ける可能性**があります。したがって protein embedding を入れるなら、それを mechanism と同一視せず、membrane/structure/dynamics が持つ追加情報が残っているかを ablation する価値があります。citeturn42search14turn42search20

### QM は「何でも精度向上」ではなく reaction-limited な問題に使う

QM/QM-MM は計算コストが高いため、「全 compound に高精度 energy を計算して tabular model に足す」という使い方より、**reactivity が bottleneck の候補間の discrimination**に使うのが筋が良いです。SMARTCyp/RS-Predictor 以来、local reactivity は SoM prediction に有用であり、Compound I の C–H activation そのものも mechanistic spectroscopy/kinetics で裏づけられています。citeturn43view0turn42search0

特に有望なのは、

1. docking/ML が複数の candidate SoM を同程度に挙げる、
2. geometry はいずれも productive、
3. それでも metabolite ratio が大きく違う、

というケースだけ QM barrier で rerank する使い方です。Feng/Friesner の reactivity + induced-fit docking は、この「化学反応性と geometry の両方が必要」という思想の実証例です。citeturn18academia38

逆に CYP3A4 で計算した一 pose に対する barrier が実験の product ratio と違ったからといって、reaction theory が間違ったとは限りません。binding-state populations が間違っている可能性があるからです。**observed product flux は barrier × population の問題**と見るべきです。この切り分けこそ docking/MD と QM を組み合わせる本当の利点です。citeturn18academia38turn31search6

## wet–dry の往復から学べる「説明」の強さ

CYP 研究では、説明を少なくとも三段階に区別すると混乱が減ります。

**予測的説明**は「この atom/fragment が model prediction に寄与した」。  
**構造的仮説**は「この pose なら residue X と相互作用し、反応し得る」。  
**機構的説明**は「X を変えると、その予測された方向へ affinity、rate、regioselectivity が変わる」。

この中で三段目に到達した研究は少なく、だからこそ価値があります。CYP2D6 の Glu216/Asp301 研究はその典型です。computational modeling が Glu216 の役割を示唆した後、site-directed mutagenesis により bufuralol/dextromethorphan affinity が大幅に変わり、さらに acidic CYP2C9 substrate の diclofenac を CYP2D6 に代謝させる方向へ specificity を re-engineer できました。citeturn40view5

CYP2C9 の Arg108 系も同様です。flurbiprofen-bound crystal structure から acidic ligand recognition の構造基盤が見え、Arg108 mutations は S-warfarin/diclofenac turnover に大きく影響しました。さらに後続研究は MD で cooperative ligand arrangements や salt bridges を仮説化し、mutant kinetics で検討しています。citeturn31search4turn31search11

重要なのは、mutation の結果ですら「その residue と ligand の direct interaction を完全証明した」とは限らない点です。mutation は loop position、protein dynamics、heme environment まで変え得ます。したがって **crystal/MD + multiple mutations + matched substrates + kinetics** が揃ったときに証拠が強くなります。CYP2D6 の研究が Asp301 の単純な salt-bridge model に修正を迫り、Glu216 の役割を明らかにしたこと自体、最初の docking explanation を intervention で更新する好例です。citeturn40view5

CYP3A4–ritonavir 系は synthesis を介した別種の強い evidence です。ritonavir-bound structure に続いて analogues を合成し、どの ligand portion が CYP3A4 interaction に効くかを調べ、CYP3A-selective pharmacokinetic enhancer の設計へつなげています。cobicistat の開発も含め、これは **structure-guided explanation が chemical intervention に耐えたか**という観点で読む価値があります。citeturn33search0turn31search2turn33search5turn33search2

この考え方を現在の ML 系へ移植するなら、最も面白い prospective wet validation は「高 confidence prediction をただ測る」ことではありません。**モデル間で機構仮説が食い違う化合物を意図的に選ぶ**ことです。

たとえば、

- 2D model は代謝されると予測するが local reactivity は低い、
- QM は非常に reactive とするが docking/MD では productive pose が少ない、
- WT では高活性だが residue mutation により強く落ちると structure model が予測する、

という discordant examples を合成・測定すれば、どの物理階層が失敗しているかを一回の wet experiment で識別できます。これは単なる prospective accuracy test より研究情報量が高い、というのが本レビューからの提案です。citeturn18academia38turn40view5turn31search11

さらに hydrogen abstraction が疑われる site では isotope substitution も有力です。ATTNSOM の論文が紹介している ivacaftor/deuterated analogue の例のように、metabolically labile position を deuterate して pathway を実際に変えることは、SoM attribution を chemical intervention で検証する一つの方法です。citeturn43view0

## 小規模・不均一・欠測データをどう扱うべきか

### 欠測は negative ではない

CYP data の最大の統計的問題の一つは、同じ化合物について全 CYP・全 assay が測られていないことです。ATTNSOM はここに対してかなり明確な設計を採り、annotation のない isoform を auxiliary loss の negative とせず mask しています。Zaretzki dataset では一つの molecule が平均して複数 isoform に annotation を持つため、cross-isoform information を共有する余地があります。citeturn43view0

この考え方は SoM よりも、むしろユーザーの現在の用途である inhibition/substrate datasets に広げると面白いと思います。つまり

\[
Y_{\text{compound,CYP,probe,condition,endpoint}}
\]

という sparse observation tensor を作り、**missingness mask を独立に持つ**設計です。

その上で shared molecular representation は使ってよいものの、

\[
z_{\rm chemistry},
z_{\rm CYP},
z_{\rm assay},
z_{\rm endpoint}
\]

を分け、IC50、single-point inhibition、substrate turnover、SoM、TDI を condition-specific heads から出す方が、「全部まとめて補助タスク」にするより negative transfer を抑えやすいと考えます。これは ATTNSOM の cross-isoform missing-label handling と qHTS の assay heterogeneity から導く提案です。citeturn43view0turn30search0

### 単一濃度スクリーニングを auxiliary label ではなく censoring information として使う

ユーザーが既に試された「単一濃度 screening value による補助学習」は、別の使い方がかなり有望です。

例えば concentration \(c\) で inhibition fraction \(r\) が観測されたとき、それを別 task の target とするのではなく、

\[
r=f(c; IC_{50},h,\text{assay context})
\]

という dose-response model の一観測点として latent \(IC_{50}\) を constrain できます。

active/inactive だけならなおさら、正確な pIC50 label に変換するより **interval-censored observation** として扱う方が自然です。例えば最大 concentration まで inhibition が弱いなら「IC50 は測定範囲より右側」という情報であって、任意の numeric ceiling を真値として回帰する必要はありません。多濃度 qHTS が 17,143 compounds に存在することは、single-point → concentration-response latent model の calibration resource として興味深いです。citeturn30search0

これは現在試された auxiliary-learning を否定するものではなく、**same data を mechanistic observation model に組み込む別用途**です。特に screening concentration が複数 project/lab で違う場合、この方がデータを無理に同じ scale に押し込まず利用できます。

### assay condition は nuisance ではなく学習可能な変数

CYP assay の heterogeneity を単に batch correction で消すのも危険です。probe dependence や enzyme source difference は noise であると同時に、生物学的情報でもあるからです。特に CYP3A4 の multiple substrate behavior を考えると、probe identity は nuisance covariate ではなく mechanistic context です。citeturn30search26turn33search34

従って assay metadata が取れる場合は、compound embedding に加えて

\[
(\text{CYP},\ \text{probe},\ [S]/K_m,\ \text{enzyme source},
\ \text{preincubation},\ \text{NADPH},\ \text{protein})
\]

を明示入力し、「標準条件へ正規化された latent potency」と「各 assay での observation」を階層モデルで分ける価値があります。

例えば

\[
\theta_i = \text{latent compound–CYP interaction}
\]

\[
y_{ij} =
g(\theta_i,\text{condition}_j)
+\epsilon_{ij}
\]

という形です。

これは単一大規模 neural model である必要はなく、Bayesian hierarchical model、mixed-effects model、Gaussian process、mechanistic dose-response head と learned representation の hybrid でも構いません。むしろデータが小さい場合、condition effect を強く regularize できる方が有利でしょう。

### cross-isoform 学習では「何が cold なのか」を決める

ATTNSOM は cross-isoform task sharing が有効であることを示す興味深い evidence ですが、公開コードを見ると同一 molecule の異なる isoform instances が fold をまたぎ得ます。auxiliary label leakage 自体は mask で対処していますが、これは **cold pair prediction** と **cold molecule prediction** の違いを浮き彫りにします。citeturn43view0 fileciteturn6file0L2-L2 fileciteturn7file0L2-L2

CYP dry research では最低でも次の評価を分けるべきです。

| split | 問うもの | 実務上の意味 |
|---|---|---|
| random molecule–assay pair | 既知 chemistry と task sharing | 補完問題 |
| molecule-group split | 一度も見ていない分子 | 新規 compound |
| scaffold split | 新規 chemotype | hit-to-lead 外挿 |
| temporal/project split | 将来 series | discovery prospective に近い |
| leave-one-CYP/task-out | 新規 assay/task | transfer/cold-start |
| assay-condition holdout | 新規 protocol/lab | robustness |

この中でユーザーの目的が「今後合成する新規化合物の CYP1A2/2C9/2D6/3A4 を予測する」なら、**同一 canonical parent を全 task まとめて fold に入れる molecule-group split が最低条件**です。その上に Bemis–Murcko/scaffold、time/project split を重ねるのがよいでしょう。これは ATTNSOM のような cross-task setting では特に重要です。citeturn43view0

### uncertainty は training-domain distance と assay uncertainty を分ける

CYP の uncertainty は少なくとも二種類あります。

一つは「この chemistry を見たことがない」という **epistemic/domain uncertainty**。もう一つは「同じ compound を測っても assay/probe で値が変わる」という **measurement/context uncertainty** です。

これらを一つの predictive variance に混ぜるより、

\[
U =
U_{\rm chemistry}
+
U_{\rm assay}
+
U_{\rm model}
\]

と概念的に分けた方が、次の wet experiment を選ぶ active-learning policy に使えます。新規 scaffold で不確実なのか、既知 scaffold だが CYP3A4 probe dependence が大きいのかでは、必要な実験が違うからです。CYP3A4 kinetic complexity と assay dependence がこの区別の生物学的根拠です。citeturn33search3turn30search26

## 評価・失敗・再現性から見えてくること

CYP ML 論文の性能値を読むときは、**dataset name より split unit を先に確認する**のがよいです。ATTNSOM はこの問題をかなり率直に示しており、自身は Zaretzki data で random 10-fold CV、XenoSite の比較値は leave-one-out であることを明示しています。そのため平均 Top-2 0.871 と 0.870 は実質的な「勝敗」を意味しません。citeturn43view0

さらに同論文は、過去手法の一部は web service のみ／現在利用不能で、isoform-resolved atom scores が得られないため MCC の同一 split 比較ができない、と明記しています。これは古典 CYP benchmark の再現性問題をよく表しています。citeturn43view0

一方、2025 年の AstraZeneca human-hepatocyte 120-compound MetID dataset の公開は重要な前進です。既存 CYP-specific dataset と違う biological layer を測っているので「純粋な iid test」ではありませんが、だからこそ **external biological validation** として価値があります。isoform model が hepatocyte soft spots にどの程度 transfer するかを調べれば、「benchmark SoM を当てているだけか」「実細胞 metabolism に通じる signal を学んでいるか」を区別できます。citeturn42search1turn42search3turn43view0

Feng/Friesner 型の hybrid method も同じ注意が必要です。CYP3A4 92.7%、CYP2D6 100% という top-2 数字は強く見えますが、dataset、candidate-site definition、test construction が ATTNSOM と同じではありません。この論文から取るべき主な知見は「100% に近い精度」ではなく、**intrinsic reactivity だけでも docking だけでもなく、その組合せが regioselectivity を説明する**という mechanistic decomposition です。citeturn18academia38

SoM metrics 自体にも問題があります。Top-2 は medicinal chemistry 的には有用ですが、一分子に候補 atom が 20–30 個あると false-positive atoms を多数許容しても高くなり得ます。ATTNSOM が MCC を併記した理由も atom-level class imbalance とこの問題です。citeturn43view0

従って SoM は少なくとも **Top-k + atom-level PR/MCC + calibration**、metabolite は **top-N metabolite recall + exact structure + pathway rank**、inhibition は **continuous/censored error + calibration + clinically relevant threshold performance**、TDI は **binary flag だけでなく \(K_I,k_{\rm inact}\) parameter recovery**を分けて評価する方がよいでしょう。この最後の評価設計は文献を踏まえた提案であり、単一の community standard が確立しているという意味ではありません。citeturn43view0turn34search7

prospective validation ではさらに、「random で 50 compounds を測る」より **pre-registration 的に化学空間を層別化**するのが有益です。例えば known scaffold / new scaffold、high/low model agreement、2D/reactivity/docking が一致する領域／対立する領域を均等に選びます。こうすると overall RMSE だけでなく、「どの regime でどの approach が役立つか」という、ユーザーが求めている問いに答えられます。

そして applicability domain を報告するときは similarity だけでなく、CYP-specific mechanistic axis を足すべきです。例として 2D6 なら protonation/basicity、2C9 なら acidic functionality、1A2 なら planarity/aromatic shape、3A4 なら size/flexibility/multiple-binding propensity を binning し、各 regime で error/calibration を出す方法です。これらの axis には構造・変異研究による実験的根拠があります。citeturn30search1turn31search4turn40view5turn31search6

## 深く読む価値のある研究と、次に検証する価値のある問い

### 深く読む優先度が高い研究

**Paine et al., CYP2D6 Glu216/Asp301 mutagenesis, JBC 2003** は最優先です。理由は「計算モデルが residue hypothesis を出す → mutation → affinity/regioselectivity が変わる → substrate class そのものが変わる」という、dry explanation の理想形に近いからです。basic substrate の \(K_m\) が 10–100 倍変化し、double neutralization により diclofenac activity が >1000-fold 上がる結果は、feature attribution と causal mechanism の違いを考える上で非常に教育的です。citeturn40view5

**Wester et al. の CYP2C9–flurbiprofen structure と、その後の Arg108 re-engineering / MD–kinetics 研究**も同じ理由で重要です。特に acidic-substrate recognition を electrostatics、pose、mutation の三方向から見ることができ、CYP2D6 と対比すると「isoform-specific chemistry を表現にどう埋めるべきか」が具体的になります。citeturn31search4turn31search11

**Sansen et al., human CYP1A2 structure** は、1A2 に対して「なぜ planarity や aromatic topology が QSAR で効きやすいのか」を protein structure から理解するために読む価値があります。CYP3A4 と並べて読むと、static structure-based prediction が isoform ごとに異なる価値を持つ理由が分かります。citeturn30search1turn30search3

**Ekroos & Sjögren, CYP3A4 ligand promiscuity** と **multisubstrate kinetic models** はセットで読むべきです。前者は structural plasticity、後者はその kinetic consequence を示すため、3A4 を単一 pocket・単一 affinity model で扱う危険性を理解できます。citeturn31search6turn33search34turn33search3

**Sevrioukova & Poulos の CYP3A4–ritonavir structure と analogue studies、cobicistat SAR** は、structure biology と medicinal chemistry の往復を学ぶ材料として優れています。特に「もっともらしい binding explanation が analog synthesis に耐えるか」という観点で読むのがおすすめです。citeturn33search0turn31search2turn33search5turn33search2

**Rittle & Green, Compound I** は drug-specific predictor 論文ではありませんが、QM/reactivity descriptors が何を近似しようとしているのか理解するための基礎です。P450 chemistry の reactive state を実験的に捉え、C–H activation kinetics を直接調べています。citeturn42search0

**Feng & Friesner, reactivity + induced-fit docking** は、「reaction chemistry と protein geometry のどちらか一方では足りない」という dry strategy の代表として読む価値があります。数値 leaderboard より architecture of reasoning が重要です。citeturn18academia38

**Veith et al., 17,143-compound CYP qHTS** は、現代的 ML 論文ではありませんが、CYP inhibition data がどう作られ、どの程度 isoform overlap があり、大規模 screening endpoint が何を意味するか理解する基礎データとして重要です。citeturn30search0

**Chen et al./AstraZeneca 2025 MetID dataset** は今後特に重要になると思います。human hepatocyte 由来の 120 compounds と exact SoM annotations は規模こそ小さいものの、public ML benchmark と実際の cellular metabolism の間をつなぐ外部検証 resource だからです。citeturn42search1turn42search3

**ATTNSOM 2026** は「最新モデルだから」ではなく、cross-isoform sharing、missing-label masking、MCC、公開コード、そして同時に random pair split の限界を一つのケースで検討できるため読む価値があります。論文だけでなく公開コードの split/masking implementation まで読むのが重要です。citeturn43view0  
公開コード: `https://github.com/dmis-lab/ATTNSOM` fileciteturn1file0L2-L6

**SyGMa** は modern ML とは逆方向ですが、metabolite prediction において reaction rules がなぜ今も有用なのか理解するために残すべき古典です。end-to-end model と rule-based generation を比較する際の baseline 以上の意味があります。citeturn39search1  
公開実装: `https://github.com/3D-e-Chem/sygma` fileciteturn2file0L2-L6

### 今後、特に検証価値が高い問い

最も優先度が高いのは、**「cross-isoform auxiliary learning の gain は、未知分子にも残るか」**です。ATTNSOM 型の gain を molecule-instance random split だけでなく、canonical molecule-group、scaffold、temporal split で再測定します。1A2/2C9/2D6/3A4 のみで構わないので、同一 molecule が train/test を跨がない split をまず作り、2C9↔2C19、3A4↔他 CYP の transfer gain を個別に見ます。これは現在の multitask/CheMeleon 系にもそのまま適用でき、モデル backbone を変更せずに答えられる重要な問いです。citeturn43view0 fileciteturn7file0L2-L2

次は、**「単一濃度値を auxiliary target ではなく dose-response の censoring constraint として使うと改善するか」**です。既に single-concentration auxiliary learning を試しているからこそ比較価値があります。exact IC50、upper/lower censored IC50、single-point fraction inhibition を一つの likelihood に統合し、従来の auxiliary head と同じ molecule split で比較します。大規模 qHTS はこの observation model の calibration source として使えます。citeturn30search0

三つ目は、**「structure information は全化合物に入れるより disagreement cases にだけ使った方が情報効率が高いか」**です。まず 2D/representation model と local-reactivity model を走らせ、uncertainty/high disagreement 化合物のみ CYP-specific docking/ensemble MD を実行します。特に 2D6/2C9 では anchor occupancy、3A4 では productive-pose fraction を特徴量にし、平均性能だけでなく disagreement subset の error reduction を測ります。CYP2D6/2C9 の mutagenesis evidence と CYP3A4 plasticity から、この isoform-dependent strategy にはかなり強い機構的根拠があります。citeturn40view5turn31search4turn31search6

四つ目は、**「intrinsic reactivity と accessibility のどちらが失敗原因かをラベル化できるか」**です。SoM ごとに QM/reactivity score と ensemble proximity score を持たせ、

\[
\text{reactive + accessible},
\quad
\text{reactive + inaccessible},
\quad
\text{unreactive + accessible}
\]

に分類します。wet MetID でどの regime が false positive/negative を作るかを測れば、単一 SoM accuracy よりはるかに mechanistic な error analysis になります。Feng/Friesner と SMARTCyp 系の核心を、現代データで直接検証する問いです。citeturn18academia38turn43view0

五つ目は、**「human hepatocyte exact SoM を別 domain として joint modeling できるか」**です。recombinant CYP-specific labels を source tasks、AstraZeneca hepatocyte annotations を target task として、shared intrinsic-reactivity encoder + isoform-specific heads + hepatocyte aggregate head を構築します。CYP-specific model の単純外部評価ではなく、

\[
P(\text{hepatocyte SoM})
=
F(P_{\rm 1A2},P_{\rm 2C9},P_{\rm 2D6},
P_{\rm 3A4},\ldots,\text{other metabolism})
\]

という latent composition として扱う方が面白いです。120 compounds と小さいからこそ、CYP data の biological prior が役立つかを試せます。citeturn42search1turn42search3turn43view0

六つ目は、**「TDI を binary task から \(K_I\)・\(k_{\rm inact}\)・bioactivation task に分解すると何が見えるか」**です。まず reversible potency、TDI IC50 shift、\(K_I\)、\(k_{\rm inact}\)、reactive-metabolite evidence を別列として保持します。shared representation は使っても output heads は分け、reaction alert/QM descriptors が \(k_{\rm inact}\) に、binding descriptors が \(K_I\) に preferentially 効くかを調べます。この hypothesis が成立すれば単なる TDI classifier より説明性が大きく上がります。progress-curve model との組み合わせも有望です。citeturn34search4turn34search7

七つ目は、**「2D6/2C9 の既知 electrostatic mechanism を model explanation が再発見できるか」**です。これは explanation benchmark として非常に強いです。Glu216/Asp301、Arg108 周辺を知らない状態で model を学習させ、matched molecular pairs の protonation/acidic functionality perturbation に対する predicted counterfactual を測ります。その後 wet mutant data と比較します。「attention が residue を指した」ではなく、「residue mutation × ligand chemistry の interaction effect を正しい方向に予測したか」を評価します。citeturn40view5turn31search4

八つ目は、**「3A4 の error は chemistry より state multiplicity で説明できるか」**です。3A4 だけ他 isoform より model disagreement/error が大きい化合物を集め、ligand flexibility、pocket occupancy、multiple-pose entropy、Hill/cooperativity behavior との関係を調べます。もし correlation が強ければ、3A4 だけ ensemble/kinetic model を追加する理由が得られます。CYP3A4 の structure と multisubstrate kinetics はこの hypothesis を強く支持する出発点です。citeturn31search6turn33search34turn33search3

最後に最も重要なのは、**prospective experiment を accuracy competition ではなく mechanism-discrimination experiment として設計すること**です。既存の 2D/tree/CheMeleon/TabPFN 系、reactivity model、docking/MD が全て同意する化合物を多数測るより、これらが違う理由で食い違う 20–50 compounds を事前に選び、MetID、IC50、time course、必要なら mutant/matched analogue を測る方が、CYP 全体への理解は速く進みます。CYP2D6、2C9、3A4 の wet–dry 成功例に共通しているのも、「prediction を当てる」より **perturbation で competing mechanistic hypotheses を潰したこと**だからです。citeturn40view5turn31search11turn31search2

総合すると、次の研究段階で最も期待値が高いのは新しい汎用 molecular encoder ではなく、**「intrinsic reactivity」「productive binding/access」「assay observation」「time-dependent bioactivation」を分離し、CYP と実験条件を明示変数にした multi-view / hierarchical modeling**です。SoM ではその分解が既にかなり成功しており、CYP2D6/2C9 では residue-level wet validation がその物理性を裏づけ、CYP3A4 は multiple binding と kinetics が単純モデルの限界を示しています。TDI と hepatocyte metabolism はまだこの分解が十分に進んでおらず、したがって今後もっとも「理解」が増える余地の大きい領域です。citeturn43view0turn40view5turn31search11turn33search3turn34search7