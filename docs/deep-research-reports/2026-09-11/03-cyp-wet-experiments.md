# CYP阻害・代謝の wet 実験を「機構の切り分け」として読む：何を測り、どこまで言えるのか

CYP の wet assay は、一見すると「基質の代謝物が何 % 減ったか」を測る単純な酵素試験に見えます。しかし、研究上の本質はそうではありません。**観測できるのは、ある実験系・ある基質・ある時間窓での生成物形成または基質消失であり、「阻害機構」そのものを直接観測しているわけではない**からです。IC₅₀、Kᵢ、k_inact/K_I、代謝物プロファイルなどは、そこから特定の反応モデルを仮定して推定した量です。Fowler & Zhang のレビューが強調している通り、CYP 阻害データを適切に利用するには、それぞれの実験系の長所と弱点を理解することが前提になります。citeturn41view1turn39search18

したがって、以下では「どう assay を回すか」よりも、

**観測 → 競合する説明 → それを分ける追加実験 → 最終的に言える範囲**

という研究者の思考過程を中心に整理します。OpenADMET CYP Challenge で既知とのことなので、単純な recombinant CYP fluorescence/MS screen や IC₅₀-shift の説明そのものではなく、**なぜ同じ化合物が別 assay では違って見えるのか、TDI をどこまで MBI と呼んでよいか、CYP3A4 の「阻害能」という概念自体がなぜ一意でないのか、代謝物や細胞系を入れると結論がどう反転し得るか**を重点的に扱います。現在の規制的枠組みでも、プローブや阻害剤の選択性・potency は「その実験条件で」確認することが重視されており、文献 IC₅₀ を普遍定数として移植する考え方とは距離があります。citeturn4search22turn4search1

## CYP assay で本当に観測しているもの

最も基本的な CYP 活性測定は、

\[
S \xrightarrow[\mathrm{CYP,\ NADPH}]{ } P
\]

について、一定時間後の **P の生成量**、あるいは **S の消失量**を測るものです。しかし、生成物形成速度が低下したという観測だけからは、「化合物が CYP の active site に結合した」とはまだ言えません。酵素自体の失活、基質の枯渇、阻害剤の枯渇、非特異的結合、反応生成物の二次代謝、分析検出への干渉などでも同じ観測になります。実際、実用的な LC-MS CYP assay の method validation では、metabolite formation が時間と microsomal protein に対して線形であるかを先に確認した上で K_M、IC₅₀、TDI を評価しています。これは単なる QC ではなく、**「終点の生成物量を初速度の代理として扱ってよい」という推論を成立させる条件**です。メーカー技術資料ですが、BD Gentest の assay validation はこの考え方を具体的に示しています。citeturn41view11turn6view0turn6view1

IC₅₀ は特に誤解されやすい量です。単純な competitive inhibition と Michaelis–Menten kinetics が成立する場合に限れば、おおよそ

\[
IC_{50}\approx K_i(1+[S]/K_M)
\]

という関係があるため、**基質濃度を変えるだけで IC₅₀ は変わります**。したがって「compound X の CYP3A4 IC₅₀ は 1 µM」という文は不完全で、本来は「midazolam を何 µM、どの HLM concentration、どの incubation time で使ったか」まで含めて一つの測定値です。さらに CYP3A4 のように複数リガンド結合・allostery・substrate inhibition が起こり得る CYP では、この単純式そのものが成立しない場合があります。複数基質・複数濃度を使った global nonlinear fitting が、Lineweaver–Burk plot のような単純化より重要になる理由です。CYP kinetics のレビュー・方法論でも、複雑な反応機構を無理に単一 Michaelis–Menten model に押し込まないことが強調されています。citeturn26search18turn26search10turn36view2

TDI ではさらに一段階モデル依存になります。複数の inhibitor concentration について preincubation time と residual activity の関係を測り、各濃度で apparent first-order inactivation rate \(k_{\mathrm{obs}}\) を求め、

\[
k_{\mathrm{obs}}=
\frac{k_{\mathrm{inact}}[I]}
{K_I+[I]}
\]

のようなモデルに fitting して \(k_{\mathrm{inact}}\) と \(K_I\) を得ます。しかしこれは、活性酵素量が時間以外の理由で変化しない、inhibitor concentration が概ね既知である、といった仮定を含みます。Stresser らの詳細な TDI review とその後の数理研究は、**preincubation 中の inhibitor depletion が大きいと IC₅₀ shift が予想より小さくなり、場合によっては preincubation 後に apparent IC₅₀ が逆方向へ動くことさえある**ため、nominal [I] ではなく実測濃度や depletion を含むモデルが必要になり得ると論じています。citeturn5view3turn37search26turn26search19

代謝物同定では別の「観測と結論の距離」があります。LC-HRMS の exact mass、isotope pattern、MS/MS fragmentation、retention time から「+16 Da の酸化体」「脱アルキル化体らしい」と推定することはできますが、**同じ elemental composition を持つ regioisomer や stereoisomer の完全同定は別問題**です。HRMS の構造解析レビューも、高分解能・高 mass accuracy が非常に強力である一方、厳密な structural identification には複数の証拠が必要だとしています。authentic standard との retention time/MS² 一致、必要なら NMR、合成標品や isotope labeling が一段強い証拠になります。citeturn39search23

この「何が言えて、何がまだ言えないか」をまとめると、次のようになります。

| 観測 | かなり強く言えること | その観測だけでは言えないこと |
|---|---|---|
| probe metabolite が低下 | その条件で reaction velocity が低下した | CYPへの直接結合、阻害様式、in vivo DDI |
| IC₅₀ concentration-response | その assay 条件での apparent potency | intrinsic Kᵢ、別 probe での potency |
| 複数 [S],[I] の initial-rate kinetics | kinetic model 間の識別が可能 | CYP3A4 で唯一の普遍的 mechanism |
| NADPH preincubation で IC₅₀ left shift | metabolism/time-dependent component がある | 不可逆的共有結合、MBI |
| 希釈後も活性低下が残る | rapidly reversible inhibition だけでは説明しにくい | covalent modification の確証 |
| dialysis で活性回復 | reversible / quasi-irreversible component を支持 | どの化学種が結合していたか |
| CYP peptide に drug-derived mass adduct | covalent protein modification の強い証拠 | それだけで臨床 DDI magnitude |
| LC-HRMS に新 metabolite peak | drug-derived transformation の候補 | 完全な位置・立体構造 |

特に重要なのは、**assay の出力値より、その値を成立させているモデルと control を読む**ことです。citeturn41view1turn5view3

## 組換え CYP、ミクロソーム、肝細胞は「精度の段階」ではなく、違う問いに答える

recombinant CYP → HLM → hepatocyte を「単純系から高級な系への一本道」と考えるより、**それぞれが別の交絡を消し、別の交絡を導入する実験**と見る方が有用です。Parkinson らが recombinant enzymes、human liver microsomes、human hepatocytes を比較した研究を “system-dependent outcomes” と呼んだのは象徴的で、特に parent compound が CYP 以外でも強く代謝される場合、どの系を使うかによって CYP TDI の見え方自体が変わり、hepatocyte を併用する価値があると結論しています。citeturn39search1turn5view3

**recombinant CYP** の最大の強みは isoform attribution です。「CYP2D6 だけ存在する状況でこの反応が起こるか」という sufficiency test ができる。一方、実際の肝小胞体とは CYP:NADPH-P450 reductase:cytochrome b₅ の比、脂質環境、他 CYP との競合が異なり得ます。そのため recombinant CYP で得られた turnover や IC₅₀ を、そのまま HLM の absolute activity と同一視することはできません。RAF/ISEF のような scaling が使われるのもそのためです。TDI の assay system dependence についても、recombinant、HLM、hepatocyte 間で outcome が異なることがレビューされています。citeturn41view1turn5view3turn39search1

**HLM** は CYP が native microsomal membrane に存在し、複数 CYP、P450 reductase、b₅、UGT などが同居する点で recombinant 系より肝臓に近く、CYP DDI の標準的 in vitro system と長く位置づけられてきました。citeturn39search18turn41view1 ただし、「複数 CYP がある」ということはそのまま解釈上の曖昧さでもあります。例えば HLM で “CYP3A4 activity” とラベルされた midazolam 1′-hydroxylation は実務上有用でも、recombinant CYP3A4 assay のように CYP3A4 だけを直接観測しているわけではありません。CYP3A5 を含む CYP3A family contribution、他経路、donor/lot differences が入り得るためです。逆に、その複雑さを利用して selective chemical inhibitor、抗体、genotyped donor、recombinant CYP を組み合わせることで reaction phenotyping が強くなります。FDA も probe/inhibitor の selectivity と potency を実験条件ごとに検証することを求めています。citeturn4search22

**primary human hepatocyte** になると、膜透過、細胞内結合、phase I/II metabolism、より生理的な cofactor environment などが同時に入ります。TDI について Stresser らは hepatocyte に利点がある証拠が増えているとレビューしており、McGinnity、Zhao、Chen らの hepatocyte TDI studies、さらに plasma 中で hepatocyte を使う試みまで整理しています。citeturn5view3 ここで重要なのは、「HLM より hepatocyte の方がいつも正しい」のではなく、**HLM の nominal inhibitor concentration と hepatocyte の intracellular inhibitor exposure は違う量**だということです。hepatocyte で弱く見えれば permeability が低いのか、細胞内で速く消えるのか、強く結合するのかもしれず、逆に強く見えれば細胞内で inhibitory metabolite が蓄積した可能性もあります。system を変えた差そのものが mechanistic information になります。citeturn39search1turn5view3

非常に教訓的なのが **gemfibrozil–CYP2C8** です。これはユーザーが主対象としている四 CYP からは外れますが、「CYP assay を NADPH-only recombinant system に閉じ込めてはいけない」ことを示す代表例です。Ogilvie らは gemfibrozil の glucuronidation によって生じる **gemfibrozil 1-O-β-glucuronide が強い metabolism-dependent CYP2C8 inhibitor になる**ことを示しました。つまり「parent compound + CYP2C8 + NADPH」だけを調べると、臨床的に重要な perpetrator species を十分再現できない可能性があります。UGT-mediated activation を許す system、生成 metabolite の直接 testing、あるいは hepatocyte が必要になるタイプです。citeturn31search2turn31search3

これは一般化すると、次のような設計になります。

> **recombinant で現象がない → HLM/hepatocyte で現象がある**なら、「recombinant assay が間違った」と考える前に、別 CYP、UGT/その他酵素による inhibitory metabolite formation、cellular accumulation などを候補にする。  
> **recombinant では強い → hepatocyte では弱い**なら、permeability、binding、competing clearance、実際の intracellular unbound exposure を疑う。

つまり system comparison は assay validation ではなく、**metabolic network を perturb する実験**として利用できます。citeturn39search1turn5view3

なお S9 fraction は microsome と cytosol の両方を含むため、CYP のみならず cytosolic metabolism を保持した中間的な選択肢です。反対に intact hepatocyte は cellular disposition を加えます。この「含まれる pathway を一つずつ増減する」という設計思想の方が、単に “more physiological” と並べるより機構研究には有用です。citeturn35search0

## 研究者が機構を切り分けた具体例

CYP 研究で最も重要な教訓の一つは、**同じ inhibitor に一つの「CYP3A4 inhibition potency」が存在するとは限らない**ことです。Racha らは human liver microsomes を用い、midazolam、testosterone、nifedipine、terfenadine という異なる CYP3A4 probe reactions に対する多数の薬物の inhibition profile を比較し、substrate によって profile が変わることを示しました。citeturn18search16 Wang らの CYP3A4 interaction studies や Galetin らの multisite kinetic analysis は、CYP3A4 active site では複数 ligand が同時に存在し得るという考えで、midazolam/testosterone/nifedipine 間の複雑な kinetics を説明しました。citeturn18search20turn36view2 後には active site 内の spatially distinct binding domains を支持する biochemical evidence が得られ、近年も「probe substrate dependence of apparent Kᵢ」は CYP3A4 allosteric inhibition の中心問題として研究されています。citeturn18search1turn18search23turn18search25

ここで研究者がしたことは、「midazolam assay の precision を上げた」ことではありません。

**観測:** inhibitor X の potency が midazolam と testosterone で大きく違う。  
**代替説明:** analytical artifact、probe metabolism の非線形性、単純な competitive inhibition の [S]/K_M 差。  
**切り分け:** 複数 probe、複数 substrate concentration、異なる inhibitor、multisite kinetic model、構造・binding evidence を組み合わせる。  
**結論:** CYP3A4 では inhibitor–enzyme interaction を substrate-independent scalar とみなせないケースがある。  
**まだ言えない:** 一つの in vitro probe から、すべての clinical CYP3A4 substrates に対する DDI magnitude を一意に決めること。

これは OpenADMET のような benchmark を読む際にも重要で、「label noise」の一部には assay error ではなく、**biological assay-context dependence** が含まれ得ます。citeturn18search16turn18search23

TDI ではさらに明瞭な例があります。市販 CRO の IC₅₀-shift assay は、0 min、30 min −NADPH、30 min +NADPH を比較し、+NADPH preincubation で left shift すれば TDI と判定します。Cyprotex/Evotec 自身も TDI と MBI は科学的には異なり、MBI は reactive metabolite による enzyme inactivation という TDI の一 subset だと説明しています。citeturn41view9 一方で同ページは IC₅₀-shift が reversible と irreversible inhibition を区別できるとも表現しており、ここは**技術サービスの説明と、厳密な mechanistic proof を区別すべき場所**です。Fowler & Zhang のレビューでは quasi-irreversible inhibition では dialysis によって CYP activity が回復し得ることが整理されており、TDI という kinetic phenotype と irreversible covalent inactivation は同義ではありません。citeturn37search0turn41view1 Helsinki 大学の学位研究でも、IC₅₀-shift で見つかった TDI の reversibility を評価するために dialysis 法を開発すること自体が研究課題になっています。citeturn39search32turn37search3

したがって mechanistic ladder は、

**NADPH-dependent IC₅₀ shift → dilution → dialysis/recovery → spectral/adduct analysis**

と進めるのが理にかなっています。

preincubation 後に大希釈しても activity loss が残れば、すぐに外れる parent inhibitor だけでは説明しにくくなります。しかし stable tight-binding metabolite や metabolic-intermediate complex はまだ候補です。dialysis で戻れば quasi-irreversible/tight-binding の証拠になり、戻らなければ irreversible mechanism の可能性は上がりますが、それでも **「共有結合ができた」ことの直接証明ではありません**。共有結合を主張するなら drug-derived adduct の intact-protein MS、peptide mapping、heme-adduct analysis などを加えるべきです。citeturn41view1turn5view3

この最終段階まで進んだ好例が **AMG 487–CYP3A4** です。Henne、VandenBrink らの研究では、AMG 487 の sequential metabolism により quinone reactive metabolite が形成され、CYP3A4 の **Cys239 を covalently modify** して TDI を生じることが示されました。citeturn40search0turn5view3 この研究が強いのは、「時間依存的に activity が落ちた」だけで終わらず、**metabolite chemistry と protein modification を接続した**点です。

ここでも、

**観測:** NADPH/time-dependent loss of CYP3A4 activity。  
**代替説明:** parent の slow-binding、単なる stable inhibitory metabolite、非特異的 enzyme deterioration。  
**切り分け:** sequential metabolites を追跡し、reactive quinone formation と CYP protein adduct を同定し、modified residue まで局在化。  
**結論:** metabolism-generated reactive species による mechanism-based CYP3A4 inactivation をかなり直接的に支持できる。  
**まだ言えない:** その chemical mechanism だけから臨床 DDI の大きさは決まらない。in vivo exposure、unbound concentration、enzyme turnover 等が別途必要。

Stresser review がこの研究を、sequential metabolism の重要性を示す例として取り上げているのもそのためです。citeturn5view3

CYP2C9 の歴史的な **tienilic acid** も別の意味で興味深い例です。yeast-expressed human CYP2C9 で mechanism-based inactivation が示され、reactive thiophene oxidation product による covalent modification が長く研究されましたが、その後 reactive intermediate の化学構造について thiophene-S-oxide か epoxide かという解釈が修正されていきました。citeturn28search0turn31search1 教訓は、「irreversible phenotype」まで強くても、**どの reactive intermediate がどこを修飾したかという chemical mechanism は別の証拠階層**だということです。

CYP2D6 については、paroxetine が実務的 TDI validation に利用される例があります。Chen らの 6-in-1 HLM cocktail study は、phenacetin/CYP1A2、diclofenac/CYP2C9、dextromethorphan/CYP2D6、testosterone/CYP3A などを LC-MS/MS で同時測定し、direct inhibition と TDI の validation に paroxetine を含む既知 inhibitor を使用しました。single-substrate assay と cocktail で概ね類似した inhibition curve/IC₅₀ shift を得ています。citeturn36view2 重要なのは、これは **「cocktail がこの validation set では単独 assay を再現した」証拠**であって、「任意の化合物で probe–probe interaction がない」ことの証明ではない点です。

CYP1A2 についても同じ考え方が必要です。実用 assay では phenacetin→acetaminophen が広く使われ、BD の LC-MS platform では tacrine を alternate probe として利用しています。citeturn41view11turn6view1 したがって、fluorescent CYP1A2 reporter で得た hit が phenacetin/tacrine LC-MS で再現するかという比較は、単なる「確認試験」ではなく、**人工 reporter chemistry に依存した interaction か、native-like drug substrate にも一般化する interaction か**を問う実験になります。

## 蛍光・発光・LC-MS と probe 条件が結論を作る

蛍光 CYP assay が HTS に強いのは、plate reader で高速・低コストに product formation を読めるからです。しかし fluorescence reporter には、test compound 自身の autofluorescence、quenching、excitation/emission overlap、reporter substrate 特有の binding mode という別の問題が加わります。Promega が P450-Glo を紹介した技術資料も、従来 fluorescence assay の optical background や spectral interference を問題として挙げています。これは有用な実務情報ですが、**メーカーによる自社技術の比較資料であり、独立した universal superiority の証明とは分けて読むべきです**。citeturn23view0turn41view10

発光 assay では CYP が luminogenic precursor を luciferin-like product に変え、その後 luciferase reaction で光を出すため、excitation light が不要という分析上の利点があります。Promega はそのため background が低く、fluorescent compounds/NADPH の spectral overlap を避けやすいとしています。citeturn23view0turn41view10 しかし mechanistic に見ると、これは「CYP reaction → luciferase detection」という **coupled assay** です。したがって test compound が luciferase 自体を阻害すれば CYP inhibition に見え得ます。独立した luciferase研究でも small molecules による firefly/Renilla luciferase inhibition が false-positive reporter effects を生じ得ることが示されています。citeturn39search4

この場合の良い counterexperiment は単純です。CYP reaction を経ずに検出可能な product/luciferin を assay に加え、test compound が detection stage を抑えるかを測ります。そこで signal が低下するなら、「CYP阻害」という説明は少なくとも単独では成立しません。このような **assay-reagent interference を独立に perturb する**という発想は、一般的な Assay Guidance Manual の assay-interference framework とも一致します。citeturn39search0

LC-MS/MS は native-like probe drug と特定 metabolite を直接分離・定量できるため、この種の optical artifact からはかなり自由です。BD/Gentest の例では、CYP1A2 phenacetin→acetaminophen、CYP2C9 diclofenac→4′-hydroxydiclofenac、CYP2D6 dextromethorphan→dextrorphan、CYP3A4 midazolam→1′-hydroxymidazolam と testosterone→6β-hydroxytestosterone を stable-isotope internal standards とともに定量し、selectivity、calibration、precision、carryover、stability、LLOQ を validation しています。citeturn6view1turn41view11

ただし LC-MS は「artifact-free」ではありません。ion suppression、co-elution、in-source fragmentation、metabolite standard の純度・response factor、chromatographic separation の不足など、artifact の場所が **optics から analytical chemistry へ移る**だけです。実際 BD 資料でも tacrine CYP1A2 assay で in-source fragmentation を method-specific issue として扱っています。citeturn41view11

読み出し法を比較すると、したがって「MS が正しい、fluorescence が間違い」という二分法ではなく、

> **異なる assay chemistry で同じ mechanistic hypothesis が生き残るか**

を問う方が強い設計です。

fluorescent reporter で hit → luminogenic reporter でも hit → drug-probe LC-MS でも hit、という一致は強いですが、三つとも同じ substrate-binding behavior を測っているとは限りません。逆に fluorescence だけ hit なら optical interference または fluorescent probe-specific binding、LC-MS の midazolam だけ hit で testosterone では弱ければ CYP3A4 probe dependence を疑う、と原因別に枝分かれできます。citeturn18search16turn18search23turn39search4

基質濃度も非常に重要です。screening assay を [S] ≈ K_M 付近で回すのは便利ですが、competitive inhibitor の apparent IC₅₀ は [S]/K_M に依存します。メーカーの BD assay でも substrate concentration、HLM protein、incubation time を isoform ごとに最適化し、linearity と K_M を先に評価しています。citeturn41view11turn6view0 よって異なるデータベースや Challenge data を比較するとき、「同じ CYP」の列を比較するだけでは不十分で、**probe identity、[S]/K_M、protein concentration、incubation time、preincubation protocol**まで assay metadata として扱うべきです。

cocktail assay も同様です。Chen らは 0.1 mg/mL HLM 中で phenacetin、coumarin、bupropion、diclofenac、dextromethorphan、testosterone の濃度と incubation time を再最適化し、single probe と cocktail の IC₅₀ が validation compounds で概ね一致することを示しました。citeturn36view2 これは multiplexing の実証として価値がありますが、cocktail の「非干渉性」は assay の永久的性質ではなく、**使用した probe concentrations と test set について検証された経験的性質**と読むべきです。

## 見かけの阻害を作る落とし穴をどう切り分けるか

最も見落としやすいのは **nominal concentration と enzyme が実際に見る free concentration の違い**です。hydrophobic compound が microsomal membrane/protein に強く非特異的結合すれば、添加した 10 µM のうち CYP にアクセス可能な濃度はずっと低いかもしれません。Foti & Fisher の bufuralol の研究などを引用しながら、TDI review は non-specific binding と incubation conditions が kinetic interpretation に影響することを扱っています。citeturn5view3 Cyprotex 自身も、preincubation→dilution protocol で microsomal protein concentration を高くすると protein binding が増え、free inhibitor が減少して potency を過小評価し得ることを説明しています。citeturn41view9

これは非常に便利な mechanistic perturbation を与えます。同じ nominal [I] で microsomal protein concentration を変えたとき apparent IC₅₀ が大きく動くなら、単純な intrinsic CYP affinity だけでは説明しにくくなります。equilibrium dialysis 等で \(f_{u,mic}\) を測る、あるいは incubation 中の aqueous/unbound concentration を確認する価値が生じます。実務家向けの microsomal protein-binding 講演でも、MPB を clearance/DDI interpretation と結びつけ、測定法と study design を独立した課題として扱っています。これは practitioner source であり、一次実証研究とは区別すべき資料です。citeturn34search0

**inhibitor depletion** は TDI で特に危険です。preincubation の目的は inhibitor に CYP を「攻撃させる」ことですが、同じ間に inhibitor 自身が代謝されます。Berry らの解析を紹介した Stresser review では、これが apparent IC₅₀ shift を小さくしたり、場合によっては想定と逆向きにしたりし得るとされています。citeturn5view3 したがって不思議な TDI curve を見たとき、すぐに複雑な enzyme mechanism を考える前に、**0、5、15、30 min の actual parent [I] と主要 metabolite [M] を LC-MS で測る**だけで説明が大きく変わることがあります。

**substrate depletion** も同じ問題です。endpoint assay は「平均速度」を測っているので、反応後半で substrate がかなり減れば initial-rate kinetics ではなくなります。metabolite がさらに代謝される場合は product formation 自体が plateau/decline することもあります。したがって time-course の progress curve を一度見ることは、単なる incubation-time optimization 以上の意味があります。substrate depletion を利用して K_M/V_max を推定する方法自体も存在しますが、その場合は depletion kinetics を明示的に model 化する必要があります。citeturn26search10turn26search6

興味深い実例として、肝 microsome/hepatocyte の depletion が一相性でないことを検討した学位研究では、化合物によって **biphasic drug depletion** が観測され、早期時間点の消失、co-eluting metabolites、モデル選択そのものが clearance interpretation に影響しました。citeturn26search21 これは「log(parent) vs time を直線 fit すれば clearance」という routine workflow の背後にも、機構的仮定があることをよく示します。

**enzyme lability** は TDI のもう一つの代替説明です。30分 preincubation 後に activity が低ければ、test compound が inactivator だったのか、CYP が単に incubation 中に失活したのかを分ける必要があります。したがって vehicle control の time dependence、−NADPH/+NADPH control、既知 reversible inhibitor/TDI inhibitor を同時に走らせることが mechanistic control になります。TDI methodology reviews が dilution、enzyme stability、protein binding、inhibitor depletion を一体として論じているのはこのためです。citeturn41view1turn5view3

**溶解性・析出・凝集**は、CYP に限らず high-concentration biochemical screening で必ず考えるべき代替仮説です。ただし CYP は membrane-associated enzyme なので、一般的な soluble-enzyme aggregation countermeasure をそのまま移植するのも危険です。ここで重要なのは「化合物は凝集するはずだ」と推測するのではなく、actual dissolved concentration、遠心前後、orthogonal detection、concentration-response の形などから仮説を試験することです。assay-interference の一般原則として、compound と assay reagent/biological component の非目的反応は本来の target mechanism と切り分けなければなりません。citeturn39search0

さらに、**急峻な concentration-response を mechanistic cooperative inhibition と即断しない**ことも重要です。precipitation、aggregation、free fraction の非線形変化、detector saturation でも Hill slope は変わり得ます。逆に CYP3A4 では本物の multisite/allosteric kinetics もあり得るため、Hill slope だけで artifact を宣告することもできません。ここで CYP3A4 の複数 probe experiment が強力になります。citeturn18search20turn18search23

落とし穴を「何を追加測定すれば分かれるか」で整理すると実用的です。

| 不思議な観測 | 有力な代替説明 | 情報量の大きい追加実験 |
|---|---|---|
| fluorescence だけ強阻害 | autofluorescence / quench | LC-MS drug probe、compound-only optical control |
| luminescence だけ強阻害 | luciferase inhibition | CYPを迂回した detection-stage control |
| protein を増やすと IC₅₀ 上昇 | nonspecific microsomal binding | \(f_{u,mic}\)、実測 free/total compound |
| preincubation後 TDI が弱くなる | inhibitor depletion | time-course parent/metabolite LC-MS |
| −NADPHでも時間依存低下 | enzyme instability、非CYP反応 | vehicle time-course、別 cofactor/system |
| midazolamだけ強く阻害 | CYP3A4 probe dependence | testosterone/nifedipine 等の orthogonal probe |
| recombinantで陰性、hepatocyteで陽性 | inhibitory metabolite/別 pathway | metabolite ID、recombinant metabolite testing |
| HLMで強く hepatocyteで弱い | uptake/free exposure/competing metabolism | intracellular exposure、parent depletion |
| dilution後も activity loss | TDI/tight binding | dialysis、recovery、adduct/spectral analysis |
| exact-mass metaboliteだけ見える | structural isomer ambiguity | authentic standard、MS/MS、必要なら NMR |

この表のポイントは、**artifact control と mechanism experiment は同じもの**だということです。良い control は単に「 assay が正常だった」と示すのではなく、競合する機構を一つ消します。citeturn41view1turn5view3turn39search0

## ハイスループットから精密機構・生理的系へどうつなぐか

HTS と mechanistic enzymology は対立するものではありません。むしろよい workflow は、「全 compound に高価な mechanistic assay をする」のでも「IC₅₀ screen を真実として終える」のでもなく、**最も不確実性を減らす assay に枝分かれさせる**ものです。

例えば BD/Gentest の RapidFire-MS 技術資料では、LC-MS/MS assay について time/protein linearity、K_M、IC₅₀、TDI、analytical performance を validation した上で高速 MS platform に移しています。43 の IC₅₀ curves を conventional LC-MS と RapidFire に split analysis した比較では \(R^2=0.97\)、大きな systematic bias がないと報告されています。citeturn41view11 これは「高速化しても同一組織内の reference workflow を概ね再現できる」という有用なメーカー validation です。ただし、**同じ assay architecture の platform comparison であり、異なる研究室・異なる HLM lot・異なる probe・異なる analysis pipeline 間の universal reproducibility を証明したものではありません**。

Chen らの cocktail study も同様に、適切に再最適化すれば single-probe assay と multiplex LC-MS assay の direct inhibition/TDI が validation compounds で概ね一致することを示しています。citeturn36view2 ここから得るべき教訓は「cocktail は絶対安全」ではなく、**high-throughput化は mechanistic assumptions を validation set に対して再確認してから行う**という設計です。

その上で、研究上は次の「evidence ladder」が有効です。

| 問いたいこと | 最も効率のよい実験 | 次の段階で除ける説明 |
|---|---|---|
| hit は assay chemistry 由来か | optical ↔ LC-MS orthogonal assay | fluorescence/luciferase interference |
| probe-specific か | 同 CYP の別構造 probe | substrate-dependent inhibition |
| isoform-specific か | recombinant panel + HLM phenotyping | 他 CYP contribution |
| parent が inhibitor か | parent と主要 metabolite を個別 assay | metabolite-mediated inhibition |
| metabolism が必要か | ±NADPH / ±必要 cofactor | direct reversible inhibition |
| TDI は回復するか | dilution → dialysis/washout | rapidly reversible vs persistent |
| 共有結合か | intact/peptide/heme adduct MS | tight binding・MI complex |
| 細胞でも成立するか | hepatocyte + intracellular exposure | microsome-only artifact |
| 臨床量として重要か | unbound exposure + IVIVE/PBPK | biochemical possibility ≠ clinical relevance |

これが一次 screen から mechanistic study に進む際の重要な違いです。**次の assay は「より精密な同じ測定」ではなく、現在残っている competing hypothesis を最もよく分離する assay を選ぶ**べきです。TDI literature でも、simple IC₅₀ shift、k_inact/K_I、hepatocyte assay にはそれぞれ misclassification の仕方があるため、一つの assay を絶対基準にしないことが議論されています。citeturn5view3turn41view1

例えば IC₅₀ shift が陽性なら、毎回ただ k_inact/K_I を機械的に測るより、

「preincubation 中に parent が残っているか」  
→「大希釈で activity loss が残るか」  
→「dialysis で戻るか」  
→「reactive metabolite が見えるか」  
→「protein/heme adduct があるか」

と進める方が、**“time-dependent” から “chemical mechanism” へ証拠を積み上げる**ことになります。Cyprotex のような商用 assay の説明でも、IC₅₀ shift 後に K_I/k_inact と機構評価へ進むことが推奨されていますが、IC₅₀ shift 自体を irreversibility の最終証明としないのが重要です。citeturn41view9turn41view1

より生理的な hepatocyte への展開でも、単に「結果が再現したか」では不十分です。HLM と hepatocyte で差が出たら、差を消すのではなく、**parent disappearance、intracellular exposure、metabolite appearance、washout 後の CYP activity recovery**などを時間軸で並べると mechanistic information になります。Parkinson らの system-comparison study と、hepatocyte TDI literature が示す重要点はまさにそこです。citeturn39search1turn5view3

reproducibility についても、「IC₅₀ が何倍以内なら再現」という一指標だけでは弱いです。最低限、

**probe identity / concentration relative to K_M、enzyme source/lot、protein concentration、preincubation/incubation time、cofactor composition、solvent、actual test-compound exposure、positive controls、curve model**

を一緒に持つ必要があります。FDA が inhibitor selectivity/potency を同一 experimental conditions で確認するよう求めているのは、CYP potency が assay context から完全には切り離せないためです。citeturn4search22

この意味では OpenADMET Challenge のデータ解釈でも、同じ CYP label の値がずれるとき、

**measurement noise**  
**assay chemistry difference**  
**probe dependence**  
**system dependence**  
**true compound chemistry**

を別々の latent source と考えるのが重要です。CYP3A4 は特に probe dependence、TDI compounds は preincubation protocol、highly lipophilic compounds は free fraction の影響を受けやすいため、assay metadata が predictive model にとって単なる付帯情報ではなく、**label-generating mechanism の一部**になります。citeturn18search16turn18search23turn5view3

## データを読むときの視点と、さらに調べる価値のある実験上の問い

最終的に CYP wet data を読む際に最も有用なのは、IC₅₀ の大小よりも **「この数字を別の機構でも説明できるか」**と問うことです。

まず、**“CYPを阻害した” と “probe product が減った” を分ける**べきです。後者は観測、前者は仮説です。probe を変える、readout を変える、protein concentration を変える、preincubation を入れる、希釈する、dialysis する、system を recombinant→HLM→hepatocyte と変える。これらはすべて「条件を最適化する操作」ではなく、特定の機構だけを敏感に変化させる perturbation です。citeturn41view1turn5view3

次に、**IC₅₀ は compound property ではなく compound × CYP system × probe × protocol の property** と考えるのが安全です。これは CYP3A4 で最も露骨ですが、他 CYP にも程度の差こそあれ当てはまります。CYP1A2 の phenacetin と alternate probe、CYP2C9 の diclofenac、CYP2D6 の dextromethorphan、CYP3A4 の midazolam/testosterone という assay identity を値と一緒に保存することが、異なるデータセットを統合する前提になります。citeturn41view11turn18search16

TDI データでは、**“NADPH-dependent” → “time-dependent” → “persistent after dilution” → “nonrecoverable” → “covalent adduct”** を一続きの同義語にしないことが最重要です。それぞれ証拠の強度が違います。AMG 487 のように metabolite chemistry と CYP3A4 residue modification まで接続した研究と、単一 IC₅₀-shift positive compound を同じ「MBI」という箱に入れると情報量を失います。citeturn40search0turn41view9turn41view1

代謝物についても、**“何ができたか” と “誰が作ったか” と “その metabolite が何をするか” は別実験**です。HRMS metabolite ID、recombinant CYP formation experiment、HLM chemical/antibody inhibition、metabolite authentic-standard inhibition assay、hepatocyte exposure を重ねて初めて、parent→metabolite→CYP inhibition という causal chain が強くなります。gemfibrozil/CYP2C8 が示すように、phase-II metabolism が CYP inhibition を「消す」のではなく、むしろ強い inhibitor を「作る」こともあります。citeturn39search23turn31search3

今後さらに掘る価値が大きい実験上の問いは、次のようなものだと思います。

| 問い | なぜ重要か | 特に関係する CYP |
|---|---|---|
| 同じ compound の inhibition landscape を 2–3 probe で測るとどれほど変わるか | assay label の substrate dependence を直接定量できる | **CYP3A4** 最重要 |
| fluorescence/luminescence hit の何割が native-probe LC-MS で残るか | HTS-specific artifact と genuine CYP interaction を分離できる | 1A2/2C9/2D6/3A4 |
| nominal IC₅₀ を unbound IC₅₀ に直すと assay間差がどこまで縮むか | microsomal binding を label variability から分離できる | lipophilic compounds 全般 |
| TDI preincubation 中の actual parent/metabolite concentration を同時測定すると k_inact/K_I がどう変わるか | inhibitor depletion bias を直接検証できる | 3A4、2D6 など |
| IC₅₀-shift positive compounds の何割が dilution/dialysis で回復するか | TDI と irreversible MBI の混同を定量できる | 全 CYP |
| HLM 陰性/hepatocyte 陽性 compounds に共通する metabolite pathway は何か | UGT等による inhibitory metabolite generation を発見できる | CYP2C8 の例を他 CYPへ |
| HLM と hepatocyte の差を intracellular unbound concentration で説明できるか | 「system difference」を exposure difference に分解できる | 全 CYP |
| CYP3A4 midazolam/testosterone disagreement は構造特徴から予測できるか | ML label noise を mechanistic signal に変えられる | **CYP3A4** |
| donor/genotype 別 HLM/hepatocyte で phenotype がどう変わるか | pooled material が消している biological heterogeneity を評価できる | **CYP2D6、2C9、3A5/3A4** |
| activity loss と CYP protein abundance を同時測定すると何が分かれるか | catalytic inactivation と protein loss/turnover を分けられる | TDI 全般 |

これらの中でも OpenADMET の観点で特に価値が高いのは、**同一化合物を複数 assay architecture で意図的に再測定する “bridge dataset”** だと思います。例えば recombinant fluorescence、recombinant LC-MS、HLM-midazolam、HLM-testosterone、HLM ± preincubation、hepatocyte を完全に同じ compound panel に対して測れば、「どの assay が正解か」ではなく、**assay間変換可能性と assay-specific biology** を学習できます。CYP3A4 の substrate dependence、TDI の system dependence、microsomal binding の影響が既に独立研究で知られていることを考えると、このようなデータは単なる replicate より情報量が大きいはずです。これは既存研究からの推論です。citeturn18search16turn18search23turn39search1turn5view3

結局、強い CYP mechanistic study には共通した型があります。**一つの readout を精密に測ることよりも、その readout を生む複数の説明を並べ、それぞれに違う予測を与える perturbation を設計する**ことです。CYP3A4 substrate dependence の研究は probe を変え、TDI 研究は時間・NADPH・希釈・回復を変え、AMG 487 の研究は activity assay から metabolite chemistry と protein adduct まで測定階層を上げ、gemfibrozil の研究は CYP assay の外側にある UGT metabolism を取り込みました。citeturn18search16turn40search0turn31search3turn5view3

その意味で、wet CYP data を読む際の最も有用な一文は、**「この実験は何を直接観測したのか。そして、著者が主張する機構以外に、その観測を説明できるものを何個消しているのか」**です。IC₅₀ の小数点以下より、この問いへの答えの方が、研究の mechanistic resolution をはるかによく表します。