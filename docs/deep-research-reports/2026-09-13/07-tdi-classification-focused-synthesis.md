# CYP TDI分類の機序・予測研究・外部データ

## 要約

OpenADMETの分類ラベルは、厳密な共有結合型mechanism-based inhibition（MBI）ではない。これは、組換えCYPを化合物と30分間preincubationしたとき、NADP(H)ありの条件でNADP(H)なしよりIC50が2倍を超えて左へ移動するか、という操作的なtime-dependent inhibition（TDI）ラベルである。したがって、共有結合型MBI、hemeとのmetabolite-intermediate complex（MIC）、強い可逆阻害代謝物の生成を区別しない。IC50-shiftは規制科学でも初期screenとして使われるが、陽性後に `K_I` と `k_inact` を測る精密評価とは異なる。[^1][^2]

この区別はモデル設計に直結する。TDIは単に「反応性官能基があるか」ではなく、少なくとも、対象CYPへの結合、対象CYPによるturnover、適切なsite of metabolism、反応性中間体または強阻害代謝物の生成、活性部位内での捕捉、30分以内という時間条件の積として成立する。よって、通常のMorgan fingerprintだけでも、構造alertだけでも機序の一部しか表現できない。

調査から最初に試す価値が高いのは次の3系統である。

1. **TDI構造alertを明示的に加えたLightGBM**：既存のMorgan + Mordredにbioactivation motif、alert周辺環境、CYP認識特徴との交互作用を追加する。
2. **challenge-native paired/multitask model**：TDI booleanだけでなく、direct pIC50、TDI-arm pIC50、両者の差、Emax、一次screenを、SMILESのみを入力として同時予測する。
3. **assay/source別の外部TDI transfer**：CYP3A4の公開TDI・`kobs` データを同一booleanへ単純結合せず、source-specific auxiliary headとして利用する。

CYP2D6の公開TDIデータはCYP3A4より著しく少ない。D6については外部ラベルの大量追加より、D6への結合・基質性とbioactivationのAND条件を表現することが現実的である。

## TDI、MDI、MBIは同義ではない

TDIは「時間とともに阻害が強くなる」という観察現象であり、機序を限定しない。Metabolism-dependent inhibition（MDI）は、その時間依存性がNADPH依存の代謝を必要とする場合を指す。MBIはさらに狭く、化合物が標的CYPの基質として代謝され、その活性部位内で生じた反応性中間体が同じCYPを不活化する場合である。[^3][^4]

MBIには、反応性中間体がアポタンパク質またはhemeへ共有結合する不可逆的経路と、nitrosoやcarbeneが還元型heme鉄に強く配位する準不可逆的MIC形成がある。[^5] 一方、親化合物のslow bindingや、代謝によって生じた強い可逆阻害剤でも時間依存的な挙動は起こり得る。OpenADMETの±NADP(H)比較は純粋なNADPH非依存slow bindingをある程度相殺するが、共有結合型MBIと強い可逆阻害代謝物を識別する設計ではない。

ICH M12は、30分間の±NADPH preincubation後に1.5倍または2倍以上のIC50左シフトがあれば、酵素不活化の可能性を示すscreening resultとして扱う。陽性後は `kobs` の時間・濃度依存性を測定し、`K_I` と `k_inact` を求めることが推奨される。[^2] したがって、外部のMBI、`kobs`、clinical DDIラベルをOpenADMETの2-fold shiftと同じ教師ラベルへ変換することは避けるべきである。

## 構造から学ぶべき機序

代表的なbioactivation motifには、methylenedioxyphenyl/benzodioxole、primaryまたはsecondary amine、N-dealkylation可能なtertiary amine、cyclopropylamine、terminal alkyne、furan、thiophene、aniline/aminophenol、hydroquinone precursor、thioamide、thiourea、hydrazineがある。生成物はcarbene、nitroso、iminium、radical、epoxide、quinone-imineなどで、heme配位、protein adduct、heme alkylationへつながり得る。[^6][^7]

ただし、alertの存在だけでは不十分である。同じ官能基でも、立体遮蔽、電子効果、対象CYPへの結合姿勢、競合する安全な代謝経路、反応性中間体が活性部位から逃げる速さで結果が変わる。FDA/Leadscope研究は58個のMBI alertを整理したが、外部検証でalert profiler単独は高感度・低特異度だった。[^7]

OpenADMET trainを用いた主催者の後ろ向き解析では、methylenedioxyphenyl、cyclopropylamine、furan、thiophene、aniline/aminophenolなどが実際にTDI陽性を濃縮した。methylenedioxyphenylのenrichmentはCYP3A4で2.16、CYP2D6で2.63だった。一方、CYP3A4陽性の半分以上は既知alertを一つも持たず、alert panelのrecallは48.2%だった。[^8] Alertは追加特徴として有望だが、rule-based classifier単体にはできない。

## CYP2D6とCYP3A4を分けて考える

CYP2D6はCYP3A4より小さく、basic nitrogenとGlu216/Asp301などの酸性残基との相互作用が基質・阻害剤認識に重要である。[^9] 既知のD6 MBIにはparoxetine、MDMA、SCH66712、EMTPPなどがある。しかしbasic amineは通常阻害や基質認識には有用でも、それ単独ではTDIを意味しない。D6では、basic center、芳香環との距離、pKa proxy、bioactivation motif、予測site of metabolismの交互作用を表す必要がある。

CYP3A4は大きく柔軟な活性部位を持ち、複数リガンド結合、複数姿勢、協同性、probe-substrate依存性を示す。Ketoconazoleを2分子同時に保持する構造も報告されている。[^10] 3A4ではlipophilicity、分子サイズ、柔軟性、芳香族表面と局所反応性の組み合わせが重要になりやすい。単一の静的構造へのdock scoreは、この可塑性とassay probe依存性を十分に表現しない可能性が高い。

## IC50-shift labelの限界

30分のscreenは、遅い逐次代謝を必要とするinactivatorを取り逃す可能性がある。AMG 487では親化合物から複数段階の代謝を経て反応性代謝物が生成し、30分ではなく長いpreincubationでshiftが現れ、probe substrateによっても結果が異なった。[^11] またpreincubation中に親化合物が非阻害性代謝物へ枯渇すると、真のTDIでもshiftが小さくなる場合がある。[^12]

組換えCYPはisoform attributionに優れるが、CYP:reductase:cytochrome b5比、膜環境、別CYPやphase II enzymeの不在により、HLMやhepatocyteとは異なる代謝物分布を作り得る。CYP3A4ではprobe substrateによる見え方の違いも大きい。したがって、外部データを利用するときはenzyme source、preincubation time、probe、inhibitor concentration、dilution条件を捨てずに保持する必要がある。

## 既存の計算予測研究

公開されたTDI/MBI予測研究は、通常のCYP阻害予測より少なく、特にCYP2D6は乏しい。

| 研究 | Endpoint | 規模 | 主な注意点 |
|---|---|---:|---|
| Xu et al. | CYP3A4 potential TDI binary | 7,914、581 positive | GraphConv等を比較。CVは非常に高いが、15化合物のwet validationではMCC 0.47まで低下。系列・sampling・random splitへの注意が必要。[^13] |
| Fluetsch et al. / Novartis | CYP3A4 `log kobs` | proprietary実測学習＋368薬test | D-MPNNの主タスクを`log kobs`、3A4/2C9/2D6 reversible pIC50を補助タスクにした。scaffold-balanced splitとprospective testを使用。[^14] |
| Faramarzi et al. / FDA | CYP3A4 heterogeneous TDI binary | 620、約49% positive | 1.5-fold shift、`kobs`、percent loss等を統合。外部性能はcoverage条件付き。58 MBI alertを整理。[^7] |

Xuデータは公開されており、外部TDI候補として大きい。しかし非常に高いrandom-CV性能をそのままOpenADMETのUMAP splitへ期待できない。[^13] Novartis研究は、TDI関連主タスクと通常阻害補助タスクを同じgraph encoderで学ぶ設計として最も参考になる。ただし公開されている約16,000化合物の値の大部分は、proprietary modelが付与したsurrogate predictionであり、実測ラベルとして扱ってはいけない。[^14]

## 外部データ候補

| 優先度 | 資源 | 利用法 | 主な制約 |
|---|---|---|---|
| 高 | FDA/Frontiers CYP3A4 TDI supplement、620化合物 | raw IC50 shiftを確認し、2-foldで再ラベルできる行は近い補助教師、それ以外はsource-specific head | Endpoint混在、承認薬重複、Leadscope本体はproprietary[^7] |
| 高 | ChEMBL TDI関連assays | D6を優先し、±NADPH、30分、paired IC50、`kobs`を条件付きrecordとして抽出 | 3A4約845、D6約154 unique moleculeという検索下限。assay heterogeneityと文献系列重複が大きい[^15] |
| 中 | Zimmerlin 400 reference drugs | CYP3A `kobs` continuous auxiliaryまたは外部検証 | IC50 shiftではなくsingle-point `kobs`。承認薬/test重複リスクが高い[^16] |
| 中 | Xu CYP3A4 TDI dataset | source-specific binary head、外部pretraining | Negative provenance、系列重複、OpenADMETとのラベル差を監査[^13] |
| 中 | Novartis public surrogate | 低weight distillation、teacher feature | 大部分が実測ではない。評価データから完全除外[^14] |
| 補助 | Veith/TDC CYP inhibition | D6/3A4 binding・通常阻害の補助タスク | TDIではない[^17] |
| 補助 | OpenADMET CYP inhibition/reactivity release | TDI-arm potencyまたはsubstrate propensity | paired direct armがなくTDI labelを作れない。同一生成元・重複監査が必要[^18] |

外部データは次の4階層で保持するのが安全である。

- A：±NADPH、約30分、paired IC50から2-fold shiftを再計算可能
- B：IC50 shiftだが閾値、時間、probe、enzyme sourceが異なる
- C：`kobs`、`k_inact/K_I`、remaining activity
- D：通常阻害、substrate、SoM、clinical narrative

A以外をOpenADMETと同じbooleanへ潰さず、別head、別loss、sample weightまたはpretrainingへ分ける。

## 計算可能な機序特徴量

最初に使うべきなのは、RDKitでversion管理したMBI SMARTSである。正確なLeadscope feature定義はproprietaryなので、FDA論文の図や公開OpenADMET実装を基に独自SMARTSとして実装し、同一物と称さず、train内enrichmentとCV増分を検証する。

次にSMARTCypが候補になる。SMARTCypは局所反応性とaccessibilityからCYP site of metabolismを順位づけし、CYP2D6/3A4向けscoreを出せる。[^19] Alert原子とtop-ranked SoMが一致するか、top-1/top-3 score、代謝候補原子数を特徴にすれば、「reactive motifが存在する」から「対象CYPがその位置を酸化し得る」へ一段近づく。

さらにBioTransformer/CyProductやSyGMaで候補代謝物を生成し、親から代謝物へのfingerprint・物性差、electrophile alert、heme-ligating motif、通常阻害scoreの変化を集約できる。[^20][^21] ただし代謝物生成誤差が上流に入るため、alert + SMARTCypで増分が出た後の実験とする。DockingはCYP3A4の可塑性を考えるとさらに低優先で、全件特徴量化よりerror analysis向けである。

## 現在のリポジトリとの照合

すでにMorgan/Mordred LightGBM、CheMeleon系表現、TabPFN/TabICL、challenge TDI headのChemprop fine-tuning、一般CYP補助学習、固定平均・Caruana ensembleは試されている。一方、次は未検証または十分に分離検証されていない。

- Version管理されたTDI bioactivation SMARTSの明示的追加
- Alert周辺の局所fingerprintとCYP認識特徴とのinteraction
- SMARTCypなどのisoform-specific SoM score
- Direct pIC50、TDI-arm pIC50、delta、booleanを同じモデルで学ぶpaired multitask
- 外部の実測TDI/`kobs`をsource-specific headとして使うtransfer
- 予測代謝物を明示的に表現するparent→metabolite model

既存のChemprop auxiliary実験は、一般回帰補助headを足した比較であり、Novartis型のchallenge-native paired potency/shift設計や、assay-aware外部TDI transferを十分に反証していない。

## 推奨実験順

### P0：データを増やさない機序特徴量ablation

既存のMorgan + Mordred LightGBMへ、12–20個のTDI alert、alert数、alert周辺Morgan、2D6 basic-center/pKa proxy、3A4 lipophilicity/size/flexibilityとのinteractionを追加する。D6/A4を別々に評価し、既存UMAP foldとhonest inner thresholdを固定する。

成功条件は、単なる全OOF fitではなくouter MCC/AUROCが一貫して改善し、少数alert系列だけに依存しないこと。失敗条件は、alert陽性群だけをmemorizeし、alertなし陽性のrecallまたはfold間安定性を悪化させることである。

### P1：challenge-native paired multitask

共有encoderとendpoint-specific headで、D6/A4 boolean、direct pIC50、TDI-arm pIC50、delta、Emax、single-concentration responseをmasked multi-task学習する。入力はSMILESのみとし、観測値をinference featureにはしない。各outer foldで全補助headからheld-out compoundを除外し、early stoppingとMCC thresholdはinner foldだけで決める。

これはNovartisの `log kobs + reversible pIC50` multi-task設計と、前回Deep Researchが提案したpaired-potency modelをchallenge固有データで検証する実験である。

### P2：外部TDI inventoryと重複監査

FDA 620、Xu 7,914、Zimmerlin/Novartis observed subset、ChEMBLのD6/3A4 TDI assaysを取得する。各rowにsource、assay type、enzyme source、probe、preincubation、concentration、endpoint、observed/surrogateを保持する。Challenge train/testとstandardized InChIKey、Murcko scaffold、Morgan最大類似度、source DOI/patentで重複監査する。

Test exact matchの外部ラベルは規約確認まで隔離する。External evaluationにも学習にも同じ系列・文献由来化合物を跨がせない。

### P3：assay-aware external transfer

OpenADMET label、FDA/Xu binary、Zimmerlin `log kobs`、ChEMBL shift/kineticsを別headとして学習し、shared encoderだけを共有する。`no external`、`experimental external only`、`+ surrogate distillation`を分ける。CYP3A4から始め、D6ではA4 TDI representation transferとD6-specific headを比較する。

### P4：SoM・代謝物特徴

P0でalertに増分があればSMARTCyp特徴を追加する。さらに増分があればtop-k予測代謝物のreactive motifとparent→metabolite descriptor deltaを試す。単一静的dockingは最後に、D6の限定診断または特定系列でのみ評価する。

## 判断

最も安価で根拠が直接的なのは、**既存LightGBMへのTDI alert + isoform-recognition interaction追加**である。次に重要なのは、**booleanだけを予測する現在の分類問題を、structure → direct potency / TDI-arm potency / shift → booleanという多目的問題へ戻すこと**である。

外部データはCYP3A4では有望だが、label定義が揃っていない。CYP2D6は公開実測TDIが少ないため、外部ラベルの量よりも、対象CYPによる代謝可能性と反応性部位の組み合わせを表現する方が優先度が高い。

## Sources

[^1]: OpenADMET. “[Lowering Inhibitions, One CYP at a Time: Building the OpenADMET CYP Inhibition and TDI Dataset](https://openadmet.github.io/octant-cyp-inhib-blog-post/).” 2026.
[^2]: ICH. “[M12 Drug Interaction Studies](https://www.fda.gov/media/161199/download).” Final guideline.
[^3]: Riley RJ, Grime K, Weaver R. “[Time-dependent CYP inhibition](https://pubmed.ncbi.nlm.nih.gov/17269894/).” 2007.
[^4]: Hollenberg PF, Kent UM, Bumpus NN. “[Mechanism-based inactivation of human cytochromes P450s](https://pubmed.ncbi.nlm.nih.gov/18052110/).” 2008.
[^5]: Isoherranen N et al. “[The nitroso–heme complex of CYP inhibition](https://pmc.ncbi.nlm.nih.gov/articles/PMC2898504/).” 2010.
[^6]: Fontana E, Dansette PM, Poli SM. “[Cytochrome P450 enzymes mechanism based inhibitors: common sub-structures and reactivity](https://pubmed.ncbi.nlm.nih.gov/16248836/).” 2005.
[^7]: Faramarzi S et al. “[Novel (Q)SAR models for prediction of reversible and time-dependent inhibition of cytochrome P450 enzymes](https://pmc.ncbi.nlm.nih.gov/articles/PMC11860084/).” Frontiers in Pharmacology, 2025; [supplementary dataset](https://pmc.ncbi.nlm.nih.gov/articles/instance/11860084/bin/DataSheet1.zip).
[^8]: Walters P. “[Do Our CYP Structural Alerts Actually Work?](https://openadmet.ghost.io/do-our-cyp-structural-alerts-actually-work/).” OpenADMET, 2026.
[^9]: Rowland P et al. “[Crystal structure of human cytochrome P450 2D6](https://pubmed.ncbi.nlm.nih.gov/16352597/).” Journal of Biological Chemistry, 2006.
[^10]: Ekroos M, Sjögren T. “[Structural basis for ligand promiscuity in cytochrome P450 3A4](https://pmc.ncbi.nlm.nih.gov/articles/PMC1564212/).” PNAS, 2006.
[^11]: PubChem. “[AID 1218817: AMG 487 CYP3A4 time-dependent inhibition](https://pubchem.ncbi.nlm.nih.gov/bioassay/1218817).”
[^12]: Berry LM et al. “[Impact of inhibitor depletion on IC50-shift assessment](https://pubmed.ncbi.nlm.nih.gov/23649703/).” 2013.
[^13]: Xu Y et al. “[Development of In Silico Models for Predicting Potential Time-Dependent Inhibitors of CYP3A4](https://doi.org/10.1021/acs.molpharmaceut.2c00571).” Molecular Pharmaceutics, 2023; [data and supplementary information](https://acs.figshare.com/articles/dataset/Development_of_In_Silico_Models_for_Predicting_Potential_Time-Dependent_Inhibitors_of_Cytochrome_P450_3A4/21666785).
[^14]: Fluetsch A et al. “[Using Machine Learning to Predict Time-Dependent Inhibition of CYP3A4](https://doi.org/10.1021/acs.chemrestox.3c00305).” Chemical Research in Toxicology, 2024; [public benchmark mirror](https://polarishub.io/datasets/novartis/novartis-cyp3a4-v1).
[^15]: ChEMBL. “[Example CYP TDI assay record](https://www.ebi.ac.uk/chembl/explore/assay/CHEMBL2185189).” Accessed 2026.
[^16]: Zimmerlin A et al. “[Utility of in vitro TDI data for 400 reference drugs](https://pubmed.ncbi.nlm.nih.gov/21383203/).” 2011.
[^17]: Therapeutics Data Commons. “[ADME benchmark datasets](https://tdcommons.ai/single_pred_tasks/overview/).”
[^18]: OpenADMET. “[Octant CYP inhibition and reactivity release](https://huggingface.co/datasets/openadmet/Octant_CYP_inhibition_reactivity_blog_release).”
[^19]: SMARTCyp. “[Open-source CYP site-of-metabolism prediction](https://github.com/cdk/smartcyp).”
[^20]: Djoumbou-Feunang Y et al. “[BioTransformer 3.0](https://pmc.ncbi.nlm.nih.gov/articles/PMC9252798/).”
[^21]: Tian S et al. “[CyProduct: CYP-specific metabolic reaction prediction](https://pmc.ncbi.nlm.nih.gov/articles/PMC9032464/).”
