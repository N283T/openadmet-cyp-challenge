# ヒトCYPの立体構造研究から見た基質認識・阻害・代謝選択性

## 結論：構造研究は何を説明でき、何をまだ説明できないか

ヒトの薬物代謝型CYPについて構造生物学が最もよく説明できるようになったのは、**「なぜある種類の分子がそのCYPに入りやすいのか」「どの官能基が認識のアンカーになるのか」「阻害剤がヘム鉄を直接配位しているか」「酵素がどの程度変形してリガンドを収容できるか」**である。一方、現在でも難しいのは、**「溶液・膜中で実際にどの結合姿勢が反応直前に選ばれるか」「複数の近接したC–H結合のうちどれが酸化されるか」「複数基質による協同性がどの占有状態から生じるか」「遺伝的変異が生体内クリアランスを何倍変えるか」**を、静的構造だけから予測することである。CYP1A2、2C9、2D6、3A4の比較は、この境界を非常にはっきり示す。citeturn5view0turn6view3turn23view0turn10view0

とくに重要なのは、**「結晶中でリガンドが見えた」ことを「その姿勢で代謝されることが証明された」と同一視しない**ことである。CYP2C9–warfarinの構造では、warfarinは予想外のポケットに観測され、著者自身が複数リガンド結合・アロステリーの可能性を提案した。一方、CYP2C9–flurbiprofenでは、Arg108によるカルボキシラート認識とヘムに対する配置が、NSAIDの選択性・位置選択性をかなり直接的に説明した。つまり同じ酵素の結晶構造でも、「結合を示す構造」と「触媒姿勢を強く支持する構造」では証拠の意味が違う。citeturn5view1turn6view3

CYP3A4はさらに極端である。2004年のprogesterone構造では基質がヘムから離れた周辺サイトに存在し、最初期の認識あるいはアロステリック部位が示唆された。ところが2006年には、ketoconazoleやerythromycinの結合によって活性部位体積が80%以上拡大し、二つのketoconazole分子や複数のerythromycin結合様式を許す構造が得られた。これは「CYP3A4の活性部位の大きさ」を一つの数値で表すこと自体が本質を外しており、**CYP3A4はポケットではなく“構造アンサンブル”として考えるべき**ことを示す代表例である。citeturn18view0turn10view0

逆に、阻害については構造だけでかなり強い結論に達する場合がある。CYP3A4–ritonavir [PDB 3NXU](https://www.rcsb.org/structure/3NXU) では、X線構造でritonavirのthiazole窒素がヘム鉄に直接配位しているだけでなく、同じ研究で結合、速度論、酸化還元電位が測定され、ritonavir結合による還元電位低下とP450 reductaseからの還元阻害が示された。これは「ヘム配位が見えた」だけでなく、**構造→電子状態→機能阻害**までつながった、非常に完成度の高い構造機能研究である。citeturn10view1

したがって、現状を大づかみに評価すると次のようになる。

| 問い | 構造研究の到達度 | 最大の注意点 |
|---|---|---|
| 基質が入るか／何を認識するか | **かなり説明できる** | 一つの構造を唯一の結合様式と思わない |
| 強い競合阻害・ヘム配位阻害 | **かなり説明できる場合がある** | 結合＝臨床DDIの大きさ、ではない |
| 位置選択性・代謝部位 | **一部は非常によく説明できる** | 近接距離だけでなく化学反応性が必要 |
| 立体選択性 | **系による** | 複数pose間の自由エネルギー差が小さい |
| CYP3A4の協同性 | **構造的基盤は見えるが完結していない** | どの多重占有状態が実際のkineticsを生むか未確定 |
| 変異による活性差 | **活性部位変異なら説明しやすい** | 遠隔変異、POR相互作用、発現・安定性は静的構造だけでは難しい |
| 生体膜中の基質進入・電子移動 | **まだ弱い** | 多くの高分解能構造が可溶化・短縮コンストラクト |
| 実際の生成物比率 | **構造だけでは不足** | 動力学、Compound I反応性、POR/b5、膜環境が必要 |

ここから先では、「どんな構造があるか」よりも、**各研究が何を問い、その問いに対してどの観測が得られ、どこまで因果的に言えるのか**を中心に整理する。

## 構造を読むための共通原理と証拠の強さ

CYPの薬物認識を考えるとき、少なくとも五つの段階を分けると文献の解釈がかなり明瞭になる。

**第一段階は単なる結合である。** 結晶電子密度に化合物が存在すれば、「この条件でこの場所を占有できる」ことは強く言える。しかし、それだけではそのposeが触媒反応中の主要poseとは限らない。CYP3A4のprogesterone [1W0F](https://www.rcsb.org/structure/1W0F) はまさにその例で、progesteroneは周辺サイトに観測され、著者らも初期認識あるいはアロステリックサイトという解釈を提示している。citeturn18view0

**第二段階はヘムに対する幾何学である。** 酸化される原子がヘム上方に適切に向いている、あるいは阻害剤のN原子などがFeに直接配位しているなら、機構的意味は強くなる。しかしP450の酸化剤は静止したFe(III)–hemeそのものではなく、触媒サイクル後半の高反応性鉄–オキソ種であるため、「最短Fe–C距離」だけで生成物選択性を決めることはできない。実際には基質poseの存在確率と、そのposeから各結合を切断する化学障壁の両方が効く。

**第三段階は分光学・電子状態との対応である。** CYP3A4–ritonavirではヘムへの直接配位という構造情報が、酸化還元電位やP450 reductaseからの電子移動への影響と対応したため、単なる“結晶中の接触”より強い阻害機構の証拠になっている。citeturn10view1

**第四段階は変異導入と生成物解析である。** ある残基を変えると結合親和性、生成物比、turnoverが予測通り変化するなら、「その残基が近い」という構造観察から「その残基が機能を担う」という因果関係に近づく。CYP2D6の最初の結晶構造 [2F9Q](https://www.rcsb.org/structure/2F9Q) はAsp301、Glu216、Phe120、Phe483などを従来の変異・モデリング結果と対応させ、後のリガンド複合体シリーズがその一部を直接可視化した。citeturn6view0turn23view0

**第五段階がnative-likeな反応系である。** ミクロソーム型CYPは本来ER膜に存在し、POR、場合によってはcytochrome b5と相互作用する。ところが代表的結晶構造の多くはE. coli発現の結晶化用コンストラクトである。例えばCYP1A2 2HI4はPDB登録配列495残基、CYP2D6の3TBG/4WNUは479残基、CYP3A4の代表構造は485残基で、結晶はnative ER膜・POR複合体そのものではない。CYP2C9の初期構造にはさらに複数の結晶化用置換が含まれる。したがって、**活性部位内部の原子配置について高精度でも、膜に対する向き、膜からの基質流入、PORとの会合を直接観測した構造ではない**という限界がある。citeturn5view0turn5view1turn6view3turn23view0turn18view0

代表的構造をこの観点で並べると次のようになる。

| CYP | 化合物・状態 | PDB | 分解能・コンストラクト上の注意 | 何がわかったか |
|---|---|---|---|---|
| CYP1A2 | α-naphthoflavone | [2HI4](https://www.rcsb.org/structure/2HI4) | 1.95 Å、PDB上mutationなし | 平面性リガンドに合う狭く閉じた空洞。citeturn5view0 |
| CYP2C9 | S-warfarin | [1OG5](https://www.rcsb.org/structure/1OG5) | 2.55 Å、PDB上7置換 | 予想外のポケット、多重結合・allosteryの可能性。citeturn5view1 |
| CYP2C9 | flurbiprofen | [1R9O](https://www.rcsb.org/structure/1R9O) | 2.0 Å、PDB上6置換 | Arg108–carboxylate認識と位置選択性。citeturn6view3 |
| CYP2D6 | apo | [2F9Q](https://www.rcsb.org/structure/2F9Q) | 3.0 Å、PDB上2置換 | Asp301/Glu216/Phe120などの認識モデルを構造化。citeturn6view0 |
| CYP2D6 | thioridazine ×1/2 | [3TBG](https://www.rcsb.org/structure/3TBG) ほか | 2.10 Å、PDB上mutationなし | open/closed、第二リガンド、access-channel antechamber。citeturn23view0 |
| CYP2D6 | quinidine | [4WNU](https://www.rcsb.org/structure/4WNU) | 2.26 Å、PDB上mutationなし | Asp/Gluだけでなく水媒介認識。citeturn23view1 |
| CYP3A4 | progesterone | [1W0F](https://www.rcsb.org/structure/1W0F) | 2.65 Å、PDB上mutationなし | 周辺認識/allosteric候補サイト。citeturn18view0 |
| CYP3A4 | metyrapone | [1W0G](https://www.rcsb.org/structure/1W0G) | 2.73 Å、PDB上mutationなし | 初期の阻害剤複合体、当時は比較的小さい活性部位像。citeturn18view1 |
| CYP3A4 | ketoconazole | [2V0M](https://www.rcsb.org/structure/2V0M) | 2.8 Å、PDB上mutationなし | 二分子結合、大規模な活性部位拡張。citeturn10view0 |
| CYP3A4 | erythromycin | [2J0D](https://www.rcsb.org/structure/2J0D) | 同研究 | 複数結合様式、柔軟な大型基質認識。citeturn10view0 |
| CYP3A4 | ritonavir | [3NXU](https://www.rcsb.org/structure/3NXU) | 2.0 Å、PDB上mutationなし | thiazole N–Fe配位とredox阻害を構造＋機能で確認。citeturn10view1 |

この表で重要なのは分解能の順位ではない。**触媒機構に対する「検証の深さ」が構造ごとに違う**ことである。

## 主要四酵素で構造と機能はどうつながったか

**CYP1A2：形状選択性はよく説明できるが、個々の代謝位置の説明はまだ間接的。**  
Sansenらの2007年構造 [2HI4](https://www.rcsb.org/structure/2HI4) は、ヒトCYP1ファミリーの薬物代謝型酵素を理解する上で大きな転機だった。α-naphthoflavoneがヘム遠位側に結合し、活性部位はCYP2・CYP3ファミリーと比べて**コンパクトで閉じ、平面性の高い芳香族化合物を収容する形**になっていた。これはCYP1A2がpolycyclic aromatic compoundのような比較的平板な分子を好むことを、単なる薬理学的pharmacophoreではなく三次元形状で説明した。PDB上では1.95 Åで、登録上アミノ酸置換なしである。citeturn5view0

ただし、ここには重要な非対称性がある。α-naphthoflavoneはこの構造では**阻害剤**であり、caffeineやphenacetinが実際に酸化される瞬間のposeが観測されたわけではない。そのため2HI4から強く言えるのは「CYP1A2の形状・芳香族面が平面性分子を選別する」というレベルであり、「caffeineのどのN-methyl基がなぜ主に酸化されるか」まで直接証明した構造ではない。CYP1A2は、**isoform selectivityの説明には構造が非常に有効だが、site-of-metabolismには結合アンサンブルと化学反応性を追加する必要がある**典型例といえる。citeturn5view0

**CYP2C9：warfarin構造は「結晶pose≠触媒pose」を教え、flurbiprofen構造はその逆を示した。**  
2003年のWilliamsらの [1OG5](https://www.rcsb.org/structure/1OG5) はS-warfarin複合体として歴史的に重要だが、単純な「warfarin 7-hydroxylationの写真」と読むべきではない。著者らが強調したのは、warfarinが予想外の結合ポケットを占めたことと、その配置が**multiple-ligand bindingやallostery**を示唆したことであった。PDBの結晶化コンストラクトには7か所の置換が登録されており、この点もnative enzymeの全状態をそのまま代表すると考えない方がよい。citeturn5view1

これに対し2004年のWesterらの flurbiprofen複合体 [1R9O](https://www.rcsb.org/structure/1R9O) は、より典型的な「構造から機能がわかる」例である。Arg108がflurbiprofenのカルボキシラートと直接相互作用し、Arg108自身はAsp293・Asn289とのネットワークに組み込まれている。その結果、基質がヘムに対して一定方向に配置され、CYP2C9が多くのアニオン性NSAIDを効率よく酸化する理由と位置選択性が説明された。またwarfarin構造との比較ではB–C領域やF–G領域の差が見られ、CYP2C9も固定ポケットではないことが明確になった。citeturn6view3

この二つを一緒に読むことが重要である。1OG5だけなら「結晶構造は触媒poseを外すことがある」、1R9Oだけなら「結晶構造から代謝選択性はかなり説明できる」と見える。両方が正しい。違いは、**電子密度でリガンドが見えるかではなく、そのposeが既知の生成物、ヘム近傍配置、残基相互作用、機能データと整合するか**である。citeturn5view1turn6view3

**CYP2D6：単一の「Asp301 pharmacophore」から、動的な酸性残基ネットワークへ。**  
最初のCYP2D6結晶構造 [2F9Q](https://www.rcsb.org/structure/2F9Q) はapo、3.0 Åで、薬物そのものは入っていなかった。それでもAsp301、Glu216、Phe120、Phe483などが既存のmutagenesisやhomology modelと対応し、特に多くのCYP2D6リガンドに共通するプロトン化塩基性窒素を酸性残基で認識するというモデルを構造的に裏付けた。Phe120は芳香環の向きを制御する候補として重要だった。citeturn6view0

その後のWangらの構造シリーズが、この単純モデルを大きく進めた。thioridazine複合体 [3TBG](https://www.rcsb.org/structure/3TBG) ではopen型とclosed型が得られ、active-site側のthioridazineのプロトン化piperidineがAsp301とcharge-assisted hydrogen bondを作る一方、open型では**第二のthioridazine分子が拡張したsubstrate-access channelのantechamberに存在し、Glu222と相互作用した**。これは「活性部位に一分子が入った完成状態」だけでなく、進入途中に相当し得る中間的占有位置まで構造として見えた非常に重要な例である。citeturn23view0

さらに同じ結晶系でthioridazineをprinomastat、quinidine、quinine、ajmalicineに置換すると、リガンドごとに酸性残基の使い方が変わった。quinineとajmalicineはGlu216と相互作用し、quinidine [4WNU](https://www.rcsb.org/structure/4WNU) のプロトン化NはAsp301とGlu216の間に置かれながら、**access channel内の水分子を介して水素結合**した。非イオン性のprinomastatは当然同じ塩橋を必要としない。さらにajmalicine結合ではGlu216の内向き移動を伴って結晶中でも酵素が閉じることが観測された。citeturn23view0turn23view1

この結果からCYP2D6について最も重要な修正は、**「塩基性NがAsp301に固定される」という一残基pharmacophoreでは不十分**ということである。実際にはAsp301、Glu216、Glu222、水、側鎖rotamer、F–G/B–C周辺の開閉が組み合わされ、同じ化学的特徴を持つ薬物でも異なる電荷ネットワークを利用できる。この構造シリーズは、rigid dockingで一つのAsp301–N距離を満たすことを高く評価するだけではCYP2D6を十分にモデル化できないことを、実験的に非常によく示している。citeturn23view0turn23view1

**CYP3A4：最大の発見は「大きなポケット」ではなく「大きく変わるポケット」である。**  
2004年のWilliamsらはapo、metyrapone、progesterone構造を報告した。[1W0F](https://www.rcsb.org/structure/1W0F) のprogesteroneは活性中心のヘム直上ではなくphenylalanine cluster上方の周辺部位に位置し、著者らはinitial substrate recognitionまたはallosteric effector siteの可能性を提示した。この時点の構造では、CYP3A4が非常に巨大な薬物を代謝する割には活性部位が意外に小さく、metyrapone/progesterone結合による構造変化も比較的小さかった。citeturn18view0turn18view1

ところがEkroosとSjögrenの2006年研究はその像を一変させた。[2V0M](https://www.rcsb.org/structure/2V0M) のketoconazole複合体と [2J0D](https://www.rcsb.org/structure/2J0D) のerythromycin複合体では、リガンドに応じて異なるopen conformationが形成され、活性部位体積が80%以上増大した。ketoconazoleは**二分子同時に収容**され、erythromycinには複数の結合様式が見いだされた。著者ら自身、この柔軟性は適切な実験構造がない状態での計算予測を難しくする、と指摘している。citeturn10view0

この2004年と2006年の“食い違い”は、むしろCYP3A4研究の最重要成果の一つである。1W0Fが間違いで2V0Mが正しいのではなく、**異なるリガンド・結晶条件が異なるエネルギー極小を選択した**と考えるべきである。つまり構造生物学が発見したのは「CYP3A4の正しい一個の構造」ではなく、「CYP3A4には複数の構造状態があり、それらの占有率をリガンドが変える」という性質である。citeturn18view0turn10view0

阻害ではritonavirが対照的に明快である。[3NXU](https://www.rcsb.org/structure/3NXU) は2.0 Åで、ritonavirのthiazole窒素がヘムFeに直接配位した。さらに同研究ではkinetics、equilibrium binding、redox measurementsが組み合わされ、結合によるP450の還元電位低下とPORからの還元抑制が確認された。したがって少なくともこの条件下の強い阻害成分については、**「分子が穴を塞ぐから阻害する」という曖昧な説明ではなく、Fe配位による電子状態変化まで構造的・物理化学的につながっている**。citeturn10view1

## 静的構造を越えて見えてきた柔軟性、水、進入経路、協同性、膜

**柔軟性は結晶学の“誤差”ではなくCYP機能そのものである。** CYP2D6では異なる結晶条件でthioridazine結合open/closed構造が得られ、さらに結晶内でリガンド交換を行っても別リガンドが活性部位へ入り、ajmalicineではタンパク質自体がよりclosedな構造へ移行した。これは少なくとも一部のCYP結晶では、結晶格子中でも基質経路と誘導適合に必要な運動性が残っていることを示す。もっとも、著者ら自身がlattice restraintの影響を議論しており、結晶中で観測される状態の存在比を溶液中の存在比とみなすことはできない。citeturn23view0

CYP3A4の結果はさらに大きい。metyrapone/progesteroneでは変化が小さかったのに、ketoconazole/erythromycinでは非常に大きく開く。したがって「apo構造を用いて10 ns程度MDをすれば自然にすべての薬物結合状態へ行く」と期待するのも危険である。実験構造間で既に大きな状態差があるからこそ、計算では一つの初期構造ではなく複数のreceptor states、十分な構造遷移sampling、あるいは実験拘束が必要になる。citeturn18view0turn10view0

**水分子は単なる結晶溶媒ではない。** CYP2D6–quinidineでは、塩基性NがAsp301に直接一対一で固定されるのではなく、水を介したネットワークとして観測された。この一例だけでも、ドッキング前処理でactive-site waterをすべて削除する操作が常に安全ではないことがわかる。一方、ritonavirのようなtype-II型の強いヘム配位では、リガンドのheteroatomがFeの遠位配位環境そのものを占有するため、水の有無は電子状態・阻害機構と直結する。citeturn23view1turn10view1

**基質進入経路も、単なる計算上の“トンネル”ではなく、少なくとも一部は実験構造として捕捉されている。** 3TBGの第二thioridazineは、expanded substrate-access channelのantechamberに位置した。この分子を必ず「基質進入中のスナップショット」と断定することはできないが、active siteの外側に化学的に意味のある酸性残基Glu222を伴う第二占有位置が存在し、そこから内部リガンドを含めて交換可能だった事実は、entry pathwayの構造モデルにかなり強い実験的制約を与える。citeturn23view0

**協同性については、構造的に“できそう”なことと、kineticsの原因が証明されたことを分ける必要がある。** CYP3A4ではperipheral progesterone siteと二分子ketoconazole結合が実際に観測されているので、「一つのCYP3A4分子に一分子しか結合できない」というモデルは明らかに不十分である。だが、例えば特定基質について観測されるsigmoidal kineticsが「1W0Fのprogesterone siteへのeffector結合」によって生じる、と一般化することまでは構造単独では証明されていない。異なる化合物では異なる多重占有状態・pose交換・構造状態が協同性を生む可能性がある。citeturn18view0turn10view0

**膜は、現在の高分解能構造を解釈するとき最大級の“欠けた環境変数”である。** ここで挙げた主要結晶構造はE. coli発現の結晶化用CYP catalytic domainを中心とするもので、native ER bilayer中でPORやcytochrome b5と複合体を形成した状態そのものではない。例えばCYP3A4 1W0F/1W0Gは485残基のコンストラクト、CYP2D6 3TBG/4WNUは479残基のコンストラクトとして登録されている。したがって結晶構造から見えるF–G領域の開口部が、ER膜中で同じ方向・同じ自由エネルギーで開くか、脂溶性基質が水相から来るのか膜相から来るのか、POR結合が活性部位アンサンブルをどの程度変えるかは別問題である。citeturn18view0turn18view1turn23view0turn23view1

ユーザーが挙げたcryo-EMについては、少なくとも今回、**CYP1A2/2C9/2D6/3A4の基質・阻害・生成物選択性を原子レベルで決定づけた代表研究として抽出できたものは依然としてX線結晶学が中心**だった。これはcryo-EMが原理的に役立たないという意味ではなく、むしろ今後はnative membrane、POR/CYP/b5の一過性複合体、複数conformationを扱える点で重要になる余地が大きい。現時点で「cryo-EM構造があるから生理的電子移動複合体が解決した」と扱うのは時期尚早であり、少なくとも主要薬物代謝CYPの成熟した構造機能エビデンスとは区別した方がよい。

## 計算手法はどこまで実験を説明し、どこで外れるか

CYPほど、**「ドッキングスコアがよい」ことと「正しい代謝機構を説明した」ことの距離が大きい標的は少ない**。理由は、結合親和性だけでなく、pose population、active-site flexibility、heme chemistry、基質固有の反応性、さらに水・膜・redox partnerが最終生成物に関与するからである。CYP2C9、2D6、3A4の構造シリーズはそのことを実験だけでかなり明確に示している。citeturn5view1turn6view3turn23view0turn10view0

**rigid dockingがうまくいきやすい条件**は、認識アンカーが強く、酵素変形が比較的小さく、一つのproductive orientationが他より明らかに有利な場合である。CYP2C9–flurbiprofenの「anion–Arg108」という相互作用や、CYP2D6の塩基性Nと酸性残基群という情報は、ドッキングに実験的constraintを与えるうえで非常に有用である。こうした条件では、構造は「どの向きを候補から落とすべきか」を強く絞れる。citeturn6view3turn23view0

ただしCYP2D6の後期構造が示したように、**正しいconstraintを一つに固定しすぎても外れる**。thioridazineではAsp301、quinine/ajmalicineではGlu216、access-channel側の第二thioridazineではGlu222、quinidineでは水媒介ネットワークが用いられる。したがって「basic N–Asp301距離が3 Å以内」という一個のpharmacophore constraintだけで候補を順位づけすれば、実験的に正しい結合様式の一部を誤って排除し得る。citeturn23view0turn23view1

**CYP3A4はrigid receptor dockingの最も厳しい反例である。** 1W0F型の比較的小さいポケットだけを受容体にすると、2V0Mで実際に観測された80%以上拡張した状態や二分子ketoconazole配置は原理的に表現しにくい。逆に、2V0Mの大きく開いた状態だけを使えば、小型リガンドに対するclosed stateを過小評価しうる。したがってCYP3A4ではensemble docking、induced-fit sampling、MDなどが「精度を少し上げるオプション」ではなく、**実験で確認された構造多型をモデルに入れるための要件**に近い。citeturn10view0turn18view0

**MDの最も説得力のある用途は、結晶構造の代用品ではなく、結晶構造間を結ぶ仮説を作ること**である。例えばCYP2D6のopen↔closed、Glu216 rotamer、access channel occupancy、CYP3A4のsmall-pocket↔expanded-pocketを初期・終状態あるいは実験拘束として使い、どの運動モードが遷移を支配するかを調べるのは理にかなう。一方、「一つのapo結晶から通常のMDを行ったら唯一の基質poseが得られた」という結果は、複数poseが実際に結晶学で観測されるCYP3A4やCYP2D6については特に慎重に解釈すべきである。citeturn23view0turn10view0

**QM/MMが必要になるのは、構造的には複数の代謝候補が残るところからである。** CYPの生成物選択性は、大まかには

\[
P(\text{生成物 }i)
\;\sim\;
\sum_j P(\text{pose}_j)\,
k_{ij}^{\text{chem}}
\]

と考えると理解しやすい。すなわち「その原子が反応中心に向いたposeがどのくらい存在するか」と「そのposeからその結合を酸化する化学障壁がどれくらい低いか」の積で決まる。一個のX線poseやdocking poseは主に前者への情報であり、後者には電子構造計算が必要になる。逆にQM/MMをどれほど精密にしても、出発poseが1OG5のような非productive/peripheral状態しか含んでいなければ、実際の生成物比を正しく説明する保証はない。CYP2C9–warfarinはこの点を考えるための非常に良いbenchmarkである。citeturn5view1

**構造予測モデルも同じ意味で“骨格仮説”と“反応機構”を区別すべきである。** 予測されたCYP foldは、未解決variantを既知WT構造へmapしたり、ポケット候補を比較したりするには有用である。しかし、ヘム鉄への配位、ordered water、複数リガンド占有、protonation/spin state、ligand-induced open state、native membrane orientationなどは、backbone predictionが高精度であることから自動的に正しくなるものではない。CYP3A4の実験構造群が示す通り、問題の中心はしばしば「foldを知ること」ではなく、「どの状態をどの条件で占有するか」にある。citeturn18view0turn10view0turn10view1

したがって計算結果の信頼度は、手法名よりも**何で検証されたか**で評価する方がよい。少なくとも、既知生成物の再現、site-directed mutantへの予測、spectral shift、binding affinity、同系統ligandの未使用データ、別の結晶構造へのprospective一致のどれかが欲しい。CYPでは「ドッキング図が化学的にもっともらしい」というだけでは、構造機能研究としては弱い。

## アイソフォーム差と変異：見た目の構造差で説明できる場合とできない場合

主要四酵素を並べると、アイソフォーム特異性にはかなり明瞭な構造的“個性”がある。

CYP1A2は**狭く、平面性化合物を好む閉じた空洞**、CYP2C9は**Arg108などを用いたアニオン性基質認識**、CYP2D6は**複数の酸性残基を使ったカチオン性窒素認識**、CYP3A4は**非常に可塑的で多重占有可能な空洞**として理解すると、多くの薬物選択性を定性的には説明できる。これは単なる配列identity比較より、実際の薬物化学との対応がよい。citeturn5view0turn6view3turn23view0turn10view0

しかし、このルールを「ポケット形状が似ていれば機能も似る」と逆向きに使うと危険である。CYP2C9のwarfarinとflurbiprofenだけでも、同じWT foldの中で異なるポケット、異なる開閉状態、異なる触媒的妥当性を持つposeが現れる。つまり**isoform間差より、同一isoform内のconformational ensemble差が大きい場合さえある**。citeturn5view1turn6view3

遺伝的variantではさらに注意が必要である。CYP2C9の代表的な*2のArg144→Cysのように活性部位中心から離れた表面側の変化と、*3のIle359→Leuのようにsubstrate-recognition領域へ近い変化では、同じ「代謝活性低下」でも構造的な意味が違う。前者のような場合は、単純に「ポケットが狭くなったから」と説明するより、タンパク質安定性・膜配置・PORとの相互作用・電子移動まで疑うべきである。後者では局所的なsubstrate orientation変化という説明がより自然になる。ただしWTの1R9Oなどへvariant位置をmapしてそう推測できることと、variantタンパク質自身の構造・turnover・POR依存性を直接測定したことは区別しなければならない。citeturn6view3

CYP2D6では同じ注意がさらに重要である。2F9Qや3TBGから活性部位variantを解釈することはできるが、臨床的なCYP2D6 phenotypeには、触媒ポケットとは別にタンパク質量、折り畳み、スプライシング、遺伝子コピー数なども関わる。したがって**「variantをAlphaFold/PDB上で可視化した」ことは「その患者の代謝表現型を構造的に説明した」ことではない**。構造が説明しやすいのは、精製variantタンパク質で結合・生成物比・電子移動の変化まで確認されたケースである。

この分野で特に有用な考え方は、variant効果を「活性部位までのÅ距離」だけで分類しないことである。P450では少なくとも、**substrate recognition、conformational dynamics、heme electronic state、POR docking、membrane orientation、protein stability**という異なる経路が同じ“活性低下”という実験結果に収束し得る。ritonavir研究が示すように、活性部位を物理的に塞がなくても電子移動能を変えるだけで触媒は止まり得るので、遠隔変異の解釈にも同じ発想が必要である。citeturn10view1

この意味で、「見た目の構造差が小さいのに活性差が大きい」という結果は構造生物学の失敗ではない。むしろ、**機能変数が静的座標以外にあることを知らせるシグナル**である。

## 優先して読むべき構造研究と、残っている重要な未解決問題

この分野を「構造一覧」ではなく「問いの進化」として読むなら、以下の順序が特に有益である。

**CYP1A2の形状選択性を理解する入口**はSansen et al., 2007, *JBC*, “Adaptations for oxidation of polycyclic aromatic hydrocarbons at the active site of human cytochrome P450 1A2”。[PDB 2HI4](https://www.rcsb.org/structure/2HI4)、[DOI](https://doi.org/10.1074/jbc.M611692200)。コンパクトな活性部位から「なぜCYP1A2らしい基質選択性が出るか」を考える原点としてよい。一方、阻害剤構造なので、これをそのままcaffeineなどのproductive poseに転用しないという訓練にもなる。citeturn5view0

**結晶poseと触媒poseの違いを学ぶ最重要論文**はWilliams et al., 2003, *Nature*, “Crystal Structure of Human Cytochrome P450 2C9 with Bound Warfarin”。[PDB 1OG5](https://www.rcsb.org/structure/1OG5)、[DOI](https://doi.org/10.1038/nature01862)。warfarinという臨床的に重要な基質の構造が得られたにもかかわらず、結果が単純な触媒配置にならず、多重結合・allosteryという新しい問いを生んだ点が重要である。citeturn5view1

それと対で読むべきなのがWester et al., 2004, *JBC*, CYP2C9–flurbiprofen。[PDB 1R9O](https://www.rcsb.org/structure/1R9O)、[DOI](https://doi.org/10.1074/jbc.M405427200)。Arg108を介したアニオン認識、基質配向、位置選択性という「構造が機能をきれいに説明した側」の代表である。citeturn6view3

**CYP2D6の基礎**にはRowland et al., 2006, *JBC*, [PDB 2F9Q](https://www.rcsb.org/structure/2F9Q)、[DOI](https://doi.org/10.1074/jbc.M511232200)。apo構造からAsp301、Glu216、Phe120などを従来のmutagenesis/pharmacophoreに結びつけた研究であり、「リガンド構造がなくても何がわかるか」を理解できる。citeturn6view0

その後に必ず読む価値があるのがWang et al., 2015, *JBC*, “Contributions of Ionic Interactions and Protein Dynamics to Cytochrome P450 2D6 Substrate and Inhibitor Binding”。[3TBG](https://www.rcsb.org/structure/3TBG)、[4WNU](https://www.rcsb.org/structure/4WNU) ほか、[DOI](https://doi.org/10.1074/jbc.M114.627661)。open/closed構造、第二リガンド、access antechamber、Asp301/Glu216/Glu222、水媒介相互作用、結晶内ligand exchangeが一つの研究にまとまっており、今回の問いに最も直接答える論文の一つである。citeturn23view0turn23view1

**CYP3A4は二本を対で読む必要がある。** まずWilliams et al., 2004, *Science*, “Crystal Structures of Human Cytochrome P450 3A4 Bound to Metyrapone and Progesterone”。[1W0F](https://www.rcsb.org/structure/1W0F)、[1W0G](https://www.rcsb.org/structure/1W0G)、[DOI](https://doi.org/10.1126/science.1099736)。progesteroneのperipheral siteと、当時の比較的小さい活性部位像が重要である。citeturn18view0turn18view1

続けてEkroos & Sjögren, 2006, *PNAS*, “Structural Basis for Ligand Promiscuity in CYP3A4”。[2V0M](https://www.rcsb.org/structure/2V0M)、[2J0D](https://www.rcsb.org/structure/2J0D)、[DOI](https://doi.org/10.1073/pnas.0603236103)。80%以上の活性部位拡張、二分子ketoconazole、erythromycinの複数poseを見ると、CYP3A4を一つのrigid pocketとして扱う考え方がなぜ破綻するかが非常によくわかる。citeturn10view0

**阻害機構を構造だけで終わらせなかった好例**としてSevrioukova & Poulos, 2010, *PNAS*, CYP3A4–ritonavir。[PDB 3NXU](https://www.rcsb.org/structure/3NXU)、[DOI](https://doi.org/10.1073/pnas.1010693107)。Fe–N配位、binding、kinetics、redoxを同じ論理で結んでいるため、「構造機能研究とは何をどこまで測れば強いのか」という観点で特に参考になる。citeturn10view1

残る未解決問題のうち、特に面白いものは以下だと思う。

**「productive poseは一つなのか、反応中に選ばれる一時的状態なのか」。** CYP2C9–warfarinとCYP3A4の複数poseが示すように、結晶で最もoccupancyの高い状態と最も速く反応する状態が同じとは限らない。反応速度は低存在率・高反応性poseに支配される可能性があり、これは静的構造だけでは解けない。citeturn5view1turn10view0

**「結合自由エネルギーと代謝選択性をどう切り分けるか」。** 強く結合するposeが必ず反応しやすいわけではない。今後重要なのは、構造ensembleの存在確率と、各状態からのCompound-I化学反応障壁を同時に扱い、実測product ratioへ接続することである。

**「CYP3A4協同性の構造単位は何か」。** 二分子ketoconazole、peripheral progesterone、複数erythromycin poseは多重占有の実在を示すが、異なる基質で得られる非Michaelis–Menten kineticsを一つの共通allosteric siteで説明できるとはまだ言えない。むしろリガンドごとに異なるoccupancy networkを持つ可能性が高い。citeturn18view0turn10view0

**「水はいつ入り、いつ抜けるのか」。** CYP2D6–quinidineの水媒介認識のように、水一分子がpharmacophoreを変える例がある。基質結合、Fe配位水の置換、O₂活性化、proton delivery、生成物放出までを同一系で時間分解して見ることができれば、現在の“結合構造”と“化学機構”の間の大きな空白が埋まる。citeturn23view1turn10view1

**「膜内での入口は結晶構造のchannelと同じなのか」。** CYP2D6のantechamberは有力な構造証拠だが、native ER bilayer中ではF–G領域やN末端anchorが膜と結合し、脂溶性基質の局所濃度も水溶液と違う。full-length CYPをnative-like lipid environmentに置き、基質の膜→酵素→活性部位経路を実験と計算の両方で追うことは、今後非常に価値が高い。citeturn23view0

**「POR/b5は電子を渡すだけなのか、それとも基質選択性を変えるallosteric partnerなのか」。** 多くの結晶構造にはredox partnerが存在しない一方、3NXU研究が示すようにCYPのredox energetics自体が阻害機構になり得る。substrate-dependentにPOR/b5との会合状態が変わり、それがconformational ensembleやcoupling efficiency、生成物比へ返ってくるかは、薬物代謝CYPを本当に生理的に理解するうえで大きな課題である。citeturn10view1

**「変異効果を静的構造以外の物理量に分解できるか」。** あるvariantの活性低下を、基質結合、productive-pose population、heme電子状態、POR electron transfer、膜配向、folding/stabilityに分けて測れれば、現在の「構造上この残基が近い／遠い」という説明から大きく進める。

総じて、この分野で最も確立した知見は「各CYPには特徴的なポケットがある」ということよりも、**そのポケットがリガンドごとに異なる形・水和状態・多重占有状態をとり、結合と触媒が必ずしも同じ構造状態で起こらない**ということである。CYP1A2では形状相補性、CYP2C9では電荷アンカーとproductive/nonproductive poseの区別、CYP2D6では酸性残基ネットワークとaccess-channel dynamics、CYP3A4では極端なplasticityとmultiple occupancyが、それぞれ構造研究によって非常によく確立された。citeturn5view0turn5view1turn6view3turn23view0turn10view0

一方で、**「どの原子が最終的に酸化されるか」を本当に予測するには、立体構造だけでなく、構造アンサンブル・水・膜・redox partner・反応性を同じ問題として扱う必要がある**。したがって今後もっとも価値が高いのは、単に新しいPDBを増やす研究より、同一CYP・同一基質について、複数構造状態、分光学、変異、turnover/product ratio、native-like membrane、必要に応じてMD/QM/MMを相互検証する研究である。CYP2D6のWangらの構造シリーズと、CYP3A4–ritonavirのSevrioukova/Poulos研究は、その方向性をすでにかなり明確に示している。citeturn23view0turn10view1