# 神原M論でのS/N

## ADU(Analog Digital Unit) の導出

| 文字             | 単位            | 意味                                           |
| ---------------- | --------------- | ---------------------------------------------- |
| $I_{s}$          | W / m$^2$ / str | 観測対象の放射輝度                             |
| $\eta$           | 無次元          | 量子効率（一つの光子を何個の電子に変換するか） |
| $\Omega$         | str             | 観測立体角                                     |
| $Gain$           | e$^-$ / DN      | システムゲイン                                 |
| $\tau _{\alpha}$ | 無次元          | 大気の透過率                                   |
| $\tau _t$        | 無次元          | 大気の透過率                                   |
| $\lambda$        | m               | 観測対象の波長                                 |
| $c$              | m / s           | 光速                                           |
|                  |                 |                                                |
| $ADU$            | DN / s          | 出力されるカウント値                           |
|                  |                 |                                                |

神原M論p19 式(1.1)より

$$
ADU
=
\frac
  {(I_s \cdot A \cdot \Omega) \cdot \tau_\alpha \cdot \tau_t \cdot \eta}
  { \left(h \cdot \cfrac{c}{\lambda} \right) \cdot Gain} \tag{1}
$$

ここで、$(I_s \cdot A \cdot \Omega) \cdot \tau_\alpha \cdot \tau_t$ は観測対象から検出器に到達するエネルギーを示し、 $(h \cdot \frac{c}{\lambda})$ はエネルギーを光子数に変換している。

## S/N の導出

| 文字       | 単位            | 意味                             |
| ---------- | --------------- | -------------------------------- |
| $t_{obs}$  | s               | 積分時間                         |
| $N_{read}$ | e$^-$ / pix     | 読み出しノイズ                   |
| $I_{dark}$ | e$^-$ / s / pix | 暗電流ノイズ                     |
| $n_{pix}$  | pix             | 観測対象が広がっているピクセル数 |
|            |                 |                                  |

神原M論p19 式(1.2)より

$$
\frac{S}{N}
=
\cfrac
  {ADU \cdot t_{obs}}
  {\sqrt{
    (ADU \cdot t_{obs})
    + \left( \cfrac{N_{read}}{Gain}\right)^2 n_{pix}
    + \left( \cfrac{I_{dark} \cdot t_{obs}}{Gain} \right) n_{pix}
    }
  } \tag{2}
$$

神原M論の原文では $ N_{read} $ と $I_{dark}$ に対して $Gain$ を処理していないが、引き継いだExcelファイル（観測可能性見積もり（神原さん作成）.xlsx）では処理しているので、執筆時の記入漏れだと考えられる

## 分からない点

- ADUの単位は [DN / s] なのか、[DN] なのか。
  - (1)式で使った文字の次元が正しければ $ADU$ [DN] になるはず
  - (2)式での分子分母を見比べると[DN / s]が正しいはず
- リードノイズの単位はExcelでの計算を見る限り [e$^-$ / pix] のはずで、撮像試験結果の解析の時に pix 当たりの値にするような処理をしていない。この実験値をそのまま式(2)に入れていいのか
- $n_{pix}$をかけるタイミングは $N_{read}$ や $I_{dark}$ を二乗した後で良いのか

## 解決した点

- なぜ分母に $ADU \cdot t_{obs}$ があるのか
  - 観測量 $S = ADU \cdot t_{obs}$ とポアソンノイズ $N$ には $N = \sqrt{S} = \sqrt{ADU \cdot t_{obs}}$ の関係
  - 互いに独立なノイズ $N_1$ と $N_2$ の和は $N_{all} = \sqrt{N_1^2 + N_2^2}$
  - よって観測量のポアソンノイズ $N$ とその他のノイズ$N_{other}$の和 $N_{all}$は以下式で表せる
  - $N_{all} = \sqrt{N^2 + N_{other}^2} = \sqrt{\sqrt{ADU \cdot t_{obs}}^2 + N_{other}^2} = \sqrt{ADU \cdot t_{obs} + N_{other}^2}$
- $I_{dark}$ の項が2乗されていないのはなぜか
  - 上記と同じ理由で、暗電流ノイズが $\sqrt{}$ の形だから。

# 野口M論でのS/N

## Signal導出

| 文字         | 単位      | 意味                                        |
| ------------ | --------- | ------------------------------------------- |
| $S_{photon}$ | e$^-$ / s | 観測対象の発光                              |
| $S_{GBT}$    | e$^-$ / s | 望遠鏡からの熱輻射（赤外発光）              |
| $S_{sky}$    | e$^-$ / s | 地球大気の赤外輻射                          |
| $I_{dark}$   | e$^-$ / s | 検出器のリーク電流（暗電流）                |
| $t_{obs}$    | s         | 露光時間                                    |
|              |           |                                             |
| $S_{all}$    | e$^-$     | ある露光時間 $t_{obs}$ での赤外輻射成分全体 |
|              |           |                                             |

検出器に入射する赤外輻射成分は野口M論p7 式(1.3)より

$$
S_{all}
=
S_{photon} \cdot t_{obs}
+ S_{sky} \cdot t_{obs}
+ S_{GBT} \cdot t_{obs} \tag{3}
$$

一方、野口M論p8では、
> なお、FPA自身のリーク電流は光信号と共に蓄積されるため、ある露光時間t(s)の全体のシグナル成分は、光学系で集光される成分にリーク電流成分を加えて以下のように表される。

とあり、野口M論p8 式(1.6)では

$$
S_{all}
=
S_{photon} \cdot t_{obs}
+ S_{sky} \cdot t_{obs}
+ S_{GBT} \cdot t_{obs}
+ I_{dark} \cdot t_{obs} \tag{4}
$$

と書かれている。また、
> Sphoton以外はノイズとなるが、ダーク画像、フラット画像、Sky画像を用いて取り除くことが可能である。

と書かれている。

## Noise導出

| 文字              | Signalとの関係                     | 単位        | 意味                        |
| ----------------- | ---------------------------------- | ----------- | --------------------------- |
| $N_{photon.shot}$ | $=\sqrt{S_{photon} \cdot t_{obs}}$ | e$^-_{rms}$ | $S_{photon}$ の統計的揺らぎ |
| $N_{sky.shot}$    | $=\sqrt{S_{sky} \cdot t_{obs}}$    | e$^-_{rms}$ | $S_{sky}$ の統計的揺らぎ    |
| $N_{GBT.shot}$    | $=\sqrt{S_{GBT} \cdot t_{obs}}$    | e$^-_{rms}$ | $S_{GBT}$ の統計的揺らぎ    |
| $N_{dark.shot}$   | $=\sqrt{I_{dark} \cdot t_{obs}}$   | e$^-_{rms}$ | $I_{dark}$ の統計的揺らぎ   |
| $N_{read}$        |                                    | e$^-_{rms}$ | 読み出しノイズ              |
|                   |                                    |             |
| $N_{FPA}$         |                                    | e$^-_{rms}$ | 検出器の全ノイズ成分        |
| $N_{all}$         |                                    | e$^-_{rms}$ | 全ノイズ成分                |
|                   |                                    |             |

それぞれの文字についての単位、意味はp7の文章を参考にまとめた。

Signalとの関係はp8の図1.7を参考にしたが、図中では
> Npshot = (Sphoton)1/2

と表記されているので、文章中の $S_{photon}$ を使って図の記載通りに解釈するなら $t_{obs}$ を掛け算する必要はない。しかし、この後p9で
> S\Nは、Sphoton・t(e)とNose all(erms)で表され、どちらも時間に対して単調増加する。（野口M論p9）

と言っているので、図1.7のSphotonとp7での $S_{photon}$ は実は単位が違うと考えるのが妥当。

> これらのノイズ成分は、平均値からのばらつき成分であり、お互いに創刊の無い複数のノイズが重畳する場合の全体のノイズは、エネルギー的に加算されるため、二乗加算平方根として求められる。（野口M論p7）

上記より、検出器のノイズはp7 式(1.4)より、
$$
N_{FPA} = \sqrt{(N_{dark.shot}) ^2 + (N_{read})^2} \tag{5}
$$

全ノイズ成分はp7 式(1.5)より、
$$
N_{all}
=
\sqrt{
  (N_{photon.shot})^2
  + (N_{sky.shot})^2
  + (N_{GBT.shot})^2
  + (N_{dark.shot})^2
  + (N_{read})^2
  } \tag{6}
$$

## S/Nの導出
> Sphoton以外はノイズとなるが、ダーク画像、フラット画像、Sky画像を用いて取り除くことが可能である。

ことから、p9 式(1.7)より、
$$
\frac{S}{N}
= \frac
  {S_{photon} \cdot t_{obs}}
  {N_{all}} \tag{7}
$$

と書かれている。（参考 p9）

よってこれを丁寧に書き下すなら

$$
\frac{S}{N}
= \frac
  {S_{photon} \cdot t_{obs}}
  {N_{all}}
= \cfrac
  {S_{photon} \cdot t_{obs}}
  {\sqrt{
    S_{photon} \cdot t_{obs}
    + S_{sky} \cdot t_{obs}
    + S_{GBT} \cdot t_{obs}
    + I_{dark} \cdot t_{obs}
    + (N_{read})^2
    }
  } \tag{8}
$$

となる。

## 分からない点

- 具体的に $N_{read}$ の中身はどうなっているのか
- ここでのNoiseはpixelの広がりをどう考慮しているのか
  - 単位を見る限り、各Noise成分はpixelの議論がすでに含まれていると考えるのが妥当
  - その場合、具体的に各Noise成分の中でpixelの議論をどう処理したのか

# 宇野M論でのS/N

## 観測対象の基本的な物理量

| 文字       | 単位                         | 意味                                        |
| ---------- | ---------------------------- | ------------------------------------------- |
| $F$        | W / m$^2$ / str / $\mu$m     | 観測対象の表面輝度                          |
| $\lambda$  | $\mu$m                       | 観測対象の波長                              |
| $h$        | m$^2$ kg / s                 | プランク定数                                |
| $\nu$      | Hz                           | photonの振動数                              |
| $c$        | m / s                        | 光速                                        |
|            |                              |                                             |
| $\epsilon$ | J                            | 波長 $\lambda$ のphoton 1個が持つエネルギー |
| $I$        | W / m$^2$ / str / $\mu$m / s | 観測対象が放出する1秒当たりのphotonの数     |
|            |                              |                                             |

文字の定義は宇野M論p74, p75を参考にまとめた。

波長 $\lambda$ のphoton 1個の持つエネルギ－ $\epsilon$ は、

$$
\epsilon = h \nu = hc /\lambda \tag{9}
$$

である。これより、ある波長 $\lambda$ で表面輝度 $F$ の対象が放出する1秒当たりのphotonの数は

$$
I = F / \epsilon \tag{10}
$$

と表せる。

## Signalの導出

| 文字             | 単位                         | 意味                                                         |
| ---------------- | ---------------------------- | ------------------------------------------------------------ |
| $A$              | m$^2$                        | 望遠鏡の有効口径                                             |
| $\Omega$         | str                          | 検出器の1pixelが見込む立体角                                 |
| $\Delta \lambda$ | $\mu$m                       | 検出器の1pixelが見込む波長幅                                 |
| $T$              | 無次元                       | 大気・望遠鏡光学系などの透過率の積算                         |
| $G$              | 無次元                       | エシェルグレーティング・クロスディスパーザーの回折効率の積算 |
| $\eta$           | e$^-$ / photon               | 検出器の量子効率                                             |
| $t$              | s                            | 観測での積分時間                                             |
| $I_{sky}$        | W / m$^2$ / str / $\mu$m / s | 地球大気や望遠鏡自身の輻射                                   |
|                  |                              |                                                              |
| $k$              | photon / s                   | 検出器の1pixelが1秒間に受け取るphotonの数                    |
| $S$              | e$^-$                        | 検出器1pixelで検出される光電子の数                           |
|                  |                              |                                                              |

文字の定義は宇野M論のp75を参考にまとめた。

式(10)の $I$ の観測において、検出器の1pixelが1秒間に受け取るphotonの数 $k$ は

$$
k = I A \Omega \Delta \lambda \tag{11}
$$

である。理想的な観測では大気や装置を通過し、検出器で電子として蓄積されるため、積分時間 $t$ の観測で蓄積される電子数 $S$ は

$$
S = k \cdot T G t \eta = I A \Omega \Delta \lambda T G t \eta \tag{12}
$$

である。実際の観測では、大気や望遠鏡自身の輻射である $I_{sky}$ が含まれるが、これは時間的に一定でありスカイ引きの処理で除去できるため、Signalとしては考慮に入れないことにする。

## S/Nの導出

| 文字       | 単位      | 意味                                   |
| ---------- | --------- | -------------------------------------- |
| $N_{read}$ | e$^-$     | 読み出しノイズ（積分時間に依存しない） |
| $I_{dark}$ | e$^-$ / s | 暗電流（積分時間に依存）               |
|            |           |                                        |
| $N$        | e$^-$     | ノイズの総和                           |
|            |           |                                        |

ノイズ成分の総和は
$$
N
= \sqrt{N_{read}^2 + I_{dark}t + I A \Omega \Delta \lambda T G t \eta + I_{sky} A \Omega \Delta \lambda T G t \eta}
= \sqrt{N_{read}^2 + \{I_{dark} + (I  + I_{sky}) A \Omega \Delta \lambda T G \eta\} t} \tag{13}
$$
となる。よって
$$
\frac{S}{N}
=
\cfrac
  {I A \Omega \Delta \lambda T G t \eta}
  {\sqrt{N_{read}^2 + \{I_{dark} + (I  + I_{sky}) A \Omega \Delta \lambda T G \eta\} t}} \tag{14}
$$
