# %%
import importlib
import math

import matplotlib.pyplot as plt
import numpy as np

import object_oriented_estimation_myclass as ooem


def mkfolder(suffix=""):
    import os
    """
    Parameters
    ----------
    suffix : str, optional
        The default is "".

    Returns
    -------
    str ( script name + suffix )
    """
    filename = os.path.basename(__file__)
    filename = filename.replace(".py", "") + suffix
    folder = "mkfolder/" + filename + "/"
    os.makedirs(folder, exist_ok=True)
    return folder


def read_csv(fpath: str) -> list[np.ndarray, np.ndarray]:
    """WebPlotDigitizer https://automeris.io/WebPlotDigitizer/ で
    グラフから抽出したデータのcsv出力を1次元arrayのlistとして読み出し

    csv内のデータ構造は
    5行目までは論文情報
    6行目はカラム名
    7行目以降はデータ
    6行目以降の構造は以下の通り

    x_column_name, y_column_name\n
    x1, y1\n
    x2, y2\n
    ...

    Parameters
    ----------
    fpath : str
        csvファイルのパス

    Returns
    -------
    list[np.ndarray, np.ndarray]
        [x座標の1次元array, y座標の1次元array]
    """

    raw = np.loadtxt(fname=fpath, delimiter=",", skiprows=6)
    x_data = raw.T[0]
    y_data = raw.T[1]

    return x_data, y_data


def valid_convolve(
        xx: np.ndarray,
        size: float) -> np.ndarray:
    """
    xxに対してsize個での移動平均を取る
    https://zenn.dev/bluepost/articles/1b7b580ab54e95 をコピペ

    Parameters
    ----------
    xx : np.ndarray
        移動平均をかける対象
    size : float
        移動平均の個数

    Returns
    -------
    np.ndarray
        移動平均後の結果
    """

    b = np.ones(size) / size
    xx_mean = np.convolve(xx, b, mode="same")

    n_conv = math.ceil(size / 2)

    # 補正部分
    xx_mean[0] *= size / n_conv
    for i in range(1, n_conv):
        xx_mean[i] *= size / (i + n_conv)
        xx_mean[-i] *= size / (i + n_conv - (size % 2))
    # size%2は奇数偶数での違いに対応するため

    return xx_mean


if __name__ == "__main__":
    importlib.reload(ooem)

    serial_name = "500Konly"
    Io_radius = 3643.2e3  # [m]

    # de Kleer et al. 2014 のfig1 に
    # WebPlotDigitizer https://automeris.io/WebPlotDigitizer/ でグラフから座標を抽出
    # その出力結果のcsvを読み出し
    # fig1は横軸：wavelength(um), 縦軸：Intensity (GW/um/sr)

    input_filepath = "raw_data/de_Kleer_etal_2014_fig2_" + serial_name + ".csv"
    rambda_um, intensity_GW_um = read_csv(fpath=input_filepath)

    # 波長を [um] から [m] に変換
    rambda = rambda_um * 1e-6

    # Fig1は黒体輻射の重ね合わせなので、本質的に滑らか
    # 抽出した座標はピクセルむらなどの理由で凸凹している
    # 移動平均をとって滑らかにする

    intensity_GW_um_smoothed = valid_convolve(
        xx=intensity_GW_um,
        size=10
    )

    # [GW/um/sr] -> [W/um/sr] の換算
    intensity_W_um = intensity_GW_um_smoothed * 1e9

    # [W/um/sr] -> [W/sr/m] の換算
    intensity_W_m = intensity_W_um * 1e6

    # plot
    fig1 = plt.figure(figsize=(10, 10))
    gs1 = fig1.add_gridspec(2, 1)

    ax11 = fig1.add_subplot(gs1[0, 0])
    ax11.plot(
        rambda,
        intensity_GW_um,
        label="raw")
    ax11.plot(
        rambda,
        intensity_GW_um_smoothed,
        label="smoothed")
    ax11.grid()
    ax11.legend()
    ax11.set_ylabel("Spectral Radiant Intensity [GW / um / sr]")

    ax12 = fig1.add_subplot(gs1[1, 0])

    fig1.tight_layout()
