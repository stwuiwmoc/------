# %%
import importlib

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

    csv内のデータ構造は以下の通り、カラム名などは付いていないものを想定

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

    raw = np.loadtxt(fname=fpath, delimiter=",")
    x_data = raw.T[0]
    y_data = raw.T[1]

    return x_data, y_data


def convert_per_sq_arcsec_to_per_sr(hoge_per_sq_arcsec: np.ndarray) -> np.ndarray:
    """[何か / arcsec^2] を [何か / sr] に変換

    計算式の根拠は以下、但しこの関数では arcsec^2 と sr は分母側に来ることに注意

    oop観測見積もり.md
    └ 近赤外装置による撮像・分光 \n
        └ pixel数関連の導出 \n

    Parameters
    ----------
    hoge_per_sq_arcsec : np.ndarray
        [何か / arcsec^2] 1次元array

    Returns
    -------
    np.ndarray
        [何か / sr] 1次元array
    """

    hoge_per_sr = hoge_per_sq_arcsec / ((np.pi / 180) * (1 / 3600))**2

    return hoge_per_sr


if __name__ == "__main__":
    importlib.reload(ooem)

    # Norwood et al. 2015 のfig2 に
    # WebPlotDigitizer https://automeris.io/WebPlotDigitizer/ でグラフから座標を抽出
    # その出力結果のcsvを読み出し
    # fig2は横軸：wavelength(um), 縦軸：Surface Brightness (Jy / sq.arcsec)

    rambda_um, Jy_per_sq_arcsec = read_csv(fpath="raw_data/Norwood_etal_2015_fig2.csv")

    # 波長の単位を [um] -> [m] に変更
    rambda = rambda_um * 1e-6

    # 強度の単位を [Jy / arcsec^2] -> [Jy / sr] に変更
    Jy_per_sr = convert_per_sq_arcsec_to_per_sr(hoge_per_sq_arcsec=Jy_per_sq_arcsec)

    # 強度の単位を [Jy / sr] -> [W / m^2 / sr / m] に変更
    spectral_radiance = ooem.convert_Jy_per_sr_to_spectral_radiance(
        rambda_=rambda,
        Jy_per_sr_=Jy_per_sr)

    # csvに保存
    save_array = np.stack([rambda, spectral_radiance]).T
    save_fpath = mkfolder() + "rambda_vs_spectral_radiance.csv"
    np.savetxt(
        fname=save_fpath,
        X=save_array,
        delimiter=", ")

    # plot
    fig1 = plt.figure(figsize=(5, 10))
    gs1 = fig1.add_gridspec(2, 1)

    fig1.suptitle("Norwood et al. 2015, fig2 (2.5 - 5 um)")

    # pngの再現
    ax11 = fig1.add_subplot(gs1[0, 0])
    ax11.plot(rambda, Jy_per_sq_arcsec)
    ax11.set_yscale("log")
    ax11.grid()
    ax11.set_ylabel("[Jy / arcsec^2] ( = 10^-26 [W / m^2 / Hz / arcsec^2])")

    # 分光放射輝度に変換したもの
    ax12 = fig1.add_subplot(gs1[1, 0])
    ax12.plot(rambda, spectral_radiance)
    ax12.set_yscale("log")
    ax12.grid()
    ax12.set_xlabel("[m]")
    ax12.set_ylabel("[W / m^2 / sr / m]")

    fig1.tight_layout()
    fig1.savefig(mkfolder() + "fig1.png")
