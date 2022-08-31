# %%
import numpy as np

import object_oriented_estimation_myclass as ooem


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


if __name__ == "__main__":
    # Norwood et al. 2015 のfig2 に
    # WebPlotDigitizer https://automeris.io/WebPlotDigitizer/ でグラフから座標を抽出
    # その出力結果のcsvを読み出し
    # fig2は横軸：wavelength(um), 縦軸：Surface Brightness (Jy / sq.arcsec)

    rambda_um, Jy_per_sq_arcsec = read_csv(fpath="raw_data/Norwood_etal_2015_fig2.csv")

    # 波長の単位を m に変更
    rambda = rambda_um * 1e-6
