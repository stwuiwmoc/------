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


if __name__ == "__main__":
    importlib.reload(ooem)
