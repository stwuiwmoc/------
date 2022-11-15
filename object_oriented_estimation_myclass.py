import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import constants as phys_consts
from scipy import interpolate


def mkhelp(instance):
    import inspect
    attr_list = list(instance.__dict__.keys())
    for attr in attr_list:
        if attr.startswith("_"):
            continue
        print(attr)
    for method in inspect.getmembers(instance, inspect.ismethod):
        if method[0].startswith("_"):
            continue
        print(method[0] + "()")


def get_latest_commit_datetime() -> list[str]:
    """現在のツリーで直近のコミット日時を文字列として取得する

    Returns
    -------
    list[str]
        ["Latest commit datetime", コミットした日付, コミットした時刻]
    """

    import subprocess

    git_command = ["git", "log", "-1", "--format='%cI'"]

    proc = subprocess.run(
        git_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True)

    latest_commit_date_text = proc.stdout.decode("utf-8")

    # タイムゾーンを表す "+0900" を削除
    latest_commit_datetime_text_without_timezone = latest_commit_date_text[1:-8]

    # 日付のみを含む文字列
    latest_commit_date_text = latest_commit_datetime_text_without_timezone[:10]

    # 時刻のみを含む文字列
    latest_commit_time_text = latest_commit_datetime_text_without_timezone[11:]

    return "Latest commit datetime", latest_commit_date_text, latest_commit_time_text


def have_some_change_in_git_status() -> bool:
    """git的な意味でのファイルの変更の有無をboolianで取得する

    Returns
    -------
    bool
        ファイルの変更箇所ありならTrue, 無しならFalse
    """

    import subprocess

    git_command = ["git", "status", "--short"]

    proc = subprocess.run(
        git_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True)

    proc_text = proc.stdout.decode("utf-8")

    if len(proc_text) == 0:
        # git status --short で出力される文字列が無い
        # つまり「ファイルの変更箇所なし」を意味する
        have_some_change = False

    else:
        # git status --short で出力された文字列がある
        # つまり「ファイルに何かしらの変更箇所あり」を意味する
        have_some_change = True

    return have_some_change


def plot_parameter_table(
        fig: matplotlib.figure.Figure,
        position: matplotlib.gridspec.GridSpec,
        parameter_table: list,
        fontsize: int) -> matplotlib.axes._subplots.Axes:
    """パラメータ表示用のtableをax内に作成

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figureオブジェクト
    position: matplotlib.gridspec.GridSpec
        fig内での配置
    parameter_table : list
        パラメータリスト。縦横の要素数がそれぞれ一致している必要がある。（正方行列のように縦横で同じ要素数になる必要はない）
    fontsize : int
        テーブル内の文字サイズ

    Returns
    -------
    matplotlib.axes._subplots.AxesSubplot
        Axesオブジェクト
    """

    ax = fig.add_subplot(position)

    table = ax.table(
        cellText=parameter_table,
        loc="center")

    # fontsizeの調整
    table.auto_set_font_size(False)
    table.set_fontsize(fontsize)

    # 軸ラベルの消去
    ax.axis("off")

    # 縦幅をaxオブジェクトの縦幅に合わせて調整
    for pos, cell in table.get_celld().items():
        cell.set_height(1 / len(parameter_table))

    return ax


def calc_Plank_law_I_prime(
        rambda: np.ndarray,
        T: float) -> np.ndarray:
    """プランクの法則から波長と温度の関数として分光放射輝度を計算

    波長の1次元arrayに対応した、分光放射輝度の1次元arrayを返す

    oop観測見積もり.md
        └ 地球大気の発光 \n
            └ プランクの法則 \n

    Parameters
    ----------
    rambda : np.ndarray
        [m] 波長
    T : float
        [K] 黒体の温度

    Returns
    -------
    np.ndarray
        [W / m^2 / sr / m] 黒体放射による分光放射輝度
    """

    h = phys_consts.h
    c = phys_consts.c
    k_B = phys_consts.k

    I_prime = (2 * h * c**2 / rambda**5) * (1 / (np.exp(h * c / (rambda * k_B * T)) - 1))

    return I_prime


def convert_Jy_per_sr_to_spectral_radiance(
        rambda_: np.ndarray,
        Jy_per_sr_: np.ndarray) -> np.ndarray:
    """単位立体角あたりのジャンスキー [Jy / sr] を 分光放射輝度 [W / m^2 / sr / m] に変換

    oop観測見積もり.md
        └ 見積もりの概略 \n
            └ ジャンスキー（Jy）から分光放射輝度への換算 \n

    Parameters
    ----------
    rambda_ : np.ndarray
        [m] 波長の1次元array
    Jy_per_sr_ : np.ndarray
        [Jy / sr] 強度の1次元array

    Returns
    -------
    np.ndarray
        [W / m^2 / sr / m] 分光放射輝度の1次元array
    """

    c = phys_consts.c

    # 1 [Jy] = 10^-26 [W / m^2 / Hz]
    F_nu_per_sr = Jy_per_sr_ * 1e-26

    # 単位周波数あたり から 単位波長あたり に変換
    F_rambda_per_sr = (c / rambda_**2) * F_nu_per_sr

    return F_rambda_per_sr


class LightGenenrator:

    def __init__(
            self,
            rambda_lower_limit: float,
            rambda_upper_limit: float,
            rambda_division_width: float) -> None:
        """波長と分光放射輝度を保持するクラス

        oop観測見積もり.md
            └ 見積もりの概略 \n
                └ シミュレーション上での分光放射輝度の扱い \n

        Parameters
        ----------
        rambda_lower_limit : float
            [m] 計算する波長の最小値
        rambda_upper_limit : float
            [m] 計算する波長の最大値
        rambda_division_width : float
            [m] 計算する波長の分割幅
        """

        # 入力されたパラメータの代入
        self.__rambda_lower_limit = rambda_lower_limit
        self.__rambda_upper_limit = rambda_upper_limit
        self.__rambda_division_width = rambda_division_width

        # rambdaの1次元arrayを生成
        self.__rambda = np.arange(
            start=self.__rambda_lower_limit,
            stop=self.__rambda_upper_limit,
            step=self.__rambda_division_width)

        # 分光放射輝度の1次元arrayを生成
        self.__I_prime = np.zeros(len(self.__rambda))

    def h(self):
        mkhelp(self)

    def get_rambda_division_width(self) -> float:
        return self.__rambda_division_width

    def get_len(self) -> int:
        """rambdaの要素数（= I_prime の要素数）を返す

        Returns
        -------
        int
            rambdaの要素数
        """
        return len(self.__rambda)

    def get_rambda(self) -> np.ndarray:
        return self.__rambda

    def get_I_prime(self) -> np.ndarray:
        return self.__I_prime

    def add_I_prime_to(
            self,
            I_prime_xx: np.ndarray) -> None:
        """現在の分光放射輝度 I' に対して他の分光放射輝度を加算する

        Parameters
        ----------
        I_prime_xx : np.ndarray
            [W / m^2 / sr / m] 分光放射輝度の1次元のarray、
            要素数は self.get_len() と一致している必要がある。
        """
        self.__I_prime = self.__I_prime + I_prime_xx

    def multiply_I_prime_to(
            self,
            magnification: np.ndarray) -> None:
        """現在の分光放射輝度 I' に対して任意の倍率を乗算する

        Parameters
        ----------
        magnification : np.ndarray
            [無次元] 倍率（透過率など）の1次元のarray、
            要素数は self.get_len() と一致している必要がある。
        """
        self.__I_prime = self.__I_prime * magnification

    def get_I(self) -> float:
        """波長方向に self.__rambda_lower_limit から self.__upper_limit まで積分して放射輝度を計算

        Returns
        -------
        float
            [W / m^2 / sr] 放射輝度
        """
        # 各波長での微小面積要素の計算
        I_d_rambda = self.__I_prime * self.__rambda_division_width

        # 微小面積要素の足し合わせ
        I_ = np.sum(I_d_rambda)
        return I_

    def show_rambda_vs_I_prime_plot(
            self,
            fig: matplotlib.figure.Figure = None,
            position: matplotlib.gridspec.GridSpec = None) -> matplotlib.axes._subplots.Axes:
        """横軸波長、縦軸分光放射輝度のグラフを作成

        Parameters
        ----------
        fig : matplotlib.figure.Figure, optional
            Figureインスタンス, by default None
        position : matplotlib.gridspec.GridSpec, optional
            GridSpecインスタンス, by default None

        Returns
        -------
        matplotlib.axes._subplots.Axes
            axインスタンス
        """

        if fig is None:
            # figに明示的に代入しない場合は、単独のfigureとして作成
            fig = plt.figure()
        else:
            # figにfigureインスタンスが代入されている場合は、それをそのまま使う
            pass

        if position is None:
            # positionに明示的に代入しない場合は、
            # 単独のfigureの中に1つのaxだけがあるものとして作成
            ax = fig.add_subplot(111)
        else:
            # positionにgridspecが代入されている場合は、それをそのまま使う
            ax = fig.add_subplot(position)

        ax.plot(
            self.__rambda,
            self.__I_prime)

        ax.grid()
        ax.set_xlabel("rambda [m]")
        ax.set_ylabel("I' [W / m^2 / sr / m]")

        fig.tight_layout()

        return ax


class VirtualOutputFileGenerator:

    def __init__(self) -> None:
        """仮想的なFITSファイルとしてインスタンスを作成し、カウント値などを保持する

        oop観測見積もり.md
            └ 見積もりの概略 \n
                └ シミュレーション上でのFITSファイルの扱い \n

        """
        self.__S_all_pix = None
        self.__t_obs = None
        self.__n_bin_rambda = None
        self.__theta_pix = None
        self.__S_dark_pix = None
        self.__S_read_pix = None
        self.__R_electron_FW = None

    def h(self):
        mkhelp(self)

    def get_S_all_pix(self) -> int:
        return self.__S_all_pix

    def get_t_obs(self) -> float:
        return self.__t_obs

    def get_n_bin_rambda(self) -> float:
        return self.__n_bin_rambda

    def get_theta_pix(self) -> float:
        return self.__theta_pix

    def get_S_dark_pix(self) -> int:
        return self.__S_dark_pix

    def get_S_read_pix(self) -> int:
        return self.__S_read_pix

    def get_R_electron_FW(self) -> float:
        return self.__R_electron_FW

    def set_S_all_pix(self, S_all_pix: int) -> None:
        self.__S_all_pix = S_all_pix

    def set_t_obs(self, t_obs: float) -> None:
        self.__t_obs = t_obs

    def set_n_bin_rambda(self, n_bin_rambda: float) -> None:
        self.__n_bin_rambda = n_bin_rambda

    def set_theta_pix(self, theta_pix: float) -> None:
        self.__theta_pix = theta_pix

    def set_S_dark_pix(self, S_dark_pix: int) -> None:
        self.__S_dark_pix = S_dark_pix

    def set_S_read_pix(self, S_read_pix: int) -> None:
        self.__S_read_pix = S_read_pix

    def set_R_electron_FW(self, R_electron_FW: float) -> None:
        self.__R_electron_FW = R_electron_FW


class GenericEmissionFromCsv:

    def __init__(
            self,
            csv_fpath: str) -> None:
        """任意の発光についての "波長[m] vs 分光放射輝度[W / m^2 / sr / m]" のcsvを読み込んで
        線形補間処理して再現

        波長 vs 分光放射輝度 のcsvは自分で数字を用意するか、
        論文などのグラフ画像を用意して https://automeris.io/WebPlotDigitizer/ で
        データ抽出することを想定

        想定するcsvのデータ構造は以下（読み込みでは冒頭5行はスキップされる）

        Main Author, hogehoge\n
        Title, hogehoge\n
        Year, hogehoge\n
        URL, hogehoge\n
        Position, hogehoge\n
        波長0, 分光放射輝度0\n
        波長1, 分光放射輝度1\n
        波長2, 分光放射輝度2\n
        ...

        Parameters
        ----------
        csv_fpath : str
            csvのファイルパス
        """

        # 入力されたパラメータの代入
        self.__csv_fpath = csv_fpath

        # 分光放射輝度関数の作成
        self.__spectral_radiance_function = self.__make_spectral_radiance_function()

    def h(self):
        mkhelp(self)

    def get_csv_fpath(self) -> str:
        return self.__csv_fpath

    def get_spectral_radiance_function(self) -> interpolate.interpolate.interp1d:
        return self.__spectral_radiance_function

    def __make_spectral_radiance_function(self) -> interpolate.interpolate.interp1d:
        """波長 - 分光放射輝度 のcsvファイルを読み込んで分光放射輝度を線形補間し、
        波長を代入すると分光放射輝度を返す関数を作成する

        Returns
        -------
        interpolate.interpolate.interp1d
            関数（波長を入力するとその波長に対する分光放射輝度を出力する）
        """

        def get_rambda_and_spectral_radiance_from_csv(
                csv_filepath_: str) -> list[np.ndarray, np.ndarray]:
            """波長 - 分光放射輝度 のcsvファイルを読み込んで、
            波長と分光放射輝度それぞれ1次元のarrayを返す

            想定するcsvの構造は以下

            波長0, 分光放射輝度0\n
            波長1, 分光放射輝度1\n
            波長2, 分光放射輝度2\n
            ...

            ここで、波長の単位は [m]、分光放射輝度の単位は [W / m^2 / sr / m] を前提に計算する

            Parameters
            ----------
            csv_filepath_ : str
                csvのファイルパス

            Returns
            -------
            list[np.ndarray, np.ndarray]
                list[0] [m] 波長の1次元array,
                list[1] [W / m^2 / sr / m] 分光放射輝度の1次元array
            """

            raw = np.loadtxt(fname=csv_filepath_, delimiter=",")
            rambda_array_ = raw[:, 0]
            spectral_radiance_array_ = raw[:, 1]

            return rambda_array_, spectral_radiance_array_

        def calc_spectral_radiance_function(
                rambda_array_: np.ndarray,
                spectral_radiance_array_: np.ndarray) -> interpolate.interpolate.interp1d:
            """離散的な波長と分光放射輝度を線形補間し、
            波長を入れると分光放射輝度を返す関数を作る

            Parameters
            ----------
            rambda_array_ : np.ndarray
                [m] 波長の1次元array
            spectral_radiance_array_ : np.ndarray
                [W / m^2 / sr / m] 分光放射輝度の1次元array

            Returns
            -------
            interpolate.interpolate.interp1d
                関数（波長を入力するとその波長に対する分光放射輝度を出力する）
            """

            spectral_radiance_function_ = interpolate.interp1d(
                x=rambda_array_,
                y=spectral_radiance_array_,
                kind="linear")

            return spectral_radiance_function_

        # 波長 vs 分光放射輝度のcsvファイルを読み込む
        rambda_array, spectral_radiance_array = get_rambda_and_spectral_radiance_from_csv(
            csv_filepath_=self.__csv_fpath)

        # 関数にする
        spectral_radiance_function = calc_spectral_radiance_function(
            rambda_array_=rambda_array,
            spectral_radiance_array_=spectral_radiance_array)

        return spectral_radiance_function

    def add_spectral_radiance_to(
            self,
            light_instance: LightGenenrator) -> None:

        # 光の波長範囲・波長幅に対応した分光放射輝度の1次元arrayを計算する
        spectral_radiance = self.__spectral_radiance_function(light_instance.get_rambda())

        # 光に対して分光放射輝度を加える
        light_instance.add_I_prime_to(I_prime_xx=spectral_radiance)


class H3plusAuroralEmission:

    def __init__(
            self,
            rambda_obj: float,
            N_H3p: float,
            g_ns: int,
            J_prime: int,
            A_if: float,
            E_prime: float,
            T_hypo: float) -> None:

        """各輝線のパラメータから発光輝線強度 I_obj を導出

        oop観測見積もり.md
            └ 観測対象の発光 \n
                ├ H3+輝線の放射輝度 \n
                └ シミュレーション上での実装 \n

        Parameters
        ----------
        rambda_obj : float
            [m] 輝線の中心波長
        N_H3p : float
            H3+のカラム密度
        g_ns : int
            [無次元] nuclear spin weight, 2 or 4
        J_prime : int
            [無次元] 回転準位
        A_if : float
            [/s] アインシュタインのA係数
        E_prime : float
            [/cm] energy of upper statement
        T_hypo : float
            [K] 想定するH3+温度
        """

        # 入力されたパラメータの代入
        self.__rambda_obj = rambda_obj
        self.__N_H3p = N_H3p
        self.__g_ns = g_ns
        self.__J_prime = J_prime
        self.__A_if = A_if
        self.__E_prime = E_prime
        self.__T_hypo = T_hypo

        # その他のパラメータの計算
        self.__omega_if = 1 / self.__rambda_obj * 1e-2  # 波数 [/cm]
        self.__Q_T = self.__calc_Q_T()

        # 発光輝線強度 I_obj の計算
        self.__I_obj = self.__calc_I_obj()

    def h(self):
        mkhelp(self)

    def __calc_Q_T(self) -> float:
        """Partition function Q(T) の導出

        Returns
        -------
        float
            partition function Q(T)
        """
        T = self.__T_hypo

        A_0 = - 1.11391
        A_1 = + 0.0581076
        A_2 = + 0.000302967
        A_3 = - 2.83724e-7
        A_4 = + 2.31119e-10
        A_5 = - 7.15895e-14
        A_6 = + 1.00150e-17

        Q_T = A_0 * T**0 + A_1 * T**1 + A_2 * T**2 + A_3 * T**3 + A_4 * T**4 + A_5 * T**5 + A_6 * T**6
        return Q_T

    def __calc_I_obj(self) -> float:
        """H3+輝線の発光強度を計算

        oop観測見積もり.md
            └ 観測対象の発光 \n
                └ H3+輝線の放射輝度 \n

        Returns
        -------
        float
            [W / m^2 / sr] 輝線の発光強度
        """

        N_H3p = self.__N_H3p
        g_ns = self.__g_ns
        J_prime = self.__J_prime
        h = phys_consts.h
        c = phys_consts.c
        omega_if = self.__omega_if
        A_if = self.__A_if
        E_prime = self.__E_prime
        k_b = phys_consts.k
        T_hypo = self.__T_hypo
        Q_T = self.__Q_T

        I_obj = N_H3p * g_ns * (2 * J_prime + 1) * h * c * (omega_if * 1e2) * A_if \
            * np.exp(- h * c * (E_prime * 1e2) / (k_b * T_hypo)) \
            / (4 * np.pi * Q_T)

        return I_obj

    def get_rambda_obj(self) -> float:
        return self.__rambda_obj

    def get_N_H3p(self) -> float:
        return self.__N_H3p

    def get_g_ns(self) -> int:
        return self.__g_ns

    def get_J_prime(self) -> int:
        return self.__J_prime

    def get_A_if(self) -> float:
        return self.__A_if

    def get_E_prime(self) -> float:
        return self.__E_prime

    def get_T_hypo(self) -> float:
        return self.__T_hypo

    def add_auroral_emission_to(
            self,
            light_instance: LightGenenrator) -> None:

        def find_rambda_obj_index(
                rambda_: np.ndarray,
                rambda_obj_: float) -> int:
            """rambda_objに最も近いrambdaのindexを探す

            Parameters
            ----------
            rambda_ : np.ndarray
                LightGenerator.get_rambda()
            rambda_obj_ : float
                self.rambda_obj

            Returns
            -------
            int
                rambda_objに最も近いrambdaのindex
            """

            # rambda_obj と最も近いrambdaを探したい
            # -> 差分の絶対値が最も小さいrambdaのインデックスを求めればよい
            rambda_diff_abs = np.abs(rambda_ - rambda_obj_)
            rambda_obj_index = np.argmin(rambda_diff_abs)
            return rambda_obj_index

        rambda_obj_index = find_rambda_obj_index(
            rambda_=light_instance.get_rambda(),
            rambda_obj_=self.__rambda_obj)

        # 実際には輝線幅≃0を想定した見積もりだが、実装上は分光放射輝度に直す必要がある
        # 輝線の中心波長での分光放射輝度を計算
        I_prime_obj_in_center_wavelength: float = self.__I_obj / light_instance.get_rambda_division_width()

        # 輝線の中心波長で上で計算した分光放射輝度、それ以外では値0の1次元arrayを作成
        I_prime_obj = np.zeros(light_instance.get_len())
        I_prime_obj[rambda_obj_index] = I_prime_obj_in_center_wavelength

        # LightGeneratorのインスタンスに、I_prime_objを加える
        light_instance.add_I_prime_to(I_prime_xx=I_prime_obj)


class EarthAtmosphere:

    def __init__(
            self,
            T_ATM: float,
            ATRAN_result_filepath: str) -> None:
        """地球大気の透過率と熱輻射の計算

        oop観測見積もり.md
            └ 地球大気の発光 \n
                ├ 地球大気のパラメータ \n
                ├ ATRANによる大気透過率の計算 \n
                └ 地球大気の熱輻射による分光放射輝度 \n

        Parameters
        ----------
        T_ATM : float
            [K] 大気の温度
        ATRAN_result_filepath : str
            ATRAN計算結果をコピペした.txtファイルパス（ATRANの出力をCtrl+Sでtxtとして保存したもの）
        """

        # 入力されたパラメータの代入
        self.__T_ATM = T_ATM
        self.__ATRAN_result_filepath = ATRAN_result_filepath

        # 透過率関数の作成
        self.__tau_ATM_function = self.__make_tau_ATM_function()

    def h(self):
        mkhelp(self)

    def get_T_ATM(self) -> float:
        return self.__T_ATM

    def get_ATRAN_result_filepath(self) -> str:
        return self.__ATRAN_result_filepath

    def get_tau_ATM_function(self) -> interpolate.interpolate.interp1d:
        return self.__tau_ATM_function

    def __make_tau_ATM_function(self) -> interpolate.interpolate.interp1d:
        """ATRAN出力結果のtxtファイルを読み込み、
        波長を代入すると大気透過率を返す関数を作成する

        oop観測見積もり.md
            └ 地球大気の発光 \n
                ├ 地球大気のパラメータ \n
                └ ATRANによる大気透過率の計算 \n

        Returns
        -------
        interpolate.interpolate.interp1d
            関数（波長を入力するとその波長に対する大気透過率を出力する）
        """

        def get_ATRAN_rambda_and_ATRAN_tau_ATM(
                ATRAN_result_filepath_: str) -> list[np.ndarray, np.ndarray]:
            """ATRANの計算結果のファイルを読み込んで波長と透過率それぞれ1次元のarrayを返す

            Parameters
            ----------
            ATRAN_result_filepath_ : str
                ATRANの計算結果のファイルパス（.txt）

            Returns
            -------
            list[np.ndarray, np.ndarray]
                list[0] [m] 波長の1次元array,
                list[1] [無次元] 透過率の1次元array
            """

            # 想定する.txtの内部構造は空白区切りで
            # idx0 波長um0 透過率0
            # idx1 波長um1 透過率1
            # idx2 波長um2 透過率2
            # ...
            # という形（ATRANの出力をCtrl+Sでtxtとして保存したもの）

            raw = np.loadtxt(fname=ATRAN_result_filepath_)

            ATRAN_rambda_array_um = raw[:, 1]  # ATRANで出力される波長はum単位
            ATRAN_tau_ATM_array_ = raw[:, 2]

            ATRAN_rambda_array_ = ATRAN_rambda_array_um * 1e-6  # 波長の単位をmに直す
            return ATRAN_rambda_array_, ATRAN_tau_ATM_array_

        def calc_tau_ATM_function(
                ATRAN_rambda_array_: np.ndarray,
                ATRAN_tau_ATM_array_: np.ndarray) -> interpolate.interpolate.interp1d:
            """ATRANの離散的な出力（波長と透過率）を線形補間し、波長を入れると透過率を返す関数を作る

            Parameters
            ----------
            ATRAN_rambda_array_ : np.ndarray
                [m] ATRANの出力した計算結果のうち、波長をum -> mに単位換算した1次元array
            ATRAN_tau_ATM_array_ : np.ndarray
                [透過率] ATRANの出力した計算結果のうち、透過率の1次元array

            Returns
            -------
            interpolate.interpolate.interp1d
                関数（波長を入力するとその波長に対する大気透過率を出力する）
            """

            tau_ATM_function_ = interpolate.interp1d(
                x=ATRAN_rambda_array_,
                y=ATRAN_tau_ATM_array_,
                kind="linear")

            return tau_ATM_function_

        # 透過率計算結果のファイルを読み込む
        ATRAN_rambda_array, ATRAN_tau_ATM_array = get_ATRAN_rambda_and_ATRAN_tau_ATM(
            ATRAN_result_filepath_=self.__ATRAN_result_filepath)

        # 関数にする
        tau_ATM_function = calc_tau_ATM_function(
            ATRAN_rambda_array_=ATRAN_rambda_array,
            ATRAN_tau_ATM_array_=ATRAN_tau_ATM_array)

        return tau_ATM_function

    def pass_through(
            self,
            light_instance: LightGenenrator) -> None:
        """光に大気透過率をかけ、大気の熱輻射による分光放射輝度を加える

        Parameters
        ----------
        light_instance : LightGenenrator
            自作インスタンス
        """

        def calc_I_prime_ATM(
                rambda_: np.ndarray,
                tau_ATM_: np.ndarray) -> np.ndarray:
            """望遠鏡の熱輻射による分光放射輝度を計算

            oop観測見積もり.md
                └ 地球大気の発光 \n
                    └ 地球大気の熱輻射による分光放射輝度 \n

            Parameters
            ----------
            rambda_ : np.ndarray
                LightGenerator.get_rambda()
            tau_ATM_ : np.ndarray
                [無次元] 大気透過率の1次元array

            Returns
            -------
            np.ndarray
                [W / m^2 / sr / m] 分光放射輝度
            """
            T_ATM = self.__T_ATM

            I_prime_ATM_ = (1 - tau_ATM_) * calc_Plank_law_I_prime(rambda=rambda_, T=T_ATM)

            return I_prime_ATM_

        # 光の波長範囲・波長幅に対応した大気透過率の1次元arrayを計算する
        tau_ATM = self.__tau_ATM_function(light_instance.get_rambda())

        # 光に対して大気の透過率をかける
        light_instance.multiply_I_prime_to(magnification=tau_ATM)

        # 光に対して大気の熱輻射による分光放射輝度を加える
        I_prime_ATM = calc_I_prime_ATM(
            rambda_=light_instance.get_rambda(),
            tau_ATM_=tau_ATM)

        light_instance.add_I_prime_to(I_prime_xx=I_prime_ATM)


class GroundBasedTelescope:

    def __init__(
            self,
            D_GBT: float,
            FNO_GBT: float,
            T_GBT: float,
            tau_GBT: float) -> None:

        """望遠鏡のパラメータを保持

        oop観測見積もり.md
            └ 望遠鏡の発光 \n
                ├ 望遠鏡光学系のパラメータ \n
                └ 望遠鏡の熱輻射による分光放射輝度 \n

        Parameters
        ----------
        D_GBT : float
            [m] 望遠鏡主鏡の口径
        FNO_GBT : float
            [無次元] 望遠鏡光学系の合成F値
        T_GBT : float
            [K] 望遠鏡光学系の温度
        tau_GBT : float
            [無次元] 望遠鏡光学系の透過率
        """

        # 入力されたパラメータの代入
        self.__D_GBT = D_GBT
        self.__FNO_GBT = FNO_GBT
        self.__T_GBT = T_GBT
        self.__tau_GBT = tau_GBT

        self.__f_GBT = self.__D_GBT * self.__FNO_GBT
        self.__A_GBT = np.pi * (self.__D_GBT / 2) ** 2

    def h(self) -> None:
        mkhelp(self)

    def get_D_GBT(self) -> float:
        return self.__D_GBT

    def get_FNO_GBT(self) -> float:
        return self.__FNO_GBT

    def get_T_GBT(self) -> float:
        return self.__T_GBT

    def get_tau_GBT(self) -> float:
        return self.__tau_GBT

    def get_f_GBT(self) -> float:
        return self.__f_GBT

    def get_A_GBT(self) -> float:
        return self.__A_GBT

    def pass_through(
            self,
            light_instance: LightGenenrator) -> None:
        """光に望遠鏡透過率をかけ、望遠鏡の熱輻射による分光放射輝度を加える

        Parameters
        ----------
        light_instance : LightGenenrator
            自作インスタンス
        """

        def calc_I_prime_GBT(
                rambda_: np.ndarray) -> np.ndarray:
            """望遠鏡の熱輻射による分光放射輝度を計算

            oop観測見積もり.md
                └ 望遠鏡の発光 \n
                    └ 望遠鏡の熱輻射による分光放射輝度 \n

            Parameters
            ----------
            rambda_ : np.ndarray
                LightGenerator.get_rambda()

            Returns
            -------
            np.ndarray
                [W / m^2 / sr / m] 分光放射輝度
            """

            T_GBT = self.__T_GBT
            tau_GBT = self.__tau_GBT

            I_prime_GBT_ = (1 - tau_GBT) * calc_Plank_law_I_prime(rambda=rambda_, T=T_GBT)

            return I_prime_GBT_

        # 光に対して望遠鏡の透過率をかける
        light_instance.multiply_I_prime_to(magnification=self.__tau_GBT)

        # 光に対して望遠鏡の熱輻射による分光放射輝度を加える
        I_prime_GBT = calc_I_prime_GBT(rambda_=light_instance.get_rambda())
        light_instance.add_I_prime_to(I_prime_xx=I_prime_GBT)


class ImagingInstrument:

    def __init__(
            self,
            is_TOPICS: bool,
            rambda_BPF_center: float,
            FWHM_BPF: float,
            tau_BPF_center: float,
            tau_i_ND: float,
            G_Amp: float,
            I_dark: float,
            N_e_read: float) -> None:
        """撮像装置のパラメータを保持

        oop観測見積もり.md
            └ 近赤外装置による撮像・分光 \n
                ├ バンドパスフィルタ透過率の導出 \n
                ├ NDフィルタ透過率の導出 \n
                ├ 装置透過率の導出（TOPICS） \n
                ├ 装置透過率の導出（ESPRIT 分光モード） \n
                ├ pixel数関連の導出 \n
                ├ システムゲインの導出 \n
                ├ 検出器に到達した段階での分光放射輝度 \n
                ├ Signalへの換算 \n
                └ フルウェルによるカウント値の制限 \n

        Parameters
        ----------
        is_TOPICS : bool
            TOPICSの場合はTrue, ESPRITの場合はFalse
        rambda_BPF_center : float
            [m] バンドパスフィルタの中心波長
        FWHM_BPF : float
            [m] バンドパスフィルタの半値全幅
        tau_BPF_center : float
            [無次元] バンドパスフィルタの中心透過率
        tau_i_ND : float
            [無次元] NDフィルタの透過率
        G_Amp : float
            [無次元] プリアンプ基板の倍率
        I_dark : float
            [e- / s / pix] 検出器暗電流
        N_e_read : float
            [e-rms / pix] 駆動回路の読み出しノイズ
        """

        # --- 入力パラメータ・固定パラメータの代入 ---
        # TPOICSかESPRITの選択パラメータ
        self.__is_TOPICS = is_TOPICS

        # バンドパスフィルタ透過率に関するパラメータ
        self.__rambda_BPF_center = rambda_BPF_center
        self.__FWHM_BPF = FWHM_BPF
        self.__tau_BPF_center = tau_BPF_center

        # NDフィルタのパラメータ
        self.__tau_i_ND = tau_i_ND

        # ピクセル数関連のパラメータ
        self.__n_bin_rambda = 1  # [pix]

        # システムゲインに関連するパラメータ
        self.__G_Amp = G_Amp
        self.__G_sys = self.__calc_G_sys()
        self.__FW = 152000  # [e- / pix]
        self.__S_FW_pix = self.__FW / self.__G_sys

        # Signalへの換算で必要なパラメータ
        self.__I_dark = I_dark
        self.__N_e_read = N_e_read

        # --- 望遠鏡への設置状況 ---
        # インスタンス生成時点では None としておく
        # インスタンス生成後にメソッドを使って書き換える
        self.__GBT_instance = None

    def h(self):
        mkhelp(self)

    def __calc_G_sys(self) -> float:
        """システムゲインの導出

        oop観測見積もり.md
            └ 近赤外装置による撮像・分光 \n
                └ システムゲインの導出 \n

        Returns
        -------
        float
            [e- / DN] システムゲイン
        """

        C_PD = 7.20e-14  # [F]
        ADU_ADC = 10 / 2**16  # [V / DN]
        G_SF = 0.699  # [無次元]
        G_Amp = self.__G_Amp
        e = phys_consts.e

        G_sys = (C_PD / (e * G_SF)) * (ADU_ADC / G_Amp)
        return G_sys

    def get_is_TOPICS(self) -> bool:
        return self.__is_TOPICS

    def get_rambda_BPF_center(self) -> float:
        return self.__rambda_BPF_center

    def get_FWHM_BPF(self) -> float:
        return self.__FWHM_BPF

    def get_tau_BPF_center(self) -> float:
        return self.__tau_BPF_center

    def get_tau_i_ND(self) -> float:
        return self.__tau_i_ND

    def get_G_Amp(self) -> float:
        return self.__G_Amp

    def get_G_sys(self) -> float:
        return self.__G_sys

    def get_FW(self) -> float:
        return self.__FW

    def get_S_FW_pix(self) -> float:
        return self.__S_FW_pix

    def get_I_dark(self) -> float:
        return self.__I_dark

    def get_N_e_read(self) -> float:
        return self.__N_e_read

    def get_n_bin_rambda(self) -> float:
        return self.__n_bin_rambda

    def set_ImagingInstrument_to(
            self,
            GBT_instance: GroundBasedTelescope) -> None:
        """望遠鏡に撮像装置を設置し、GBTのインスタンスを内部で呼び出せるようにする

        Parameters
        ----------
        GBT_instance : GroundBasedTelescope
            GroundBasedTelescopeクラスのインスタンス
        """

        def calc_theta_pix(f_GBT: float) -> float:
            """1pixelが見込む角度（プレートスケール）を計算

            oop観測見積もり.md
                └ 近赤外装置による撮像・分光 \n
                    └ pixel数関連の導出 \n

            Parameters
            ----------
            f_GBT : float
                GroundBasedTelescope.get_f_GBT()

            Returns
            -------
            float
                [arcsec / pix] プレートスケール
            """

            # 装置内部光学系の倍率
            if self.__is_TOPICS:
                # TOPICSの場合
                m_i_all = 2  # [無次元]
            else:
                # ESPRITの場合
                m_i_all = 1  # [無次元]

            s_pix = 30e-6  # [m / pix]

            theta_pix = np.arctan(s_pix / (m_i_all * f_GBT)) * (180 / np.pi) * 3600
            return theta_pix

        self.__GBT_instance = GBT_instance

        # pixel数関連の導出
        self.__theta_pix = calc_theta_pix(
            f_GBT=self.__GBT_instance.get_f_GBT())
        self.__Omega_pix = ((self.__theta_pix / 3600) * (np.pi / 180)) ** 2

    def is_set_to_GBT_instance(self) -> bool:
        if self.__GBT_instance is None:
            return False
        else:
            return True

    def get_theta_pix(self):
        if self.is_set_to_GBT_instance() is False:
            print("Method 'set_ImagingInstrument_to()' is not used yet.")
            return None

        else:
            return self.__theta_pix

    def get_Omega_pix(self):
        if self.is_set_to_GBT_instance() is False:
            print("Method 'set_ImagingInstrument_to()' is not used yet.")
            return None

        else:
            return self.__Omega_pix

    def shoot_light_and_save_to_fits(
            self,
            light_instance: LightGenenrator,
            virtual_output_file_instance: VirtualOutputFileGenerator,
            t_obs: float) -> None:
        """光を受け取って、カウント値を計算し、仮想的なfitsファイルに値を保存

        oop観測見積もり.md
            └ 近赤外装置による撮像・分光 \n
                ├ 検出器に到達した段階での分光放射輝度 \n
                └ Signalへの換算 \n

        Parameters
        ----------
        light_instance : LightGenenrator
            自作インスタンス
        virtual_output_file_instance : VirtualOutputFileGenerator
            自作インスタンス
        t_obs : float
            [s] 積分時間、floatでも1次元のarrayでもよい
        """

        def calc_gaussian(
                y_max: float,
                x_array: np.ndarray,
                FWHM_: float,
                x_center: float) -> np.ndarray:
            """横軸x、縦軸y としてガウシアンを計算

            oop観測見積もり.md
                └ 近赤外装置による撮像・分光 \n
                    └ バンドパスフィルタ透過率の導出 \n

            Parameters
            ----------
            y_max : float
                yの最大値（ガウシアンのピークでの値）
            x_array : np.ndarray
                x方向の1次元array
            FWHM_ : float
                ガウシアンの半値全幅
            x_center : float
                ガウシアンがピークをとるときの x の値

            Returns
            -------
            np.ndarray
                横軸x_arrayに対するガウシアンの値の1次元array
            """

            y_gaussian = y_max * np.exp(
                - (x_array - x_center)**2 / (2 * (FWHM_ / (2 * np.sqrt(2 * np.log(2))))**2)
            )

            return y_gaussian

        def calc_tau_i_thermal(rambda_array_: np.ndarray) -> np.ndarray:
            """サーマルカットフィルターの透過率をcsvファイルから読み出して補完し、
            入力された波長に対して計算して透過率の1次元arrayとして出力

            Parameters
            ----------
            rambda_array_ : np.ndarray
                [m] 波長の1次元array

            Returns
            -------
            np.ndarray
                [無次元] サーマルカットフィルターの透過率の1次元array
            """
            # 神原M論 p94より、WG31050の透過率特性を
            # https://www.thorlabs.co.jp/NewGroupPage9.cfm?ObjectGroup_ID=3982
            # からダウンロードし、波長[um] と透過率[%] のみをtxtに保存
            i_thermal_filepath = "raw_data/tau_i_thermal.txt"

            # txtの読み出し
            raw = np.loadtxt(
                fname=i_thermal_filepath,
                delimiter=" ")

            # 読み出した波長[um]と透過率[%]を、波長[m]と透過率[無次元]に換算
            rambda_data_array_um, tau_i_thermal_data_array_percent = raw.T
            rambda_data_array = rambda_data_array_um * 1e-6
            tau_i_thermal_data_array = tau_i_thermal_data_array_percent * 1e-2

            # 透過率関数の作成
            tau_i_thermal_func = interpolate.interp1d(
                x=rambda_data_array,
                y=tau_i_thermal_data_array,
                kind="linear"
            )

            # light_instanceの波長に対応した透過率の1次元arrayの作成
            tau_i_thermal_array_ = tau_i_thermal_func(rambda_array_)

            return tau_i_thermal_array_

        def calc_S_photon_pix(
                light_instance_: LightGenenrator,
                Omega_pix_: float,
                A_GBT_: float,
                t_obs_: float,
                G_sys_: float) -> int:
            """光によるカウント値 S_photon_pixを計算

            oop観測見積もり.md
                └ 近赤外装置による撮像・分光 \n
                    └ Signalへの換算 \n

            Parameters
            ----------
            light_instance_ : LightGenenrator
                自作インスタンス,
                I_all_primeの状態になっている（装置透過率までかけてある）もの
            A_GBT_ : float
                [m^2] 望遠鏡の開口面積
            Omega_pix_ : float
                [sr / pix] 1pixelが見込む立体角
            t_obs_ : float
                [s] 積分時間
            G_sys_ : float
                [e- / DN]システムゲイン

            Returns
            -------
            int
                [DN / pix] 光による1pixelあたりのカウント値
            """

            # 入力パラメータ以外の文字の定義
            I_prime_all_ = light_instance_.get_I_prime()
            h_ = phys_consts.h
            c_ = phys_consts.c
            rambda_ = light_instance_.get_rambda()
            eta_ = 0.889  # [e- / photon] 検出器の量子効率

            # 波長積分での被積分関数の導出
            integrand_ = ((I_prime_all_ * A_GBT_ * Omega_pix_) / (h_ * (c_ / rambda_))) * eta_

            # 波長方向積分
            # 被積分関数の1次元arrayの各要素に対して波長方向の分割幅をかけてから総和をとる
            integration_result = np.sum(integrand_ * light_instance_.get_rambda_division_width())

            # S_photon.pixの導出
            S_photon_pix_ = integration_result * t_obs_ / G_sys_

            # カウント値は整数なので小数は切り捨て
            return np.floor(S_photon_pix_).astype(int)

        # バンドパスフィルタの定義
        tau_i_BPF = calc_gaussian(
            y_max=self.__tau_BPF_center,
            x_array=light_instance.get_rambda(),
            FWHM_=self.__FWHM_BPF,
            x_center=self.__rambda_BPF_center)

        # 装置透過率の導出
        if self.__is_TOPICS:
            # TOPICSの場合
            tau_i_lens = 0.9 ** 3
            tau_i_mirror = 1
            tau_i_ND = self.__tau_i_ND
            tau_i = tau_i_lens * tau_i_mirror * tau_i_BPF * tau_i_ND

        else:
            # ESPRITの場合
            tau_i_mirror = 0.96 ** 9
            tau_i_thermal = calc_tau_i_thermal(rambda_array_=light_instance.get_rambda())
            tau_i_ND = self.__tau_i_ND
            tau_i = tau_i_mirror * tau_i_BPF * self.__tau_i_ND * tau_i_thermal

        # 装置透過率を光にかける
        light_instance.multiply_I_prime_to(magnification=tau_i)

        # 光によるカウント値の計算
        S_photon_pix = calc_S_photon_pix(
            light_instance_=light_instance,
            Omega_pix_=self.get_Omega_pix(),
            A_GBT_=self.__GBT_instance.get_A_GBT(),
            t_obs_=t_obs,
            G_sys_=self.__G_sys)

        # 暗電流によるカウント値の計算（カウント値は整数）
        S_dark_pix = np.floor(self.__I_dark * t_obs / self.__G_sys).astype(int)

        # 読み出しノイズによるカウント値の計算（カウント値は整数）
        S_read_pix = np.floor((self.__N_e_read / self.__G_sys)**2).astype(int)

        # カウント値の総和を計算
        S_all_pix = S_photon_pix + S_dark_pix + S_read_pix

        # PDに蓄積された電荷とフルウェルの比を計算
        R_electron_FW = (S_photon_pix + S_dark_pix) / self.__S_FW_pix

        # fitsへの保存
        virtual_output_file_instance.set_S_all_pix(S_all_pix=S_all_pix)
        virtual_output_file_instance.set_t_obs(t_obs=t_obs)
        virtual_output_file_instance.set_n_bin_rambda(n_bin_rambda=self.__n_bin_rambda)
        virtual_output_file_instance.set_theta_pix(theta_pix=self.__theta_pix)
        virtual_output_file_instance.set_S_read_pix(S_read_pix=S_read_pix)
        virtual_output_file_instance.set_S_dark_pix(S_dark_pix=S_dark_pix)
        virtual_output_file_instance.set_R_electron_FW(R_electron_FW=R_electron_FW)


class SNRCalculator:

    def __init__(
            self,
            all_image_instance: VirtualOutputFileGenerator,
            sky_image_instance: VirtualOutputFileGenerator,
            n_bin_spatial: int) -> None:
        """観測対象のvirtual_output_fileとsky画像のvirtual_output_fileを保持し、
        任意のbinning数に対するSNRを計算

        oop観測見積もり.md
        └ SNRの計算 \n
            ├ 各Signalの導出 \n
            ├ binning後のSignalの導出 \n
            ├ 各Noiseの導出 \n
            └ SNRの導出 \n

        Parameters
        ----------
        all_image_instance : VirtualOutputFileGenerator
            観測対象とsky backgroundを両方含む virtual_output_file インスタンス
        sky_image_instance : VirtualOutputFileGenerator
            sky backgroundを両方含む virtual_output_file インスタンス
        n_bin_spatial : int
            [pix] 空間方向1辺のBinning数
        """

        # パラメータ設定
        self.__all_image_instance = all_image_instance
        self.__sky_image_instance = sky_image_instance
        self.__n_bin_spatial = n_bin_spatial

        # 空間・波長方向合わせた総binning数の計算
        self.__n_bin_total = self.__calc_n_bin_total(
            n_bin_spatial_=self.__n_bin_spatial,
            n_bin_rambda_=self.__all_image_instance.get_n_bin_rambda()
        )

        # SNRの計算
        S_all = self.__all_image_instance.get_S_all_pix() * self.__n_bin_total
        S_sky = self.__sky_image_instance.get_S_all_pix() * self.__n_bin_total
        self.__S_obj = S_all - S_sky
        self.__N_all = np.sqrt(S_all)
        self.__SNR = self.__S_obj / self.__N_all

        # 空間分解能の計算
        self.__spatial_resolution = 2 * self.__all_image_instance.get_theta_pix() * self.__n_bin_spatial

    def h(self):
        mkhelp(self)

    def get_n_bin_spatial(self):
        return self.__n_bin_spatial

    def get_n_bin_total(self):
        return self.__n_bin_total

    def get_S_obj(self):
        return self.__S_obj

    def get_N_all(self):
        return self.__N_all

    def get_SNR(self):
        return self.__SNR

    def get_spatial_resolution(self):
        return self.__spatial_resolution

    def __calc_n_bin_total(
            self,
            n_bin_spatial_: int,
            n_bin_rambda_: float) -> float:
        """撮像と分光で場合分けしてn_bin_totalを計算

        oop観測見積もり.md
            └ SNRの計算 \n
                └ binning後のSignalの導出 \n

        Parameters
        ----------
        n_bin_spatial_ : int
            [pix] 空間方向1辺のbinning数
        n_bin_rambda_ : float
            [pix] 波長方向のbinning数

        Returns
        -------
        float
            _description_
        """

        if n_bin_rambda_ == 1:
            # 波長方向のbinning = 1, つまり、撮像
            n_bin_total_ = n_bin_spatial_**2
            return n_bin_total_

        else:
            # 波長方向にbinningあり, つまり、分光
            n_bin_total_ = n_bin_spatial_ * n_bin_rambda_
            return n_bin_total_
