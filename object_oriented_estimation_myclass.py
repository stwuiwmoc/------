import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import constants as phys_consts


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

    def show_rambda_vs_I_prime_plot(self) -> None:
        """横軸波長、縦軸分光放射輝度のグラフを表示
        """
        fig = plt.figure()
        gs = fig.add_gridspec(1, 1)

        ax = fig.add_subplot(gs[0, 0])

        ax.plot(
            self.__rambda,
            self.__I_prime)

        ax.grid()
        ax.set_xlabel("rambda [m]")
        ax.set_ylabel("I' [W / m^2 / sr / m]")

        fig.tight_layout()


class VirtualOutputFileGenerator:

    def __init__(self) -> None:
        """仮想的なFITSファイルとしてインスタンスを作成し、カウント値などを保持する

        oop観測見積もり.md
            └ 見積もりの概略 \n
                └ シミュレーション上でのFITSファイルの扱い \n

        """
        self.__S_all_pix = None
        self.__S_FW_pix = None
        self.__t_obs = None
        self.__n_bin_rambda = None

    def h(self):
        mkhelp(self)

    def get_S_all_pix(self) -> float:
        return self.__S_all_pix

    def get_S_FW_pix(self) -> float:
        return self.__S_FW_pix

    def get_t_obs(self) -> float:
        return self.__t_obs

    def get_n_bin_rambda(self) -> float:
        return self.__n_bin_rambda

    def set_S_all_pix(self, S_all_pix) -> None:
        self.__S_all_pix = S_all_pix

    def set_S_FW_pix(self, S_FW_pix) -> None:
        self.__S_FW_pix = S_FW_pix

    def set_t_obs(self, t_obs) -> None:
        self.__t_obs = t_obs

    def set_n_bin_rambda(self, n_bin_rambda) -> None:
        self.__n_bin_rambda = n_bin_rambda


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
                └ 望遠鏡の熱輻射による放射輝度 \n

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
                    └ 望遠鏡の熱輻射による放射輝度 \n

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
            rambda_fl_center: float,
            FWHM_fl: float,
            tau_fl_center: float,
            G_Amp: float,
            I_dark: float,
            N_read: float) -> None:
        """撮像装置のパラメータを保持

        oop観測見積もり.md
            └ 近赤外装置による撮像・分光 \n
                ├ 干渉フィルター透過率の導出 \n
                ├ 装置透過率の導出（TOPICS） \n
                ├ pixel数関連の導出 \n
                ├ システムゲインの導出 \n
                ├ 検出器に到達した段階での分光放射輝度 \n
                └ Signalへの換算 \n

        Parameters
        ----------
        rambda_fl_center : float
            [m] 干渉フィルターの中心波長
        FWHM_fl : float
            [m] 干渉フィルターの半値全幅
        tau_fl_center : float
            [無次元] 干渉フィルターの中心透過率
        G_Amp : float
            [無次元] プリアンプ基板の倍率
        I_dark : float
            [e- / s / pix] 検出器暗電流
        N_read : float
            [e-rms / pix] 駆動回路の読み出しノイズ
        """

        # --- 入力パラメータ・固定パラメータの代入 ---
        # 干渉フィルター透過率に関するパラメータ
        self.__rambda_fl_center = rambda_fl_center
        self.__FWHM_fl = FWHM_fl
        self.__tau_fl_center = tau_fl_center

        # ピクセル数関連のパラメータ
        self.__n_bin_rambda = 1  # [pix]

        # システムゲインに関連するパラメータ
        self.__G_Amp = G_Amp
        self.__G_sys = self.__calc_G_sys()
        self.__FW = 152000  # [e- / pix]
        self.__S_FW_pix = self.__FW / self.__G_sys

        # Signalへの換算で必要なパラメータ
        self.__I_dark = I_dark
        self.__N_read = N_read

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

    def get_rambda_fl_center(self) -> float:
        return self.__rambda_fl_center

    def get_FWHM_fl(self) -> float:
        return self.__FWHM_fl

    def get_tau_fl_center(self) -> float:
        return self.__tau_fl_center

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

    def get_N_read(self) -> float:
        return self.__N_read

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
            m_i_all = 2  # [無次元]
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
            virtual_output_file_instance: VirtualOutputFileGenerator) -> None:

        def calc_gaussian(
                y_max: float,
                x_array: np.ndarray,
                FWHM_: float,
                x_center: float) -> np.ndarray:
            """横軸x、縦軸y としてガウシアンを計算

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

        # 干渉フィルターの定義
        tau_i_filter = calc_gaussian(
            y_max=self.__tau_fl_center,
            x_array=light_instance.get_rambda(),
            FWHM_=self.__FWHM_fl,
            x_center=self.__rambda_fl_center)

        # テスト：フィルター透過率だけをlight_instanceにかける
        light_instance.multiply_I_prime_to(magnification=tau_i_filter)

        # カウント値への変換

        # fitsへの保存
