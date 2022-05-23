import numpy as np
import matplotlib.pyplot as plt


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


class PhysicalConstants:

    def __init__(self) -> None:
        self.c = 299792458  # [m/s] 光速
        self.h = 6.62606957e-34  # [J・s] プランク定数
        self.k_b = 1.3806488e-23  # [J/K] ボルツマン定数
        self.e = 1.60217663e-19  # [C/e-] 電気素量
        self.Q_T_dict = {  # [無次元] partition function
            100: 7.36, 200: 20.726, 300: 37.608, 400: 57.649, 500: 80.579,
            600: 106.393, 700: 135.33, 800: 167.812, 900: 204.388, 1000: 245.683, 1200: 345.179}

    def h(self):
        mkhelp(self)


class EmissionLineParameters:

    def __init__(
            self, lambda_um: float, g_ns: int, J_prime: int, A_if: float, E_prime: float) -> None:
        """__init__ 各輝線のパラメータ

        Parameters
        ----------
        lambda_um : float
            輝線の中心波長 [um]
        g_ns : int
            nuclear spin weight, 2 or 4 [無次元]
        J_prime : int
            回転準位 [無次元]
        A_if : float
            アインシュタインのA係数 [/s]
        E_prime : float
            energy of upper statement [/cm]
        """
        self.lambda_um = lambda_um
        self.omega_if = 1 / (self.lambda_um * 1e-6) * 1e-2  # 波数 [/cm]
        self.g_ns = g_ns
        self.J_prime = J_prime
        self.A_if = A_if
        self.E_prime = E_prime

    def h(self):
        mkhelp(self)
