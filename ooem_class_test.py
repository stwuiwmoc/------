# %%
import importlib

import object_oriented_estimation_myclass as ooem


if __name__ == "__main__":
    importlib.reload(ooem)

    # グローバル変数の定義
    column_density_H3plus = 2.0e+16  # [/m^2] H3+カラム密度
    T_thermospheric_H3plus = 1200  # [K] H3+熱圏温度

    # 各インスタンス生成
    light = ooem.LightGenenrator(
        rambda_division_width=0.1e-9,
        rambda_lower_limit=3.3e-6,
        rambda_upper_limit=3.5e-6)

    R_3_0 = ooem.H3plusAuroralEmission(
        rambda_obj=3.4128e-6,
        N_H3p=column_density_H3plus,
        g_ns=4,
        J_prime=4,
        A_if=177.6,
        E_prime=3382.9299,
        T_hypo=T_thermospheric_H3plus)

    R_3_1 = ooem.H3plusAuroralEmission(
        rambda_obj=3.4149e-6,
        N_H3p=column_density_H3plus,
        g_ns=2,
        J_prime=4,
        A_if=110.4,
        E_prime=3359.002,
        T_hypo=T_thermospheric_H3plus)

    R_3_2 = ooem.H3plusAuroralEmission(
        rambda_obj=3.4207e-6,
        N_H3p=column_density_H3plus,
        g_ns=2,
        J_prime=4,
        A_if=86.91,
        E_prime=3287.2629,
        T_hypo=T_thermospheric_H3plus)

    R_3_3 = ooem.H3plusAuroralEmission(
        rambda_obj=3.4270e-6,
        N_H3p=column_density_H3plus,
        g_ns=4,
        J_prime=4,
        A_if=42.94,
        E_prime=3169.252,
        T_hypo=T_thermospheric_H3plus)

    R_4_3 = ooem.H3plusAuroralEmission(
        rambda_obj=3.4547e-6,
        N_H3p=column_density_H3plus,
        g_ns=4,
        J_prime=5,
        A_if=62.91,
        E_prime=3489.2151,
        T_hypo=T_thermospheric_H3plus)

    R_4_4 = ooem.H3plusAuroralEmission(
        rambda_obj=3.4548e-6,
        N_H3p=column_density_H3plus,
        g_ns=2,
        J_prime=5,
        A_if=122.9,
        E_prime=3332.4121,
        T_hypo=T_thermospheric_H3plus)

    T60 = ooem.GroundBasedTelescope(
        D_GBT=0.6,
        FNO_GBT=12,
        T_GBT=280,
        tau_GBT=0.66)

    # GBTのテスト
    # プランクの法則の確認用
    light_test1 = ooem.LightGenenrator(
        rambda_lower_limit=1e-6,
        rambda_upper_limit=20e-6,
        rambda_division_width=1e-9)

    T60.pass_through(light_instance=light_test1)
    light_test1.show_rambda_vs_I_prime_plot()

    # 絶対値の概算用
    light_test2 = ooem.LightGenenrator(
        rambda_lower_limit=3.4128e-6,
        rambda_upper_limit=3.4128e-6 + 0.1e-9,
        rambda_division_width=0.1e-9)

    T60.pass_through(light_instance=light_test2)
    print(light_test2.get_I_prime() * 17e-9)
    print(light_test2.get_I_prime() * 17e-9 * T60.get_tau_GBT())  # calc_S_xxの計算の中で間違ってI_GBTにも望遠鏡透過率を掛けて過小評価してる

    # 輝線発光を加える
    light.show_rambda_vs_I_prime_plot()

    R_3_0.add_auroral_emission_to(light_instance=light)
    light.show_rambda_vs_I_prime_plot()
    print(light.get_I())

    R_3_1.add_auroral_emission_to(light_instance=light)
    light.show_rambda_vs_I_prime_plot()

    R_3_2.add_auroral_emission_to(light_instance=light)
    light.show_rambda_vs_I_prime_plot()

    R_3_3.add_auroral_emission_to(light_instance=light)
    light.show_rambda_vs_I_prime_plot()

    R_4_3.add_auroral_emission_to(light_instance=light)
    light.show_rambda_vs_I_prime_plot()

    R_4_4.add_auroral_emission_to(light_instance=light)
    light.show_rambda_vs_I_prime_plot()

    T60.pass_through(light_instance=light)
    light.show_rambda_vs_I_prime_plot()
