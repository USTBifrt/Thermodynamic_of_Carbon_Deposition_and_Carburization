import sympy as sp
import numpy as np
from PhaseD_logKf import *
import pandas as pd
import plotly.figure_factory as ff

Ptot = 1
inputH2O = 0
inputCO2 = 0

T = 1173
# for T in range (773,1373,100):
# for Ptot in np.arange (1.5,8.5,0.5):
for inputH2O in np.arange (0.1,0.6,0.1):
    V = 1-inputH2O-inputCO2
# # CO + 3H2 = CH4 +H2O  …  （1）
# logK1 = logKf("CO2", T) - 2 *  logKf("CO", T)
# k1 = 10 ** logK1
# # CO + H2O = H2 +CO2   …  （2）
# logK2 = logKf("H2O", T) - logKf("CO", T) - logKf("H2", T)
# k2 = 10 ** logK2
# # C+CO2 = 2CO          …  （3）
# logK3 = 2 * logKf("H2", T) - logKf("CH4", T)
# k3 = 10 ** logK3
# print(k1,k2,k3)

    # 2CO = C + CO2
    deltG1 = -166500 + 171 * T
    k1 = math.exp(- deltG1 / (8.314 * T))

    # H2 + CO = C + H2O
    deltG2 = -133100 + 141.65 * T
    k2 = math.exp(- deltG2 / (8.314 * T))

    # CH4 = C + 2H2
    deltG3 = 91044 - 110.67 * T
    k3 = math.exp(- deltG3 / (8.314 * T))
    print(k1, k2, k3)

    nCO_array = []
    nH2_array = []
    nCH4_array = []
    solution_C_array = []

    for nCO in np.arange(0.1, 1, 0.1):
        nCO = round(nCO, 1)
        for nH2 in np.arange(0.1, round(1 - nCO, 1), 0.1):
            nH2 = round(nH2, 1)
            nCH4 = abs(round(1 - nCO - nH2, 1))
            iCO = nCO * V
            iH2 = nH2 * V
            iCH4 = nCH4 * V
            XC = iCO + iCH4 + inputCO2
            XH = 2 * iH2 + 4*iCH4 + 2 * inputH2O
            XO = iCO + 2 * inputCO2 + inputH2O
            found_solution = False

            for in_xCO in np.arange(0, XC, 0.05):
                for in_xC in np.arange(0, XC-in_xCO, 0.05):
                    for in_xCO2 in np.arange(0, XC - in_xCO - in_xC, 0.05):
                        in_xCH4 = XC - in_xCO - in_xCO2 - in_xC
                        in_xH2O = XO - in_xCO - 2 * in_xCO2
                        in_xH2 = (XH - 4 * in_xCH4 - 2 * in_xH2O) / 2

                        initial_guess = [in_xCO, in_xH2, in_xCH4, in_xCO2, in_xH2O, in_xC]
                        # print(initial_guess)
                        if found_solution:
                            break  # 如果已找到解，退出外层循环
                        # elif in_xCO == nCO:
                        #     print("%d, %f, %f, %f: 没有找到解" % (T, nCO, nH2, nCH4))
                        #     break
                        if all(element >= 0 if not isinstance(element, sp.Expr) else True for element in initial_guess):
                            # print(initial_guess)

                            xCO, xH2, xCH4, xCO2, xH2O, xC = sp.symbols('xCO xH2 xCH4 xCO2 xH2O xC')

                            PCO = Ptot * xCO/(xCO + xH2 + xCH4 + xCO2 + xH2O)
                            PCO2 = Ptot * xCO2 / (xCO + xH2 + xCH4 + xCO2 + xH2O)
                            PH2O = Ptot * xH2O / (xCO + xH2 + xCH4 + xCO2 + xH2O)
                            PH2 = Ptot * xH2 / (xCO + xH2 + xCH4 + xCO2 + xH2O)
                            PCH4 = Ptot * xCH4 / (xCO + xH2 + xCH4 + xCO2 + xH2O)

                            # 方程组
                            eq1 = k1 - (PCO2) / (PCO ** 2)
                            eq2 = k2 - (PH2O) / (PH2 * PCO)
                            eq3 = k3 - (PH2) ** 2 / (PCH4)

                            eq4 = xCO + xCO2 + xCH4 + xC - XC
                            eq5 = xCO + 2 * xCO2 + xH2O - XO
                            eq6 = 2 * xH2 + 2 * xH2O + 4 * xCH4 - XH

                            vars = [xCO, xH2, xCH4, xCO2, xH2O, xC]
                            equations = [eq1, eq2, eq3, eq4, eq5, eq6]

                            try:
                                solution = sp.nsolve(equations, vars, initial_guess)

                                if all(element >= 0 for element in solution):
                                    print("%s, %f, %f, %f" % (inputH2O, nCO, nH2, nCH4), solution[0:])
                                    solution_C_array.append(float(solution[5]))
                                    nCO_array.append(nCO)
                                    nH2_array.append(nH2)
                                    nCH4_array.append(nCH4)
                                    # print("Ptot:%f"%(nCO + nH2 + nCH4 - xH2 + xCH4))
                                    # 设置标志，表示已找到解
                                    found_solution = True


                            except ValueError as e:
                                pass
                                # print(f"Error: {e}. Skipping initial_guess.")

                            except ZeroDivisionError as e:
                                pass
                                # print(f"Error: {e}. Skipping initial_guess.")

    # solution_C_array = np.array(solution_C_array)
    # print(solution_C_array)
    # # nCO_array = [0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    # # nCO_array = np.array(nCO_array)
    # # nH2_array = np.array(nH2_array)
    # # nCH4_array = np.array(nCH4_array)
    #
    # fig = ff.create_ternary_contour(np.array([nCO_array, nH2_array, nCH4_array]), solution_C_array,
    #                                 pole_labels=['CO', 'H2', 'CH4'],
    #                                 interp_mode='cartesian',
    #                                 coloring=None,
    #                                 ncontours=12,
    #                                 colorscale="Jet",
    #                                 showmarkers=True,
    #                                 showscale=True,
    #                                 range=dict(A=(0.1, 0.8), B=(0.1, 0.8), C=(0.1, 0.8)))
    # fig.write_image("E://Garbage//Desktop//Phd.FeOPhase//refs - Iron oxides Reduction//Phase D//fig_%d.png" % (Ptot), scale=5)


Ptot = 1
inputH2O = 0
inputCO2 = 0

T = 1173
# for T in range (773,1373,100):
# for Ptot in np.arange (1.5,8.5,0.5):
for inputCO2 in np.arange (0.1,0.6,0.1):
# # CO + 3H2 = CH4 +H2O  …  （1）
# logK1 = logKf("CO2", T) - 2 *  logKf("CO", T)
# k1 = 10 ** logK1
# # CO + H2O = H2 +CO2   …  （2）
# logK2 = logKf("H2O", T) - logKf("CO", T) - logKf("H2", T)
# k2 = 10 ** logK2
# # C+CO2 = 2CO          …  （3）
# logK3 = 2 * logKf("H2", T) - logKf("CH4", T)
# k3 = 10 ** logK3
# print(k1,k2,k3)
    V = 1-inputH2O-inputCO2
    # 2CO = C + CO2
    deltG1 = -166500 + 171 * T
    k1 = math.exp(- deltG1 / (8.314 * T))

    # H2 + CO = C + H2O
    deltG2 = -133100 + 141.65 * T
    k2 = math.exp(- deltG2 / (8.314 * T))

    # CH4 = C + 2H2
    deltG3 = 91044 - 110.67 * T
    k3 = math.exp(- deltG3 / (8.314 * T))
    print(k1, k2, k3)

    nCO_array = []
    nH2_array = []
    nCH4_array = []
    solution_C_array = []

    for nCO in np.arange(0.1, 1, 0.1):
        nCO = round(nCO, 1)
        for nH2 in np.arange(0.1, round(1 - nCO, 1), 0.1):
            nH2 = round(nH2, 1)
            nCH4 = abs(round(1 - nCO - nH2, 1))
            iCO = nCO * V
            iH2 = nH2 * V
            iCH4 = nCH4 * V
            XC = iCO + iCH4 + inputCO2
            XH = 2 * iH2 + 4*iCH4 + 2 * inputH2O
            XO = iCO + 2 * inputCO2 + inputH2O
            found_solution = False

            for in_xCO in np.arange(0, XC, 0.05):
                for in_xC in np.arange(0, XC-in_xCO, 0.05):
                    for in_xCO2 in np.arange(0, XC - in_xCO - in_xC, 0.05):
                        in_xCH4 = XC - in_xCO - in_xCO2 - in_xC
                        in_xH2O = XO - in_xCO - 2 * in_xCO2
                        in_xH2 = (XH - 4 * in_xCH4 - 2 * in_xH2O) / 2

                        initial_guess = [in_xCO, in_xH2, in_xCH4, in_xCO2, in_xH2O, in_xC]
                        # print(initial_guess)
                        if found_solution:
                            break  # 如果已找到解，退出外层循环
                        # elif in_xCO == nCO:
                        #     print("%d, %f, %f, %f: 没有找到解" % (T, nCO, nH2, nCH4))
                        #     break
                        if all(element >= 0 if not isinstance(element, sp.Expr) else True for element in initial_guess):
                            # print(initial_guess)

                            xCO, xH2, xCH4, xCO2, xH2O, xC = sp.symbols('xCO xH2 xCH4 xCO2 xH2O xC')

                            PCO = Ptot * xCO/(xCO + xH2 + xCH4 + xCO2 + xH2O)
                            PCO2 = Ptot * xCO2 / (xCO + xH2 + xCH4 + xCO2 + xH2O)
                            PH2O = Ptot * xH2O / (xCO + xH2 + xCH4 + xCO2 + xH2O)
                            PH2 = Ptot * xH2 / (xCO + xH2 + xCH4 + xCO2 + xH2O)
                            PCH4 = Ptot * xCH4 / (xCO + xH2 + xCH4 + xCO2 + xH2O)

                            # 方程组
                            eq1 = k1 - (PCO2) / (PCO ** 2)
                            eq2 = k2 - (PH2O) / (PH2 * PCO)
                            eq3 = k3 - (PH2) ** 2 / (PCH4)

                            eq4 = xCO + xCO2 + xCH4 + xC - XC
                            eq5 = xCO + 2 * xCO2 + xH2O - XO
                            eq6 = 2 * xH2 + 2 * xH2O + 4 * xCH4 - XH

                            vars = [xCO, xH2, xCH4, xCO2, xH2O, xC]
                            equations = [eq1, eq2, eq3, eq4, eq5, eq6]

                            try:
                                solution = sp.nsolve(equations, vars, initial_guess)

                                if all(element >= 0 for element in solution):
                                    print("%s, %f, %f, %f" % (inputCO2, nCO, nH2, nCH4), solution[0:])
                                    solution_C_array.append(float(solution[5]))
                                    nCO_array.append(nCO)
                                    nH2_array.append(nH2)
                                    nCH4_array.append(nCH4)
                                    # print("Ptot:%f"%(nCO + nH2 + nCH4 - xH2 + xCH4))
                                    # 设置标志，表示已找到解
                                    found_solution = True


                            except ValueError as e:
                                pass
                                # print(f"Error: {e}. Skipping initial_guess.")

                            except ZeroDivisionError as e:
                                pass
                                # print(f"Error: {e}. Skipping initial_guess.")
