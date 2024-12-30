import sympy as sp
import numpy as np
from PhaseD_logKf import *
import pandas as pd

n = 0
for T in range (873,1373,100):

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

    for XO in np.arange(0.05,0.55,0.05):
        in_XO = XO
        in_XC = XO
        in_XH = 1-in_XC-in_XO
        nC = 0

        solutions = []
        initial_guesses = []
        found_solution = False
        for in_nCO in np.arange(0,in_XC,0.0002):
            for in_nCO2 in np.arange(0,in_XC,0.0002):
                if found_solution:
                    break  # 如果已找到解，退出外层循环
                in_nCH4 = in_XC -in_nCO - in_nCO2
                in_nH2O = in_XO - in_nCO - 2*in_nCO2
                in_nH2 = (in_XH - 4*in_nCH4 - 2*in_nH2O)/2
                initial_guess = [in_nCO, in_nCO2, in_nH2O, in_nH2, in_nCH4, in_XC, in_XH]
                # 检查所有元素是否大于0
                if all(element > 0 if not isinstance(element, sp.Expr) else True for element in initial_guess):
                    # print("initial_guess:", initial_guess)

                    nCO, nCO2, nH2O, nH2, nCH4, XC, XH = sp.symbols('nCO nCO2 nH2O nH2 nCH4 XC XH')
                    # 方程组
                    eq1 = k1 - (nCO2) / (nCO ** 2)
                    eq2 = k2 - (nH2O) / (nH2 * nCO)
                    eq3 = k3 - nH2 ** 2 / nCH4

                    # 其他方程
                    eq4 = nCO + nCO2 + nCH4 + nC - XC
                    eq5 = nCO + 2 * nCO2 + nH2O - XO
                    eq6 = 2 * nH2 + 2 * nH2O + 4 * nCH4 - XH
                    eq7 = XO + XC + XH - 1

                    vars = [nCO, nCO2, nH2O, nH2, nCH4, XC, XH]
                    equations = [eq1, eq2, eq3, eq4, eq5, eq6, eq7]

                    try:
                        # Attempt to find a root
                        solution = sp.nsolve(equations, vars, initial_guess)
                        # print("Solution:", solution)

                        # Check if all elements of the solution are greater than 0
                        if all(element > 0 for element in solution):
                            # print("XC:", solution[5])
                            # solutions.append(solution)
                            n += 1
                            # 设置标志，表示已找到解
                            found_solution = True


                    except ValueError as e:
                        pass
                        # print(f"Error: {e}. Skipping initial_guess.")

                    except ZeroDivisionError as e:
                        pass
                        # print(f"Error: {e}. Skipping initial_guess.")
        print("%d, %f,"%(T,XO),solution[0:])

