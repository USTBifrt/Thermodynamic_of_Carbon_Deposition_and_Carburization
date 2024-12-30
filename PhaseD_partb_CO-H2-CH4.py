import sympy as sp
import numpy as np
from PhaseD_logKf import *
import pandas as pd

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

    for nCO in np.arange(0.1,0.6,0.1):
        for nH2 in np.arange(0.8,1-nCO,0.1):
            n = 1
            nCH4 = 1-nCO-nH2
            XC = nCO + nCH4
            XH = 2 * nH2 + 4*nCH4
            XO = nCO
            found_solution = False
            for in_xCO in np.arange(0,nCO,0.01):
                for in_xH2 in np.arange(0, nH2, 0.01):
                    for in_xCH4 in np.arange(-nCH4, nCH4, 0.01):
                        if n != 1:
                            break  # 如果已找到解，退出外层循环
                        if found_solution:
                            break  # 如果已找到解，退出外层循环
                        elif in_xCO == nCO-0.01:
                            print("%d, %f, %f, %f: 没有找到解" % (T, nCO, nH2, nCH4))
                            break

                        initial_guess = [in_xCO,in_xH2,in_xCH4]
                        # print(initial_guess)

                        xCO, xH2, xCH4 = sp.symbols('xCO xH2 xCH4')
                        Ptot = nCO + nH2 + nCH4 - xH2 + xCH4

                        # 方程组
                        eq1 = k1 - (xCO/Ptot) / (((nCO-xCO-xH2)/Ptot) ** 2)
                        eq2 = k2 - (xH2/Ptot) / (((nH2-xH2+2*xCH4)/Ptot) * ((nCO-xCO-xH2)/Ptot))
                        eq3 = k3 - ((nH2-xH2+2*xCH4)/Ptot) ** 2 / ((nCH4-xCH4)/Ptot)

                        vars = [xCO, xH2, xCH4]
                        equations = [eq1, eq2, eq3]

                        try:
                            solution = sp.nsolve(equations, vars, initial_guess)
                            # print(solution)

                            if all(element >= 0 for element in solution):
                                print("%d, %f, %f, %f" % (T, nCO, nH2, nCH4), solution[0:])

                                # print("Ptot:%f"%(nCO + nH2 + nCH4 - xH2 + xCH4))
                                # 设置标志，表示已找到解
                                found_solution = True


                        except ValueError as e:
                            pass
                            # print(f"Error: {e}. Skipping initial_guess.")

                        except ZeroDivisionError as e:
                            pass
                            # print(f"Error: {e}. Skipping initial_guess.")
