# -*- coding: utf-8 -*-
"""
Created on Tue May 29 10:24:49 2018

@author: msb
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import math
import os


def P_onedraw(w, start_conc, draw_size, back=False):

    Q3_s = start_conc[0]
    Q2_s = start_conc[1]
    Q1_s = start_conc[2]
    Q0_s = start_conc[3]

    p3 = Q3_s * w[0] / ((Q3_s * w[0]) + (Q2_s * w[1]) + (Q1_s * w[2]))
    p2 = Q2_s * w[1] / ((Q3_s * w[0]) + (Q2_s * w[1]) + (Q1_s * w[2]))
    p1 = Q1_s * w[2] / ((Q3_s * w[0]) + (Q2_s * w[1]) + (Q1_s * w[2]))

    p3 = p3 * draw_size
    p2 = p2 * draw_size
    p1 = p1 * draw_size

    if Q3_s - p3 < 0:
        next_Q3 = 0
    else:
        next_Q3 = Q3_s - p3

    if Q2_s + p3 - p2 < 0:
        next_Q2 = 0
    else:
        next_Q2 = Q2_s + p3 - p2

    if Q1_s + p2 - p1 < 0:
        next_Q1 = 0
    else:
        next_Q1 = Q1_s + p2 - p1

    if Q0_s + p1 < 0:
        next_Q0 = 0
    else:
        next_Q0 = Q0_s + p1

    return next_Q3, next_Q2, next_Q1, next_Q0


def P_draw(H1, tg, frac=None, s_plt=False, s_dat=False, p=False):
    """
       This function will plot the SRO scale structural evolution of silicate
       glasses by accounting for the enthalpic and entropic contributons to
       modifier-former interactions.

    =============================================================================
       model(H1, H2 = None, frac = None, s_plt = False, s_dat = False)
    =============================================================================

       where H1 is the necessary enthalpic contribution in a bunary glass.
       Examples are provided: "module.HNaSi", "module.HKSi", "module.HLiSi".

       H2 may be set to enthalpy values for a second modifier,
       where frac defines the fraction of the first to second modifier (0-1).

       s_plt and s_dat may be set to "True" to save the plot and data as png
       and csv files


       Example:

       >>> model(HNaSi, H2 = HLiSi, frac = 0.6, s_plt = True, s_dat = True)
    """
    draw_nr = list(range(300))
    draw_ar = np.array(draw_nr)

    M2O = []

    for i in draw_ar:
        next_mod = draw_ar[i] / (100 + draw_ar[i]) * 100
        M2O.append(next_mod)

    M2O.append(75)

    Tg = np.array(tg)

    w_Q3 = []
    w_Q2 = []
    w_Q1 = []

    if frac is None:

        H = [0, H1[0], H1[1]]

        for i in draw_ar:
            next_w_Q3 = math.exp(-H[0] / (Tg[i] * 0.00831))
            w_Q3.append(next_w_Q3)

            next_w_Q2 = math.exp(-H[1] / (Tg[i] * 0.00831))
            w_Q2.append(next_w_Q2)

            next_w_Q1 = math.exp(-H[2] / (Tg[i] * 0.00831))
            w_Q1.append(next_w_Q1)

        Q3 = [
            100,
        ]
        Q2 = [
            0,
        ]
        Q1 = [
            0,
        ]
        Q0 = [
            0,
        ]

        for i in draw_ar:

            p3 = (
                Q3[-1]
                * w_Q3[i]
                / ((Q3[-1] * w_Q3[i]) + (Q2[-1] * w_Q2[i]) +
                   (Q1[-1] * w_Q1[i]))
            )
            p2 = (
                Q2[-1]
                * w_Q2[i]
                / ((Q3[-1] * w_Q3[i]) + (Q2[-1] * w_Q2[i]) +
                   (Q1[-1] * w_Q1[i]))
            )
            p1 = (
                Q1[-1]
                * w_Q1[i]
                / ((Q3[-1] * w_Q3[i]) + (Q2[-1] * w_Q2[i]) +
                   (Q1[-1] * w_Q1[i]))
            )

            if Q3[-1] - p3 < 0:
                next_Q3 = 0
            else:
                next_Q3 = Q3[-1] - p3

            if Q2[-1] + p3 - p2 < 0:
                next_Q2 = 0
            else:
                next_Q2 = Q2[-1] + p3 - p2

            if Q1[-1] + p2 - p1 < 0:
                next_Q1 = 0
            else:
                next_Q1 = Q1[-1] + p2 - p1

            if Q0[-1] + p1 < 0:
                next_Q0 = 0
            else:
                next_Q0 = Q0[-1] + p1

            Q3.append(next_Q3)
            Q2.append(next_Q2)
            Q1.append(next_Q1)
            Q0.append(next_Q0)

    elif type(H1) is tuple:
        w_Na_Q3 = []
        w_Na_Q2 = []
        w_Na_Q1 = []

        w_K_Q3 = []
        w_K_Q2 = []
        w_K_Q1 = []

        for i in draw_ar:

            next_w_Na_Q3 = math.exp(-H1[1][0] / (Tg[i] * 0.00831))
            w_Na_Q3.append(next_w_Na_Q3)

            next_w_Na_Q2 = math.exp(-H1[1][1] / (Tg[i] * 0.00831))
            w_Na_Q2.append(next_w_Na_Q2)

            next_w_Na_Q1 = math.exp(-H1[1][2] / (Tg[i] * 0.00831))
            w_Na_Q1.append(next_w_Na_Q1)

            next_w_K_Q3 = math.exp(-H1[0][0] / (Tg[i] * 0.00831))
            w_K_Q3.append(next_w_K_Q3)

            next_w_K_Q2 = math.exp(-H1[0][1] / (Tg[i] * 0.00831))
            w_K_Q2.append(next_w_K_Q2)

            next_w_K_Q1 = math.exp(-H1[0][2] / (Tg[i] * 0.00831))
            w_K_Q1.append(next_w_K_Q1)

        Q3 = [
            100,
        ]
        Q2 = [
            0,
        ]
        Q1 = [
            0,
        ]
        Q0 = [
            0,
        ]

        for i in draw_ar:

            p3 = (
                (
                    Q3[-1]
                    * w_K_Q3[i]
                    / (
                        (Q3[-1] * w_K_Q3[i])
                        + (Q2[-1] * w_K_Q2[i])
                        + (Q1[-1] * w_K_Q1[i])
                    )
                )
                * frac[0]
            ) + (
                (
                    Q3[-1]
                    * w_Na_Q3[i]
                    / (
                        (Q3[-1] * w_Na_Q3[i])
                        + (Q2[-1] * w_Na_Q2[i])
                        + (Q1[-1] * w_Na_Q1[i])
                    )
                )
                * frac[1]
            )
            p2 = (
                (
                    Q2[-1]
                    * w_K_Q2[i]
                    / (
                        (Q3[-1] * w_K_Q3[i])
                        + (Q2[-1] * w_K_Q2[i])
                        + (Q1[-1] * w_K_Q1[i])
                    )
                )
                * frac[0]
            ) + (
                (
                    Q2[-1]
                    * w_Na_Q2[i]
                    / (
                        (Q3[-1] * w_Na_Q3[i])
                        + (Q2[-1] * w_Na_Q2[i])
                        + (Q1[-1] * w_Na_Q1[i])
                    )
                )
                * frac[1]
            )
            p1 = (
                (
                    Q1[-1]
                    * w_K_Q1[i]
                    / (
                        (Q3[-1] * w_K_Q3[i])
                        + (Q2[-1] * w_K_Q2[i])
                        + (Q1[-1] * w_K_Q1[i])
                    )
                )
                * frac[0]
            ) + (
                (
                    Q1[-1]
                    * w_Na_Q1[i]
                    / (
                        (Q3[-1] * w_Na_Q3[i])
                        + (Q2[-1] * w_Na_Q2[i])
                        + (Q1[-1] * w_Na_Q1[i])
                    )
                )
                * frac[1]
            )

            if Q3[-1] - p3 < 0:
                next_Q3 = 0
            else:
                next_Q3 = Q3[-1] - p3

            if Q2[-1] + p3 - p2 < 0:
                next_Q2 = 0
            else:
                next_Q2 = Q2[-1] + p3 - p2

            if Q1[-1] + p2 - p1 < 0:
                next_Q1 = 0
            else:
                next_Q1 = Q1[-1] + p2 - p1

            if Q0[-1] + p1 < 0:
                next_Q0 = 0
            else:
                next_Q0 = Q0[-1] + p1

            Q3.append(next_Q3)
            Q2.append(next_Q2)
            Q1.append(next_Q1)
            Q0.append(next_Q0)

    else:
        return print("Wrong H format")
    if s_plt is False and p is True:
        plt.plot(
            M2O,
            Q3,
            "r-",
            M2O,
            Q2,
            "k-",
            M2O,
            Q1,
            "b-",
            M2O,
            Q0,
            "g-",
        )
        plt.axis([0, 75, 0, 100])
        plt.legend(["$Q^3$", "$Q^2$", "$Q^1$", "$Q^0$"])
        plt.xlabel("Modifier mol %")
        plt.ylabel("Qn species concentration")
        plt.title("Qn distribution")
        plt.show()

    if s_plt is True:
        if not os.path.exists("P2O5_Structure"):
            os.mkdir("P2O5_Structure")
        plt.plot(
            M2O,
            Q3,
            "r-",
            M2O,
            Q2,
            "k-",
            M2O,
            Q1,
            "b-",
            M2O,
            Q0,
            "g-",
        )
        plt.axis([0, 75, 0, 100])
        plt.legend(["$Q^3$", "$Q^2$", "$Q^1$", "$Q^0$"])
        plt.xlabel("Modifier mol %")
        plt.ylabel("Qn species concentration")
        plt.title("Qn distribution")
        plt.savefig(os.path.join("P2O5_Structure", "Qn_distribution.png"))
        plt.show()

    if s_dat is True:
        if not os.path.exists("P2O5_Structure"):
            os.mkdir("P2O5_Structure")
        m_data = np.column_stack([M2O, Q3, Q2, Q1, Q0])
        np.savetxt(os.path.join("P2O5_Structure", "Model_data.csv"), m_data)

    elif s_plt is False and s_dat is False and p is False:
        return M2O, Q3, Q2, Q1, Q0


def P_SSE(H1, data, tg, frac=None, s_plt=False, s_dat=False, p=False):
    """
       This function will plot the SRO scale structural evolution of silicate
       glasses by accounting for the enthalpic and entropic contributons to
       modifier-former interactions.

    =============================================================================
       model(H1, H2 = None, frac = None, s_plt = False, s_dat = False)
    =============================================================================

       where H1 is the necessary enthalpic contribution in a bunary glass.
       Examples are provided: "module.HNaSi", "module.HKSi", "module.HLiSi".

       H2 may be set to enthalpy values for a second modifier, where frac
       defines the fraction of the first to second modifier (0-1).

       s_plt and s_dat may be set to "True" to save the plot and data as png
       and csv files


       Example:

       >>> model(HNaSi, H2 = HLiSi, frac = 0.6, s_plt = True, s_dat = True)
    """

    mod_data = data[0]
    Q3_data = data[1]
    Q2_data = data[2]
    Q1_data = data[3]
    Q0_data = data[4]

    M2O, Q3, Q2, Q1, Q0 = P_draw(H1, tg, frac)

    if s_plt is False and p is True:
        plt.plot(
            M2O,
            Q3,
            "r-",
            M2O,
            Q2,
            "k-",
            M2O,
            Q1,
            "b-",
            M2O,
            Q0,
            "g-",
        )
        plt.plot(
            mod_data,
            Q3_data,
            "rd",
            mod_data,
            Q2_data,
            "kd",
            mod_data,
            Q1_data,
            "bd",
            mod_data,
            Q0_data,
            "gd",
        )
        plt.axis([0, 75, 0, 100])
        plt.legend(["$Q^3$", "$Q^2$", "$Q^1$", "$Q^0$"])
        plt.xlabel("Modifier mol %")
        plt.ylabel("Qn species concentration")
        plt.title("Qn distribution")
        plt.show()

    if s_plt is True:
        if not os.path.exists("P2O5_Structure"):
            os.mkdir("P2O5_Structure")
        plt.plot(
            M2O,
            Q3,
            "r-",
            M2O,
            Q2,
            "k-",
            M2O,
            Q1,
            "b-",
            M2O,
            Q0,
            "g-",
        )
        plt.plot(
            mod_data,
            Q3_data,
            "rd",
            mod_data,
            Q2_data,
            "kd",
            mod_data,
            Q1_data,
            "bd",
            mod_data,
            Q0_data,
            "gd",
        )
        plt.axis([0, 75, 0, 100])
        plt.legend(["$Q^3$", "$Q^2$", "$Q^1$", "$Q^0$"])
        plt.xlabel("Modifier mol %")
        plt.ylabel("Qn species concentration")
        plt.title("Qn distribution")
        plt.savefig(os.path.join("P2O5_Structure", "Qn_distribution.png"))
        plt.show()

    if s_dat is True:
        if not os.path.exists("P2O5_Structure"):
            os.mkdir("P2O5_Structure")
        m_data = np.column_stack([M2O, Q3, Q2, Q1, Q0])
        np.savetxt(os.path.join("P2O5_Structure", "Model_data.csv"), m_data)

    if p is False:
        mod_m = []
        Q3_m = []
        Q2_m = []
        Q1_m = []
        Q0_m = []

        for i in mod_data:
            next_mod_m = min(M2O, key=lambda x: abs(x - i))
            mod_m.append(next_mod_m)

        for i in mod_m:
            ind = M2O.index(i)
            next_Q3_m = Q3[ind]
            Q3_m.append(next_Q3_m)

            next_Q2_m = Q2[ind]
            Q2_m.append(next_Q2_m)

            next_Q1_m = Q1[ind]
            Q1_m.append(next_Q1_m)

            next_Q0_m = Q0[ind]
            Q0_m.append(next_Q0_m)

        Q3_m = np.array(Q3_m)
        Q2_m = np.array(Q2_m)
        Q1_m = np.array(Q1_m)
        Q0_m = np.array(Q0_m)

        SSE = sum(
            ((Q3_data - Q3_m) ** 2)
            + ((Q2_data - Q2_m) ** 2)
            + ((Q1_data - Q1_m) ** 2)
            + ((Q0_data - Q0_m) ** 2)
        )

        return SSE


def P_engine(fil, data, tg, it=10):
    dat = data
    w0 = [20, 30]

    minimizer_kwargs = {
        "method": "COBYLA",
        "args": (
            dat,
            tg,
        ),
    }
    res = scipy.optimize.basinhopping(
        P_SSE,
        w0,
        niter=it,
        T=2.0,
        stepsize=1,
        minimizer_kwargs=minimizer_kwargs,
        take_step=None,
        accept_test=None,
        callback=None,
        interval=50,
        disp=True,
        niter_success=None,
        seed=None,
    )

    return res.x
