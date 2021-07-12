# -*- coding: utf-8 -*-
"""
Created on Tue May 29 10:24:49 2018

@author: msb
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import scipy.optimize


def Si_onedraw(w, start_conc, draw_size, back=False):

    Q4_s = start_conc[0]
    Q3_s = start_conc[1]
    Q2_s = start_conc[2]
    Q1_s = start_conc[3]
    Q0_s = start_conc[4]

    if back:

        wQ0 = 1 / w[3]
        wQ1 = 1 / w[2]
        wQ2 = 1 / w[1]
        wQ3 = 1 / w[0]

        p0 = Q0_s * wQ0 / ((Q0_s * wQ0) + (Q1_s * wQ1) +
                           (Q2_s * wQ2) + (Q3_s * wQ3))
        p1 = Q1_s * wQ1 / ((Q0_s * wQ0) + (Q1_s * wQ1) +
                           (Q2_s * wQ2) + (Q3_s * wQ3))
        p2 = Q2_s * wQ2 / ((Q0_s * wQ0) + (Q1_s * wQ1) +
                           (Q2_s * wQ2) + (Q3_s * wQ3))
        p3 = Q3_s * wQ3 / ((Q0_s * wQ0) + (Q1_s * wQ1) +
                           (Q2_s * wQ2) + (Q3_s * wQ3))

        p0 = p0 * draw_size
        p3 = p3 * draw_size
        p2 = p2 * draw_size
        p1 = p1 * draw_size

        if Q0_s + p0 < 0:
            next_Q0 = 0
            p1 = p1 + Q0_s + p0
        else:
            next_Q0 = Q0_s + p0

        if Q1_s - p0 + p1 < 0:
            next_Q1 = 0
            p2 = p2 + Q1_s + p1
        else:
            next_Q1 = Q1_s - p0 + p1

        if Q2_s - p1 + p2 < 0:
            next_Q2 = 0
            p3 = p3 + Q2_s + p2
        else:
            next_Q2 = Q2_s - p1 + p2

        if Q3_s - p2 + p3 < 0:
            next_Q3 = 0
            p3 = p3 + Q3_s + p3
        else:
            next_Q3 = Q3_s - p2 + p3

        if Q4_s - p3 < 0:
            next_Q4 = 0
        else:
            next_Q4 = Q4_s - p3

    else:

        p4 = (
            Q4_s
            * w[0]
            / ((Q4_s * w[0]) + (Q3_s * w[1]) + (Q2_s * w[2]) + (Q1_s * w[3]))
        )
        p3 = (
            Q3_s
            * w[1]
            / ((Q4_s * w[0]) + (Q3_s * w[1]) + (Q2_s * w[2]) + (Q1_s * w[3]))
        )
        p2 = (
            Q2_s
            * w[2]
            / ((Q4_s * w[0]) + (Q3_s * w[1]) + (Q2_s * w[2]) + (Q1_s * w[3]))
        )
        p1 = (
            Q1_s
            * w[3]
            / ((Q4_s * w[0]) + (Q3_s * w[1]) + (Q2_s * w[2]) + (Q1_s * w[3]))
        )

        p4 = p4 * draw_size
        p3 = p3 * draw_size
        p2 = p2 * draw_size
        p1 = p1 * draw_size

        if Q4_s - p4 < 0:
            next_Q4 = 0
        else:
            next_Q4 = Q4_s - p4

        if Q3_s + p4 - p3 < 0:
            next_Q3 = 0
        else:
            next_Q3 = Q3_s + p4 - p3

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

    return next_Q4, next_Q3, next_Q2, next_Q1, next_Q0


def Si_draw(H1, tg, frac=None, s_plt=False, s_dat=False, p=False):
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
    draw_nr = list(range(400))
    draw_ar = np.array(draw_nr)

    M2O = []

    for i in draw_ar:
        next_mod = draw_ar[i] * 0.5 / (100 + draw_ar[i] * 0.5) * 100
        M2O.append(next_mod)

    M2O.append(67)

    Tg = np.array(tg)

    w_Q4 = []
    w_Q3 = []
    w_Q2 = []
    w_Q1 = []

    if frac is None:

        H = [0, H1[0], H1[1], H1[2]]

        for i in draw_ar:
            next_w_Q4 = math.exp(-H[0] / (Tg[i] * 0.00831))
            w_Q4.append(next_w_Q4)

            next_w_Q3 = math.exp(-H[1] / (Tg[i] * 0.00831))
            w_Q3.append(next_w_Q3)

            next_w_Q2 = math.exp(-H[2] / (Tg[i] * 0.00831))
            w_Q2.append(next_w_Q2)

            next_w_Q1 = math.exp(-H[3] / (Tg[i] * 0.00831))
            w_Q1.append(next_w_Q1)

        Q4 = [
            100,
        ]
        Q3 = [
            0,
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

            p4 = (
                Q4[-1]
                * w_Q4[i]
                / (
                    (Q4[-1] * w_Q4[i])
                    + (Q3[-1] * w_Q3[i])
                    + (Q2[-1] * w_Q2[i])
                    + (Q1[-1] * w_Q1[i])
                )
            )
            p3 = (
                Q3[-1]
                * w_Q3[i]
                / (
                    (Q4[-1] * w_Q4[i])
                    + (Q3[-1] * w_Q3[i])
                    + (Q2[-1] * w_Q2[i])
                    + (Q1[-1] * w_Q1[i])
                )
            )
            p2 = (
                Q2[-1]
                * w_Q2[i]
                / (
                    (Q4[-1] * w_Q4[i])
                    + (Q3[-1] * w_Q3[i])
                    + (Q2[-1] * w_Q2[i])
                    + (Q1[-1] * w_Q1[i])
                )
            )
            p1 = (
                Q1[-1]
                * w_Q1[i]
                / (
                    (Q4[-1] * w_Q4[i])
                    + (Q3[-1] * w_Q3[i])
                    + (Q2[-1] * w_Q2[i])
                    + (Q1[-1] * w_Q1[i])
                )
            )

            if Q4[-1] - p4 < 0:
                next_Q4 = 0
            else:
                next_Q4 = Q4[-1] - p4

            if Q3[-1] + p4 - p3 < 0:
                next_Q3 = 0
            else:
                next_Q3 = Q3[-1] + p4 - p3

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

            Q4.append(next_Q4)
            Q3.append(next_Q3)
            Q2.append(next_Q2)
            Q1.append(next_Q1)
            Q0.append(next_Q0)

    elif type(H1) is tuple:
        w_Na_Q4 = []
        w_Na_Q3 = []
        w_Na_Q2 = []
        w_Na_Q1 = []

        w_K_Q4 = []
        w_K_Q3 = []
        w_K_Q2 = []
        w_K_Q1 = []

        for i in draw_ar:
            next_w_Na_Q4 = math.exp(-H1[1][0] / (Tg[i] * 0.00831))
            w_Na_Q4.append(next_w_Na_Q4)

            next_w_Na_Q3 = math.exp(-H1[1][1] / (Tg[i] * 0.00831))
            w_Na_Q3.append(next_w_Na_Q3)

            next_w_Na_Q2 = math.exp(-H1[1][2] / (Tg[i] * 0.00831))
            w_Na_Q2.append(next_w_Na_Q2)

            next_w_Na_Q1 = math.exp(-H1[1][3] / (Tg[i] * 0.00831))
            w_Na_Q1.append(next_w_Na_Q1)

            next_w_K_Q4 = math.exp(-H1[0][0] / (Tg[i] * 0.00831))
            w_K_Q4.append(next_w_K_Q4)

            next_w_K_Q3 = math.exp(-H1[0][1] / (Tg[i] * 0.00831))
            w_K_Q3.append(next_w_K_Q3)

            next_w_K_Q2 = math.exp(-H1[0][2] / (Tg[i] * 0.00831))
            w_K_Q2.append(next_w_K_Q2)

            next_w_K_Q1 = math.exp(-H1[0][3] / (Tg[i] * 0.00831))
            w_K_Q1.append(next_w_K_Q1)

        Q4 = [
            100,
        ]
        Q3 = [
            0,
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

            p4 = (
                (
                    Q4[-1]
                    * w_K_Q4[i]
                    / (
                        (Q4[-1] * w_K_Q4[i])
                        + (Q3[-1] * w_K_Q3[i])
                        + (Q2[-1] * w_K_Q2[i])
                        + (Q1[-1] * w_K_Q1[i])
                    )
                )
                * frac[0]
            ) + (
                (
                    Q4[-1]
                    * w_Na_Q4[i]
                    / (
                        (Q4[-1] * w_Na_Q4[i])
                        + (Q3[-1] * w_Na_Q3[i])
                        + (Q2[-1] * w_Na_Q2[i])
                        + (Q1[-1] * w_Na_Q1[i])
                    )
                )
                * frac[1]
            )
            p3 = (
                (
                    Q3[-1]
                    * w_K_Q3[i]
                    / (
                        (Q4[-1] * w_K_Q4[i])
                        + (Q3[-1] * w_K_Q3[i])
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
                        (Q4[-1] * w_Na_Q4[i])
                        + (Q3[-1] * w_Na_Q3[i])
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
                        (Q4[-1] * w_K_Q4[i])
                        + (Q3[-1] * w_K_Q3[i])
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
                        (Q4[-1] * w_Na_Q4[i])
                        + (Q3[-1] * w_Na_Q3[i])
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
                        (Q4[-1] * w_K_Q4[i])
                        + (Q3[-1] * w_K_Q3[i])
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
                        (Q4[-1] * w_Na_Q4[i])
                        + (Q3[-1] * w_Na_Q3[i])
                        + (Q2[-1] * w_Na_Q2[i])
                        + (Q1[-1] * w_Na_Q1[i])
                    )
                )
                * frac[1]
            )

            if Q4[-1] - p4 < 0:
                next_Q4 = 0
            else:
                next_Q4 = Q4[-1] - p4

            if Q3[-1] + p4 - p3 < 0:
                next_Q3 = 0
            else:
                next_Q3 = Q3[-1] + p4 - p3

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

            Q4.append(next_Q4)
            Q3.append(next_Q3)
            Q2.append(next_Q2)
            Q1.append(next_Q1)
            Q0.append(next_Q0)

    else:
        return print("Wrong H format")
    if s_plt is False and p is True:
        plt.plot(
            M2O,
            Q4,
            "r-",
            M2O,
            Q3,
            "k-",
            M2O,
            Q2,
            "b-",
            M2O,
            Q1,
            "g-",
            M2O,
            Q0,
            "y-",
        )
        plt.axis([0, 67, 0, 100])
        plt.legend(["$Q^4$", "$Q^3$", "$Q^2$", "$Q^1$", "$Q^0$"])
        plt.xlabel("Modifier mol %")
        plt.ylabel("Qn species concentration")
        plt.title("Qn distribution")
        plt.show()

    if s_plt is True:
        if not os.path.exists("SiO2_Structure"):
            os.mkdir("SiO2_Structure")
        plt.plot(
            M2O,
            Q4,
            "r-",
            M2O,
            Q3,
            "k-",
            M2O,
            Q2,
            "b-",
            M2O,
            Q1,
            "g-",
            M2O,
            Q0,
            "y-",
        )
        plt.axis([0, 67, 0, 100])
        plt.legend(["$Q^4$", "$Q^3$", "$Q^2$", "$Q^1$", "$Q^0$"])
        plt.xlabel("Modifier mol %")
        plt.ylabel("Qn species concentration")
        plt.title("Qn distribution")
        plt.savefig(os.path.join("SiO2_Structure", "Qn_distribution.png"))
        plt.show()

    if s_dat is True:
        if not os.path.exists("SiO2_Structure"):
            os.mkdir("SiO2_Structure")
        m_data = np.column_stack([M2O, Q4, Q3, Q2, Q1, Q0])
        np.savetxt(os.path.join("SiO2_Structure", "Model_data.csv"), m_data)

    elif s_plt is False and s_dat is False and p is False:
        return M2O, Q4, Q3, Q2, Q1, Q0


def Si_SSE(H1, data, tg, frac=None, s_plt=False, s_dat=False, p=False):
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
    Q4_data = data[1]
    Q3_data = data[2]
    Q2_data = data[3]
    Q1_data = data[4]
    Q0_data = data[5]

    M2O, Q4, Q3, Q2, Q1, Q0 = Si_draw(H1, tg, frac)

    if s_plt is False and p is True:
        plt.plot(
            M2O,
            Q4,
            "r-",
            M2O,
            Q3,
            "k-",
            M2O,
            Q2,
            "b-",
            M2O,
            Q1,
            "g-",
            M2O,
            Q0,
            "y-",
        )
        plt.plot(
            mod_data,
            Q4_data,
            "rd",
            mod_data,
            Q3_data,
            "kd",
            mod_data,
            Q2_data,
            "bd",
            mod_data,
            Q1_data,
            "gd",
            mod_data,
            Q0_data,
            "yd",
        )
        plt.axis([0, 67, 0, 100])
        plt.legend(["$Si^4$", "$Si^3$", "$Si^2$", "$Si^1$", "$Si^0$"])
        plt.xlabel("Modifier mol %")
        plt.ylabel("$Si^n$ species concentration")
        plt.title("$Si^n$ distribution")
        plt.show()

    if s_plt is True:
        if not os.path.exists("SiO2_Structure"):
            os.mkdir("SiO2_Structure")
        plt.plot(
            M2O,
            Q4,
            "r-",
            M2O,
            Q3,
            "k-",
            M2O,
            Q2,
            "b-",
            M2O,
            Q1,
            "g-",
            M2O,
            Q0,
            "y-",
        )
        plt.plot(
            mod_data,
            Q4_data,
            "rd",
            mod_data,
            Q3_data,
            "kd",
            mod_data,
            Q2_data,
            "bd",
            mod_data,
            Q1_data,
            "gd",
            mod_data,
            Q0_data,
            "yd",
        )
        plt.axis([0, 67, 0, 100])
        plt.legend(["$Si^4$", "$Si^3$", "$Si^2$", "$Si^1$", "$Si^0$"])
        plt.xlabel("Modifier mol %")
        plt.ylabel("$Si^n$ species concentration")
        plt.title("$Si^n$ distribution")
        plt.savefig(os.path.join("SiO2_Structure", "Qn_distribution.png"))
        plt.show()

    if s_dat is True:
        if not os.path.exists("SiO2_Structure"):
            os.mkdir("SiO2_Structure")
        m_data = np.column_stack([M2O, Q4, Q3, Q2, Q1, Q0])
        np.savetxt(os.path.join("SiO2_Structure", "Model_data.csv"), m_data)

    if p is False:
        mod_m = []
        Q4_m = []
        Q3_m = []
        Q2_m = []
        Q1_m = []
        Q0_m = []

        for i in mod_data:
            next_mod_m = min(M2O, key=lambda x: abs(x - i))
            mod_m.append(next_mod_m)

        for i in mod_m:
            ind = M2O.index(i)
            next_Q4_m = Q4[ind]
            Q4_m.append(next_Q4_m)

            next_Q3_m = Q3[ind]
            Q3_m.append(next_Q3_m)

            next_Q2_m = Q2[ind]
            Q2_m.append(next_Q2_m)

            next_Q1_m = Q1[ind]
            Q1_m.append(next_Q1_m)

            next_Q0_m = Q0[ind]
            Q0_m.append(next_Q0_m)

        Q3_m = np.array(Q3_m)
        Q4_m = np.array(Q4_m)
        Q2_m = np.array(Q2_m)
        Q1_m = np.array(Q1_m)
        Q0_m = np.array(Q0_m)

        SSE = sum(
            ((Q4_data - Q4_m) ** 2)
            + ((Q3_data - Q3_m) ** 2)
            + ((Q2_data - Q2_m) ** 2)
            + ((Q1_data - Q1_m) ** 2)
            + ((Q0_data - Q0_m) ** 2)
        )

        return SSE


def Si_engine(fil, data, tg, it=10):
    dat = data
    w0 = [10, 20, 30]

    minimizer_kwargs = {
        "method": "COBYLA",
        "args": (
            dat,
            tg,
        ),
    }
    res = scipy.optimize.basinhopping(
        Si_SSE,
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


# minimizer_kwargs = {"method": "BFGS"}
# w0 = [10, 20, 30]
# mybounds = MyBounds()

# res = scipy.optimize.basinhopping(model, w0, niter=50, T=2.0, stepsize=1,
# minimizer_kwargs=minimizer_kwargs, take_step=None,
# accept_test=None, callback=None, interval=50,
# disp=True, niter_success=None, seed=None)

# model(res.x,p=True)
