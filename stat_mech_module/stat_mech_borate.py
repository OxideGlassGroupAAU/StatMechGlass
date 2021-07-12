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


def B_onedraw(w, start_conc, draw_size, back=False):

    B3_s = start_conc[0]
    B4_s = start_conc[1]
    B2_s = start_conc[2]
    B1_s = start_conc[3]
    B0_s = start_conc[4]

    if (
        (B4_s + B2_s + B1_s * 2 + B0_s * 3) /
        (B3_s + B4_s + B2_s + B1_s * +B0_s)
    ) < 0.428:
        B4_B2 = 1
    else:
        B4_B2 = 0

    if back:

        wB0 = 1 / w[3]
        wB1 = 1 / w[2]
        wB2 = 1 / w[1]
        wB4 = 1 / w[0]

        p_B0 = B0_s * wB0 / (B0_s * wB0 + B1_s * wB1 + B2_s * wB2 + B4_s * wB4)
        p_B1 = B1_s * wB1 / (B0_s * wB0 + B1_s * wB1 + B2_s * wB2 + B4_s * wB4)
        p_B2 = B2_s * wB2 / (B0_s * wB0 + B1_s * wB1 + B2_s * wB2 + B4_s * wB4)
        p_B4 = B4_s * wB4 / (B0_s * wB0 + B1_s * wB1 + B2_s * wB2 + B4_s * wB4)

        p_B0 = p_B0 * draw_size
        p_B1 = p_B1 * draw_size
        p_B2 = p_B2 * draw_size
        p_B4 = p_B4 * draw_size

        # Contribution to N4 from B
        if p_B0 == 0 and p_B1 == 0 and p_B4 == 0:
            CB0 = 0
            CB1 = 0
            CB4 = 0
        else:
            CB0 = p_B0 / (p_B0 + p_B4 + p_B1)
            CB1 = p_B1 / (p_B0 + p_B4 + p_B1)
            CB4 = p_B4 / (p_B0 + p_B4 + p_B1)

        # Evolution of borate Qn units
        if B0_s + p_B0 + (p_B2 * B4_B2 * CB0) < 0:
            next_B0 = 0

        else:
            next_B0 = B0_s + p_B0 + (p_B2 * B4_B2 * CB0)

        if (B1_s - p_B0 + p_B1 + (CB1 * B4_B2 * p_B2) -
                (CB0 * B4_B2 * p_B2) < 0):
            next_B1 = 0
        else:
            next_B1 = (B1_s - p_B0 + p_B1 + (CB1 * B4_B2 * p_B2) -
                       (CB0 * B4_B2 * p_B2))

        if B2_s - p_B1 + p_B2 - (CB1 * B4_B2 * p_B2) < 0:
            next_B2 = 0
        else:
            next_B2 = B2_s - p_B1 + p_B2 - (CB1 * B4_B2 * p_B2)

        if B4_s - (p_B2 * B4_B2) + p_B4 + (CB4 * B4_B2 * p_B2) < 0:
            next_B4 = 0
            p_B4 = -B4_s
        else:
            next_B4 = B4_s - (p_B2 * B4_B2) + p_B4 + (CB4 * B4_B2 * p_B2)

        if B3_s - p_B4 - (p_B2 * (1 - B4_B2)) - (CB4 * B4_B2 * p_B2) < 0:
            next_B3 = 0
        else:
            next_B3 = B3_s - p_B4 - (p_B2 * (1 - B4_B2)) - (CB4 * B4_B2 * p_B2)

    else:

        p_B3 = B3_s * w[0] / (B3_s * w[0] + B4_s * w[1] +
                              B2_s * w[2] + B1_s * w[3])
        p_B4 = B4_s * w[1] / (B3_s * w[0] + B4_s * w[1] +
                              B2_s * w[2] + B1_s * w[3])
        p_B2 = B2_s * w[2] / (B3_s * w[0] + B4_s * w[1] +
                              B2_s * w[2] + B1_s * w[3])
        p_B1 = B1_s * w[3] / (B3_s * w[0] + B4_s * w[1] +
                              B2_s * w[2] + B1_s * w[3])

        p_B3 = p_B3 * draw_size
        p_B4 = p_B4 * draw_size
        p_B2 = p_B2 * draw_size
        p_B1 = p_B1 * draw_size

        # Contribution to N4 from B
        CB3 = p_B3 / (p_B3 + p_B2 + p_B1)
        CB2 = p_B2 / (p_B3 + p_B2 + p_B1)
        CB1 = p_B1 / (p_B3 + p_B2 + p_B1)

        # Evolution of borate Qn units
        if B3_s - p_B3 - (p_B4 * CB3) < 0:
            next_B3 = 0
        else:
            next_B3 = B3_s - p_B3 - (p_B4 * CB3)

        if B4_s + (p_B3 * B4_B2) - p_B4 < 0:
            next_B4 = 0
        else:
            next_B4 = B4_s + (p_B3 * B4_B2) - p_B4

        if (B2_s + (p_B3 * (1 - B4_B2)) + (p_B4 * CB3) + p_B4 -
                p_B2 - (p_B4 * CB2) < 0):
            next_B2 = 0
        else:
            next_B2 = (
                B2_s + (p_B3 * (1 - B4_B2)) + (p_B4 * CB3) +
                p_B4 - p_B2 - (p_B4 * CB2)
            )

        if B1_s + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1) < 0:
            next_B1 = 0
        else:
            next_B1 = B1_s + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1)

        if B0_s + p_B1 + (p_B4 * CB1) < 0:
            next_B0 = 0
        else:
            next_B0 = B0_s + p_B1 + (p_B4 * CB1)

    return next_B3, next_B4, next_B2, next_B1, next_B0


def B_back_onedraw(w, start_conc, draw_size):

    B3_s = start_conc[0]
    B4_s = start_conc[1]
    B2_s = start_conc[2]
    B1_s = start_conc[3]
    B0_s = start_conc[4]

    sum_S = sum(start_conc)

    if (
        (B4_s + B2_s + B1_s * 2 + B0_s * 3) / (B3_s + B4_s +
                                               B2_s + B1_s * +B0_s)
    ) < 0.428:
        B4_B2 = 1
    else:
        B4_B2 = 0

    B3_res = sum_S - B3_s
    B4_res = sum_S - B4_s
    B2_res = sum_S - B2_s
    B1_res = sum_S - B1_s

    p_B3 = (
        B3_res * w[0] / (B3_res * w[0] + B4_res * w[1] + B2_res *
                         w[2] + B1_res * w[3])
    )
    p_B4 = (
        B4_res * w[1] / (B3_res * w[0] + B4_res * w[1] + B2_res *
                         w[2] + B1_res * w[3])
    )
    p_B2 = (
        B2_res * w[2] / (B3_res * w[0] + B4_res * w[1] + B2_res *
                         w[2] + B1_res * w[3])
    )
    p_B1 = (
        B1_res * w[3] / (B3_res * w[0] + B4_res * w[1] + B2_res *
                         w[2] + B1_res * w[3])
    )

    p_B3 = p_B3 * draw_size
    p_B4 = p_B4 * draw_size
    p_B2 = p_B2 * draw_size
    p_B1 = p_B1 * draw_size

    # Contribution to N4 from B
    CB3 = p_B3 / (p_B3 + p_B2 + p_B1)
    CB2 = p_B2 / (p_B3 + p_B2 + p_B1)
    CB1 = p_B1 / (p_B3 + p_B2 + p_B1)

    # Evolution of borate Qn units
    if B3_s - p_B3 - (p_B4 * CB3) < 0:
        next_B3 = 0
    else:
        next_B3 = B3_s - p_B3 - (p_B4 * CB3)

    if B4_s + (p_B3 * B4_B2) - p_B4 < 0:
        next_B4 = 0
    else:
        next_B4 = B4_s + (p_B3 * B4_B2) - p_B4

    if (B2_s + (p_B3 * (1 - B4_B2)) + (p_B4 * CB3) +
            p_B4 - p_B2 - (p_B4 * CB2) < 0):
        next_B2 = 0
    else:
        next_B2 = (
            B2_s + (p_B3 * (1 - B4_B2)) + (p_B4 * CB3) + p_B4
            - p_B2 - (p_B4 * CB2)
        )

    if B1_s + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1) < 0:
        next_B1 = 0
    else:
        next_B1 = B1_s + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1)

    if B0_s + p_B1 + (p_B4 * CB1) < 0:
        next_B0 = 0
    else:
        next_B0 = B0_s + p_B1 + (p_B4 * CB1)

    return next_B3, next_B4, next_B2, next_B1, next_B0


def B_draw(w1, frac=None, s_plt=False, s_dat=False, p=False):

    draw_nr = list(range(300))
    draw_ar = np.array(draw_nr)

    M2O = []

    for i in draw_ar:
        next_mod = draw_ar[i] / (100 + (draw_ar[i])) * 100
        M2O.append(next_mod)

    M2O.append(75)

    M2Onp = np.array(M2O)

    Tg = np.array(0.0014 * M2Onp ** 3 - 0.3315 * M2Onp **
                  2 + 16.459 * M2Onp + 508.94)

    H = np.array([0, abs(w1[1]), abs(w1[2]), abs(w1[3])])

    w_B3 = []
    w_B4 = []
    w_B2 = []
    w_B1 = []

    for i in draw_ar:
        next_w_B3 = math.exp(-H[0] / (Tg[i] * 0.008314462))
        w_B3.append(next_w_B3)

    for i in draw_ar:
        next_w_B4 = math.exp(-H[1] / (Tg[i] * 0.008314462))
        w_B4.append(next_w_B4)

    for i in draw_ar:
        next_w_B2 = math.exp(-H[2] / (Tg[i] * 0.008314462))
        w_B2.append(next_w_B2)

    for i in draw_ar:
        next_w_B1 = math.exp(-H[2] / (Tg[i] * 0.008314462))
        w_B1.append(next_w_B1)

    B4_B2 = []
    for i in draw_ar:
        if M2O[i] < w1[0]:
            B4_B2.append(1)
        else:
            B4_B2.append(0)

    B4_B2 = np.array(B4_B2)

    B3 = [
        100,
    ]
    B4 = [
        0,
    ]
    B2 = [
        0,
    ]
    B1 = [
        0,
    ]
    B0 = [
        0,
    ]

    # The function that makes each iteration at a time

    for i in draw_ar:

        p_B3 = (
            B3[-1]
            * w_B3[i]
            / (
                B3[-1] * w_B3[i]
                + B4[-1] * w_B4[i]
                + B2[-1] * w_B2[i]
                + B1[-1] * w_B1[i]
            )
        )
        p_B4 = (
            B4[-1]
            * w_B4[i]
            / (
                B3[-1] * w_B3[i]
                + B4[-1] * w_B4[i]
                + B2[-1] * w_B2[i]
                + B1[-1] * w_B1[i]
            )
        )
        p_B2 = (
            B2[-1]
            * w_B2[i]
            / (
                B3[-1] * w_B3[i]
                + B4[-1] * w_B4[i]
                + B2[-1] * w_B2[i]
                + B1[-1] * w_B1[i]
            )
        )
        p_B1 = (
            B1[-1]
            * w_B1[i]
            / (
                B3[-1] * w_B3[i]
                + B4[-1] * w_B4[i]
                + B2[-1] * w_B2[i]
                + B1[-1] * w_B1[i]
            )
        )

        # Contribution to N4 from B
        CB3 = p_B3 / (p_B3 + p_B2 + p_B1)
        CB2 = p_B2 / (p_B3 + p_B2 + p_B1)
        CB1 = p_B1 / (p_B3 + p_B2 + p_B1)

        # Evolution of borate Qn units
        if B3[-1] - p_B3 - (p_B4 * CB3) < 0:
            next_B3 = 0
        else:
            next_B3 = B3[-1] - p_B3 - (p_B4 * CB3)

        if B4[-1] + (p_B3 * B4_B2[i]) - p_B4 < 0:
            next_B4 = 0
        else:
            next_B4 = B4[-1] + (p_B3 * B4_B2[i]) - p_B4

        if (
            B2[-1] + (p_B3 * (1 - B4_B2[i])) + (p_B4 * CB3) +
            p_B4 - p_B2 - (p_B4 * CB2)
            < 0
        ):
            next_B2 = 0
        else:
            next_B2 = (
                B2[-1]
                + (p_B3 * (1 - B4_B2[i]))
                + (p_B4 * CB3)
                + p_B4
                - p_B2
                - (p_B4 * CB2)
            )

        if B1[-1] + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1) < 0:
            next_B1 = 0
        else:
            next_B1 = B1[-1] + p_B2 - p_B1 + (p_B4 * CB2) - (p_B4 * CB1)

        if B0[-1] + p_B1 + (p_B4 * CB1) < 0:
            next_B0 = 0
        else:
            next_B0 = B0[-1] + p_B1 + (p_B4 * CB1)

        B3.append(next_B3)
        B4.append(next_B4)
        B2.append(next_B2)
        B1.append(next_B1)
        B0.append(next_B0)

    return M2O, B3, B4, B2, B1, B0


def B_SSE(w1, data, frac=None, s_plt=False, s_dat=False, p=False):

    mod_data = data[0]
    Q4_data = data[1]
    M2O, B3, B4, B2, B1, B0 = B_draw(w1)

    if s_plt is False and p is True:
        plt.plot(
            M2O,
            B3,
            "r-",
            M2O,
            B4,
            "k-",
            M2O,
            B2,
            "b-",
            M2O,
            B1,
            "g-",
            M2O,
            B0,
            "y-",
        )
        plt.plot(mod_data, Q4_data, "kd")
        plt.axis([0, 75, 0, 100])
        plt.legend(["$B^3$", "$B^4$", "$B^2$", "$B^1$", "$B^0$"])
        plt.xlabel("Modifier mol %")
        plt.ylabel("Bn species concentration")
        plt.title("Bn distribution")
        plt.show()

    if s_plt is True:
        if not os.path.exists("B2O3_Structure"):
            os.mkdir("B2O3_Structure")
        plt.plot(
            M2O,
            B3,
            "r-",
            M2O,
            B4,
            "k-",
            M2O,
            B2,
            "b-",
            M2O,
            B1,
            "g-",
            M2O,
            B0,
            "y-",
        )
        plt.plot(mod_data, Q4_data, "kd")
        plt.axis([0, 75, 0, 100])
        plt.legend(["$B^3$", "$B^4$", "$B^2$", "$B^1$", "$B^0$"])
        plt.xlabel("Modifier mol %")
        plt.ylabel("Bn species concentration")
        plt.title("Bn distribution")
        plt.savefig(os.path.join("B2O3_Structure", "Qn_distribution.png"))
        plt.show()

    if s_dat is True:
        if not os.path.exists("B2O3_Structure"):
            os.mkdir("B2O3_Structure")
        m_data = np.column_stack([M2O, B3, B4, B2, B1, B0])
        np.savetxt(os.path.join("SiO2_Structure", "Model_data.csv"), m_data)

    if p is False:

        mod_m = []
        B3_m = []
        B4_m = []
        B2_m = []
        B1_m = []
        B0_m = []

        for i in mod_data:
            next_mod_m = min(M2O, key=lambda x: abs(x - i))
            mod_m.append(next_mod_m)

        for i in mod_m:
            ind = M2O.index(i)
            next_B3_m = B3[ind]
            B3_m.append(next_B3_m)

            next_B4_m = B4[ind]
            B4_m.append(next_B4_m)

            next_B2_m = B2[ind]
            B2_m.append(next_B2_m)

            next_B1_m = B1[ind]
            B1_m.append(next_B1_m)

            next_B0_m = B0[ind]
            B0_m.append(next_B0_m)

        B3_m = np.array(B3_m)
        B4_m = np.array(B4_m)
        B2_m = np.array(B2_m)
        B1_m = np.array(B1_m)
        B0_m = np.array(B0_m)

        SSE = sum(((Q4_data - B4_m) ** 2))

        return SSE


def B_engine(fil, data, it=10):
    dat = data
    w0 = [35, 10, 20, 30]

    minimizer_kwargs = {"method": "COBYLA", "args": (dat,)}
    res = scipy.optimize.basinhopping(
        B_SSE,
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
