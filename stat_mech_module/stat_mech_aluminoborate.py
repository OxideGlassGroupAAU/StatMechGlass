import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import os


def AlB_first_draw(w1, start_conc, former):

    if former == "Si":
        w = np.array([1, abs(w1[0]), abs(w1[1]), abs(w1[2]), abs(w1[3])])

        draw_nr = list(range(401))
        draw_ar = np.array(draw_nr)

        M2O = []

        for i in draw_ar:
            next_mod = draw_ar[i] * 0.5 / (100 + draw_ar[i] * 0.5) * 100
            M2O.append(next_mod)

        M2O.append(67)

        # Startværdier

        # for ind, i in enumerate(r_in, start=0):
        r = start_conc[5] / start_conc[0]

        Q4 = [
            start_conc[0],
        ]
        Al5 = [
            start_conc[5],
        ]

        Al4A = [
            start_conc[6],
        ]

        Q3A = [
            start_conc[1],
        ]
        Q2A = [
            start_conc[2],
        ]
        Q1AA = [
            start_conc[3],
        ]
        Q0AAA = [
            start_conc[4],
        ]

        Al_draw = [
            0,
        ]

        # Draw of all Al
        if Al5[-1] > 0:
            while (Al_draw[-1] / 3) * 100 / (100 + (Al_draw[-1] / 3)) < Al5[
                -1
            ] * 100 / (100 - Al4A[-1]):
                gQ4 = (
                    (Q4[-1] * w[0])
                    * (1 / (r + 1))
                    / (
                        Q4[-1] * w[0]
                        + Q3A[-1] * w[1]
                        + Q2A[-1] * w[2]
                        + Q1AA[-1] * w[3]
                        + Al5[-1] * w[4]
                    )
                )
                gQ3 = (
                    (Q3A[-1] * w[1])
                    * (1 / (r + 1))
                    / (
                        Q4[-1] * w[0]
                        + Q3A[-1] * w[1]
                        + Q2A[-1] * w[2]
                        + Q1AA[-1] * w[3]
                        + Al5[-1] * w[4]
                    )
                )
                gQ2 = (
                    (Q2A[-1] * w[2])
                    * (1 / (r + 1))
                    / (
                        Q4[-1] * w[0]
                        + Q3A[-1] * w[1]
                        + Q2A[-1] * w[2]
                        + Q1AA[-1] * w[3]
                        + Al5[-1] * w[4]
                    )
                )
                gQ1 = (
                    (Q1AA[-1] * w[3])
                    * (1 / (r + 1))
                    / (
                        Q4[-1] * w[0]
                        + Q3A[-1] * w[1]
                        + Q2A[-1] * w[2]
                        + Q1AA[-1] * w[3]
                        + Al5[-1] * w[4]
                    )
                )
                gAl5 = (
                    (Al5[-1] * w[4])
                    * (1 / (r + 1))
                    / (
                        Q4[-1] * w[0]
                        + Q3A[-1] * w[1]
                        + Q2A[-1] * w[2]
                        + Q1AA[-1] * w[3]
                        + Al5[-1] * w[4]
                    )
                )

                # overgang

                # Draws
                if Q4[-1] - gQ4 > 0:
                    next_Q4 = Q4[-1] - gQ4
                else:
                    next_Q4 = 0

                if Q3A[-1] + gQ4 - gQ3 > 0:
                    next_Q3A = Q3A[-1] + gQ4 - gQ3
                else:
                    next_Q3A = 0

                if Q2A[-1] + gQ3 - gQ2 > 0:
                    next_Q2A = Q2A[-1] + gQ3 - gQ2

                else:
                    next_Q2A = 0

                if Q1AA[-1] + gQ2 - gQ1 < 0:
                    next_Q1AA = 0
                else:
                    next_Q1AA = Q1AA[-1] + gQ2 - gQ1

                if Q0AAA[-1] + gQ1 < 0:
                    next_Q0AAA = 0
                else:
                    next_Q0AAA = Q0AAA[-1] + gQ1

                if Al5[-1] - gAl5 > 0:
                    next_Al5 = Al5[-1] - gAl5
                else:
                    next_Al5 = 0

                if Al4A[-1] + gAl5 > 0:
                    next_Al4A = Al4A[-1] + gAl5
                else:
                    next_Al4A = 0

                next_Al_draw = Al_draw[-1] + 1
                Al_draw.append(next_Al_draw)

                Q4.append(next_Q4)
                Q3A.append(next_Q3A)
                Q2A.append(next_Q2A)
                Q1AA.append(next_Q1AA)
                Q0AAA.append(next_Q0AAA)

                Al5.append(next_Al5)
                Al4A.append(next_Al4A)

            # Residual Al draw ##
            aa = (
                (-300 * (Al5[-1] * 100 / (100 - Al4A[-1])))
                / (-100 + (Al5[-1] * 100 / (100 - Al4A[-1])))
            ) - (-300 * (Al5[-2] * 100 / (100 - Al4A[-2]))) / (
                -100 + (Al5[-2] * 100 / (100 - Al4A[-2]))
            )
            b = -aa * Al_draw[-1] + (-300 * (Al5[-1] * 100 /
                                             (100 - Al4A[-1]))) / (
                -100 + (Al5[-1] * 100 / (100 - Al4A[-1]))
            )
            balance_draw = (-b / (aa - 1)) - Al_draw[-2]

            gQ4 = (
                balance_draw
                * (Q4[-2] * w[0])
                * (1 / (r + 1))
                / (
                    Q4[-2] * w[0]
                    + Q3A[-2] * w[1]
                    + Q2A[-2] * w[2]
                    + Q1AA[-2] * w[3]
                    + Al5[-2] * w[4]
                )
            )
            gQ3 = (
                balance_draw
                * (Q3A[-2] * w[1])
                * (1 / (r + 1))
                / (
                    Q4[-2] * w[0]
                    + Q3A[-2] * w[1]
                    + Q2A[-2] * w[2]
                    + Q1AA[-2] * w[3]
                    + Al5[-2] * w[4]
                )
            )
            gQ2 = (
                balance_draw
                * (Q2A[-2] * w[2])
                * (1 / (r + 1))
                / (
                    Q4[-2] * w[0]
                    + Q3A[-2] * w[1]
                    + Q2A[-2] * w[2]
                    + Q1AA[-2] * w[3]
                    + Al5[-2] * w[4]
                )
            )
            gQ1 = (
                balance_draw
                * (Q1AA[-2] * w[3])
                * (1 / (r + 1))
                / (
                    Q4[-2] * w[0]
                    + Q3A[-2] * w[1]
                    + Q2A[-2] * w[2]
                    + Q1AA[-2] * w[3]
                    + Al5[-2] * w[4]
                )
            )
            gAl5 = (
                balance_draw
                * (Al5[-2] * w[4])
                * (1 / (r + 1))
                / (
                    Q4[-2] * w[0]
                    + Q3A[-2] * w[1]
                    + Q2A[-2] * w[2]
                    + Q1AA[-2] * w[3]
                    + Al5[-2] * w[4]
                )
            )

            return (
                    Q4[-1], Q3A[-1], Q2A[-1], Q1AA[-1],
                    Q0AAA[-1], Al5[-1], Al4A[-1]
                    )

    if former == "B":
        w = np.array([1, abs(w1[0]), abs(w1[1]), abs(w1[2]), abs(w1[3]), 30])

        draw_nr = list(range(301))
        draw_ar = np.array(draw_nr)

        M2O = []

        for i in draw_ar:
            next_mod = draw_ar[i] / (100 + draw_ar[i]) * 100
            M2O.append(next_mod)

        M2O.append(75)

        # Startværdier

        # for ind, i in enumerate(r_in, start=0):
        r = start_conc[5] / start_conc[0]

        Q3 = [
            start_conc[0],
        ]
        Al5 = [
            start_conc[5],
        ]

        Al4A = [
            start_conc[6],
        ]

        Q4A = [
            start_conc[1],
        ]
        Q2A = [
            start_conc[2],
        ]
        Q1AA = [
            start_conc[3],
        ]
        Q0AAA = [
            start_conc[4],
        ]

        Al_draw = [
            0,
        ]

        # Draw of all Al
        if Al5[-1] > 0:
            while (Al_draw[-1] / 3) * 100 / (100 + (Al_draw[-1] / 3)) < Al5[
                -1
            ] * 100 / (100 - Al4A[-1]):
                gQ3 = (
                    (Q3[-1] * w[0])
                    * (1 / (r + 1))
                    / (
                        Q3[-1] * w[0]
                        + Q4A[-1] * w[1]
                        + Q2A[-1] * w[2]
                        + Q1AA[-1] * w[3]
                        + Al5[-1] * w[4]
                    )
                )
                gQ4 = (
                    (Q4A[-1] * w[1])
                    * (1 / (r + 1))
                    / (
                        Q3[-1] * w[0]
                        + Q4A[-1] * w[1]
                        + Q2A[-1] * w[2]
                        + Q1AA[-1] * w[3]
                        + Al5[-1] * w[4]
                    )
                )
                gQ2 = (
                    (Q2A[-1] * w[2])
                    * (1 / (r + 1))
                    / (
                        Q3[-1] * w[0]
                        + Q4A[-1] * w[1]
                        + Q2A[-1] * w[2]
                        + Q1AA[-1] * w[3]
                        + Al5[-1] * w[4]
                    )
                )
                gQ1 = (
                    (Q1AA[-1] * w[3])
                    * (1 / (r + 1))
                    / (
                        Q3[-1] * w[0]
                        + Q4A[-1] * w[1]
                        + Q2A[-1] * w[2]
                        + Q1AA[-1] * w[3]
                        + Al5[-1] * w[4]
                    )
                )
                gAl5 = (
                    (Al5[-1] * w[4])
                    * (1 / (r + 1))
                    / (
                        Q3[-1] * w[0]
                        + Q4A[-1] * w[1]
                        + Q2A[-1] * w[2]
                        + Q1AA[-1] * w[3]
                        + Al5[-1] * w[4]
                    )
                )

                rgQ3 = gQ3 / (gQ3 + gQ2 + gQ1 + gAl5)
                rgQ2 = gQ2 / (gQ3 + gQ2 + gQ1 + gAl5)
                rgQ1 = gQ1 / (gQ3 + gQ2 + gQ1 + gAl5)
                rgAl5 = gAl5 / (gQ3 + gQ2 + gQ1 + gAl5)

                if (
                        (Q4A[-1] + Q2A[-1] + 2 * Q1AA[-1] + 3 * Q0AAA[-1])
                        < Q3[0] * w[5]):
                    P = 1
                else:
                    P = 0

                # Draws
                if Q3[-1] - gQ3 + (-rgQ3) * gQ4 > 0:
                    next_Q3 = Q3[-1] - gQ3 + (-rgQ3) * gQ4

                else:
                    next_Q3 = 0

                if Q4A[-1] + gQ3 * P - gQ4 + (rgQ3 * P) * gQ4 > 0:
                    next_Q4A = Q4A[-1] + gQ3 * P - gQ4 + (rgQ3 * P) * gQ4

                else:
                    next_Q4A = 0

                if (
                    Q2A[-1] + gQ4 + gQ3 * (1 - P) - gQ2 +
                    (rgQ3 * (1 - P) - rgQ2) * gQ4 > 0
                ):
                    next_Q2A = (
                        Q2A[-1]
                        + gQ4
                        + gQ3 * (1 - P)
                        - gQ2
                        + (rgQ3 * (1 - P) - rgQ2) * gQ4
                    )

                else:
                    next_Q2A = 0

                if Q1AA[-1] + gQ2 - gQ1 + (rgQ2 - rgQ1) * gQ4 < 0:
                    next_Q1AA = 0

                else:
                    next_Q1AA = Q1AA[-1] + gQ2 - gQ1 + (rgQ2 - rgQ1) * gQ4

                if Q0AAA[-1] + gQ1 + (rgQ1) * gQ4 < 0:
                    next_Q0AAA = 0

                else:
                    next_Q0AAA = Q0AAA[-1] + gQ1 + (rgQ1) * gQ4

                if Al5[-1] - gAl5 + (-rgAl5) * gQ4 > 0:
                    next_Al5 = Al5[-1] - gAl5 + (-rgAl5) * gQ4
                else:
                    next_Al5 = 0

                if Al4A[-1] + gAl5 + (rgAl5) * gQ4 > 0:
                    next_Al4A = Al4A[-1] + gAl5 + (rgAl5) * gQ4
                else:
                    next_Al4A = 0

                next_Al_draw = Al_draw[-1] + 1
                Al_draw.append(next_Al_draw)

                Q3.append(next_Q3)
                Q4A.append(next_Q4A)
                Q2A.append(next_Q2A)
                Q1AA.append(next_Q1AA)
                Q0AAA.append(next_Q0AAA)

                Al5.append(next_Al5)
                Al4A.append(next_Al4A)

            # Residual Al draw
            aa = (
                (-300 * (Al5[-1] * 100 / (100 - Al4A[-1])))
                / (-100 + (Al5[-1] * 100 / (100 - Al4A[-1])))
            ) - (-300 * (Al5[-2] * 100 / (100 - Al4A[-2]))) / (
                -100 + (Al5[-2] * 100 / (100 - Al4A[-2]))
            )
            b = -aa * Al_draw[-1] + (-300 * (Al5[-1] * 100 /
                                             (100 - Al4A[-1]))) / (
                -100 + (Al5[-1] * 100 / (100 - Al4A[-1]))
            )
            balance_draw = (-b / (aa - 1)) - Al_draw[-2]

            gQ3 = (
                balance_draw
                * (Q3[-2] * w[0])
                * (1 / (r + 1))
                / (
                    Q3[-2] * w[0]
                    + Q4A[-2] * w[1]
                    + Q2A[-2] * w[2]
                    + Q1AA[-2] * w[3]
                    + Al5[-2] * w[4]
                )
            )
            gQ4 = (
                balance_draw
                * (Q4A[-2] * w[1])
                * (1 / (r + 1))
                / (
                    Q3[-2] * w[0]
                    + Q4A[-2] * w[1]
                    + Q2A[-2] * w[2]
                    + Q1AA[-2] * w[3]
                    + Al5[-2] * w[4]
                )
            )
            gQ2 = (
                balance_draw
                * (Q2A[-2] * w[2])
                * (1 / (r + 1))
                / (
                    Q3[-2] * w[0]
                    + Q4A[-2] * w[1]
                    + Q2A[-2] * w[2]
                    + Q1AA[-2] * w[3]
                    + Al5[-2] * w[4]
                )
            )
            gQ1 = (
                balance_draw
                * (Q1AA[-2] * w[3])
                * (1 / (r + 1))
                / (
                    Q3[-2] * w[0]
                    + Q4A[-2] * w[1]
                    + Q2A[-2] * w[2]
                    + Q1AA[-2] * w[3]
                    + Al5[-2] * w[4]
                )
            )
            gAl5 = (
                balance_draw
                * (Al5[-2] * w[4])
                * (1 / (r + 1))
                / (
                    Q3[-2] * w[0]
                    + Q4A[-2] * w[1]
                    + Q2A[-2] * w[2]
                    + Q1AA[-2] * w[3]
                    + Al5[-2] * w[4]
                )
            )

            rgQ3 = gQ3 / (gQ3 + gQ2 + gQ1 + gAl5)
            rgQ2 = gQ2 / (gQ3 + gQ2 + gQ1 + gAl5)
            rgQ1 = gQ1 / (gQ3 + gQ2 + gQ1 + gAl5)
            rgAl5 = gAl5 / (gQ3 + gQ2 + gQ1 + gAl5)

            return (Q3[-1], Q4A[-1], Q2A[-1], Q1AA[-1],
                    Q0AAA[-1], Al5[-1], Al4A[-1])


def AlB_one_draw(w, start_conc, draw_size):

    Al5_s = start_conc[0]
    Al4_s = start_conc[1]

    p5 = draw_size

    if Al5_s - p5 < 0:
        next_Al5 = 0
    else:
        next_Al5 = Al5_s - p5

    next_Al4 = Al4_s + p5

    return next_Al5, next_Al4


def AlB_draw(w1, r, mod):
    w = np.array(
        [
            1,
            abs(w1[0]),
            abs(w1[1]),
            abs(w1[2]),
            abs(w1[3]),
            abs(w1[4]),
            abs(w1[5]),
            abs(w1[6]),
            abs(w1[7]),
            abs(w1[8]),
            abs(w1[9]),
            abs(w1[10]),
        ]
    )
    iw = 1 / w

    draw_nr = list(range(301))
    draw_ar = np.array(draw_nr)

    M2O = []

    for i in draw_ar:
        next_mod = draw_ar[i] / (100 + draw_ar[i]) * 100
        M2O.append(next_mod)

    M2O.append(75)

    # Startværdier

    # for ind, i in enumerate(r_in, start=0):

    Q3 = [
        (100 / (r + 1)),
    ]
    Al5 = [
        (100 * r / (r + 1)),
    ]

    Al4A = [
        0,
    ]
    Al4N = [
        0,
    ]

    # Q4
    Q4A = [
        0,
    ]
    Q4N = [
        0,
    ]

    # Q2
    Q2A = [
        0,
    ]
    Q2N = [
        0,
    ]

    # Q1
    Q1NN = [
        0,
    ]
    Q1AN = [
        0,
    ]
    Q1AA = [
        0,
    ]

    # Q0
    Q0AAA = [
        0,
    ]
    Q0AAN = [
        0,
    ]
    Q0ANN = [
        0,
    ]
    Q0NNN = [
        0,
    ]

    Al_draw = [
        0,
    ]

    Q3Na = []
    Al5Na = []
    Al4Na = []
    Q4Na = []
    Q2Na = []
    Q1Na = []
    Q0Na = []

    Q3shinanigens = []
    # Draw of all Al
    if Al5[-1] > 0:
        while ((Al_draw[-1] / 3) * 100 / (100 + (Al_draw[-1] / 3)) <
                Al5[-1] * 100 / (100 - Al4A[-1])):
            gQ3 = (
                (Q3[-1] * w[6])
                * (1 / (r + 1))
                / (
                    Q3[-1] * w[6]
                    + Q4A[-1] * w[7]
                    + Q2A[-1] * w[8]
                    + Q1AA[-1] * w[9]
                    + Al5[-1] * w[10]
                )
            )
            gQ4 = (
                (Q4A[-1] * w[7])
                * (1 / (r + 1))
                / (
                    Q3[-1] * w[6]
                    + Q4A[-1] * w[7]
                    + Q2A[-1] * w[8]
                    + Q1AA[-1] * w[9]
                    + Al5[-1] * w[10]
                )
            )
            gQ2 = (
                (Q2A[-1] * w[8])
                * (1 / (r + 1))
                / (
                    Q3[-1] * w[6]
                    + Q4A[-1] * w[7]
                    + Q2A[-1] * w[8]
                    + Q1AA[-1] * w[9]
                    + Al5[-1] * w[10]
                )
            )
            gQ1 = (
                (Q1AA[-1] * w[9])
                * (1 / (r + 1))
                / (
                    Q3[-1] * w[6]
                    + Q4A[-1] * w[7]
                    + Q2A[-1] * w[8]
                    + Q1AA[-1] * w[9]
                    + Al5[-1] * w[10]
                )
            )
            gAl5 = (
                (Al5[-1] * w[10])
                * (1 / (r + 1))
                / (
                    Q3[-1] * w[6]
                    + Q4A[-1] * w[7]
                    + Q2A[-1] * w[8]
                    + Q1AA[-1] * w[9]
                    + Al5[-1] * w[10]
                )
            )

            rgQ3 = gQ3 / (gQ3 + gQ2 + gQ1 + gAl5)
            rgQ2 = gQ2 / (gQ3 + gQ2 + gQ1 + gAl5)
            rgQ1 = gQ1 / (gQ3 + gQ2 + gQ1 + gAl5)
            rgAl5 = gAl5 / (gQ3 + gQ2 + gQ1 + gAl5)

            # overgang
            if ((Q4A[-1] + Q2A[-1] + 2 * Q1AA[-1] + 3 *
                 Q0AAA[-1]) < Q3[0] * w[5]):
                P = 1
            else:
                P = 0

            # Draws
            if Q3[-1] - gQ3 + (-rgQ3) * gQ4 > 0:
                next_Q3 = Q3[-1] - gQ3 + (-rgQ3) * gQ4

            else:
                next_Q3 = 0

            if Q4A[-1] + gQ3 * P - gQ4 + (rgQ3 * P) * gQ4 > 0:
                next_Q4A = Q4A[-1] + gQ3 * P - gQ4 + (rgQ3 * P) * gQ4

            else:
                next_Q4A = 0

            if (Q2A[-1] + gQ4 + gQ3 * (1 - P) - gQ2 +
                    (rgQ3 * (1 - P) - rgQ2) * gQ4 > 0):
                next_Q2A = (
                    Q2A[-1] + gQ4 + gQ3 * (1 - P) - gQ2 +
                    (rgQ3 * (1 - P) - rgQ2) * gQ4
                )

            else:
                next_Q2A = 0

            if Q1AA[-1] + gQ2 - gQ1 + (rgQ2 - rgQ1) * gQ4 < 0:
                next_Q1AA = 0

            else:
                next_Q1AA = Q1AA[-1] + gQ2 - gQ1 + (rgQ2 - rgQ1) * gQ4

            if Q0AAA[-1] + gQ1 + (rgQ1) * gQ4 < 0:
                next_Q0AAA = 0

            else:
                next_Q0AAA = Q0AAA[-1] + gQ1 + (rgQ1) * gQ4

            if Al5[-1] - gAl5 + (-rgAl5) * gQ4 > 0:
                next_Al5 = Al5[-1] - gAl5 + (-rgAl5) * gQ4
            else:
                next_Al5 = 0

            if Al4A[-1] + gAl5 + (rgAl5) * gQ4 > 0:
                next_Al4A = Al4A[-1] + gAl5 + (rgAl5) * gQ4
            else:
                next_Al4A = 0

            next_Al_draw = Al_draw[-1] + 1
            Al_draw.append(next_Al_draw)

            Q3.append(next_Q3)
            Q4A.append(next_Q4A)
            Q2A.append(next_Q2A)
            Q1AA.append(next_Q1AA)
            Q0AAA.append(next_Q0AAA)

            Al5.append(next_Al5)
            Al4A.append(next_Al4A)

        # Residual Al draw
        aa = (
            (-300 * (Al5[-1] * 100 / (100 - Al4A[-1])))
            / (-100 + (Al5[-1] * 100 / (100 - Al4A[-1])))
        ) - (-300 * (Al5[-2] * 100 / (100 - Al4A[-2]))) / (
            -100 + (Al5[-2] * 100 / (100 - Al4A[-2]))
        )
        b = -aa * Al_draw[-1] + (-300 * (Al5[-1] * 100 / (100 - Al4A[-1]))) / (
            -100 + (Al5[-1] * 100 / (100 - Al4A[-1]))
        )
        balance_draw = (-b / (aa - 1)) - Al_draw[-2]

        gQ3 = (
            balance_draw
            * (Q3[-2] * w[6])
            * (1 / (r + 1))
            / (
                Q3[-2] * w[6]
                + Q4A[-2] * w[7]
                + Q2A[-2] * w[8]
                + Q1AA[-2] * w[9]
                + Al5[-2] * w[10]
            )
        )
        gQ4 = (
            balance_draw
            * (Q4A[-2] * w[7])
            * (1 / (r + 1))
            / (
                Q3[-2] * w[6]
                + Q4A[-2] * w[7]
                + Q2A[-2] * w[8]
                + Q1AA[-2] * w[9]
                + Al5[-2] * w[10]
            )
        )
        gQ2 = (
            balance_draw
            * (Q2A[-2] * w[8])
            * (1 / (r + 1))
            / (
                Q3[-2] * w[6]
                + Q4A[-2] * w[7]
                + Q2A[-2] * w[8]
                + Q1AA[-2] * w[9]
                + Al5[-2] * w[10]
            )
        )
        gQ1 = (
            balance_draw
            * (Q1AA[-2] * w[9])
            * (1 / (r + 1))
            / (
                Q3[-2] * w[6]
                + Q4A[-2] * w[7]
                + Q2A[-2] * w[8]
                + Q1AA[-2] * w[9]
                + Al5[-2] * w[10]
            )
        )
        gAl5 = (
            balance_draw
            * (Al5[-2] * w[10])
            * (1 / (r + 1))
            / (
                Q3[-2] * w[6]
                + Q4A[-2] * w[7]
                + Q2A[-2] * w[8]
                + Q1AA[-2] * w[9]
                + Al5[-2] * w[10]
            )
        )

        rgQ3 = gQ3 / (gQ3 + gQ2 + gQ1 + gAl5)
        rgQ2 = gQ2 / (gQ3 + gQ2 + gQ1 + gAl5)
        rgQ1 = gQ1 / (gQ3 + gQ2 + gQ1 + gAl5)
        rgAl5 = gAl5 / (gQ3 + gQ2 + gQ1 + gAl5)

        # overgang
        if (Q4A[-2] + Q2A[-2] + 2 * Q1AA[-2] + 3 * Q0AAA[-2]) < Q3[0] * w[5]:
            P = 1
        else:
            P = 0

        # Draws
        if Q3[-2] - gQ3 + (-rgQ3) * gQ4 > 0:
            next_Q3 = Q3[-2] - gQ3 + (-rgQ3) * gQ4

        else:
            next_Q3 = 0

        if Q4A[-2] + gQ3 * P - gQ4 + (rgQ3 * P) * gQ4 > 0:
            next_Q4A = Q4A[-2] + gQ3 * P - gQ4 + (rgQ3 * P) * gQ4

        else:
            next_Q4A = 0

        if (Q2A[-2] + gQ4 + gQ3 * (1 - P) - gQ2 +
                (rgQ3 * (1 - P) - rgQ2) * gQ4 > 0):
            next_Q2A = (
                Q2A[-2] + gQ4 + gQ3 * (1 - P) - gQ2 +
                (rgQ3 * (1 - P) - rgQ2) * gQ4
            )

        else:
            next_Q2A = 0

        if Q1AA[-2] + gQ2 - gQ1 + (rgQ2 - rgQ1) * gQ4 < 0:
            next_Q1AA = 0

        else:
            next_Q1AA = Q1AA[-2] + gQ2 - gQ1 + (rgQ2 - rgQ1) * gQ4

        if Q0AAA[-2] + gQ1 + (rgQ1) * gQ4 < 0:
            next_Q0AAA = 0

        else:
            next_Q0AAA = Q0AAA[-2] + gQ1 + (rgQ1) * gQ4

        if Al5[-2] - gAl5 + (-rgAl5) * gQ4 > 0:
            next_Al5 = Al5[-2] - gAl5 + (-rgAl5) * gQ4
        else:
            next_Al5 = 0

        if Al4A[-2] + gAl5 + (rgAl5) * gQ4 > 0:
            next_Al4A = Al4A[-2] + gAl5 + (rgAl5) * gQ4
        else:
            next_Al4A = 0

        Q3.append(next_Q3)
        Q4A.append(next_Q4A)
        Q2A.append(next_Q2A)
        Q1AA.append(next_Q1AA)
        Q0AAA.append(next_Q0AAA)
        Al5.append(next_Al5)
        Al4A.append(next_Al4A)

        Q3Na.append(next_Q3)
        Al5Na.append(next_Al5)
        Al4Na.append(next_Al4A)
        Q4Na.append(next_Q4A)
        Q2Na.append(next_Q2A)
        Q1Na.append(next_Q1AA)
        Q0Na.append(next_Q0AAA)

    # For loop for udregning af g og Q_Al5

    for i in draw_ar:
        if Q4A[-1] + Q4N[-1] > 0:
            Q4 = Q4A[-1] + Q4N[-1]
        else:
            Q4 = 0

        if Q2A[-1] + Q2N[-1] > 0:
            Q2 = Q2A[-1] + Q2N[-1]
        else:
            Q2 = 0

        if Q1AA[-1] + Q1AN[-1] + Q1NN[-1] > 0:
            Q1 = Q1AA[-1] + Q1AN[-1] + Q1NN[-1]
        else:
            Q1 = 0

        if Q0AAA[-1] + Q0AAN[-1] + Q0ANN[-1] + Q0NNN[-1] > 0:
            Q0 = Q0AAA[-1] + Q0AAN[-1] + Q0ANN[-1] + Q0NNN[-1]
        else:
            Q0 = 0

        if Al4A[-1] + Al4N[-1] > 0:
            Al4 = Al4A[-1] + Al4N[-1]
        else:
            Al4 = 0

        # M2O DRAW ##
        ipQ4 = (Q4A[-1] * iw[6]) / (
            (Q4A[-1] * iw[6])
            + (Q2A[-1] * iw[6])
            + ((Q1AA[-1] * 2 + Q1AN[-1]) * iw[8])
            + ((Q0AAA[-1] * 3 + Q0AAN[-1] * 2 + Q0ANN[-1]) * iw[9])
            + (Al4A[-1] * w[11])
        )
        ipQ2 = (Q2A[-1] * iw[6]) / (
            (Q4A[-1] * iw[6])
            + (Q2A[-1] * iw[6])
            + ((Q1AA[-1] * 2 + Q1AN[-1]) * iw[8])
            + ((Q0AAA[-1] * 3 + Q0AAN[-1] * 2 + Q0ANN[-1]) * iw[9])
            + (Al4A[-1] * w[11])
        )
        ipQ1 = ((Q1AA[-1] * 2 + Q1AN[-1]) * iw[8]) / (
            (Q4A[-1] * iw[6])
            + (Q2A[-1] * iw[6])
            + ((Q1AA[-1] * 2 + Q1AN[-1]) * iw[8])
            + ((Q0AAA[-1] * 3 + Q0AAN[-1] * 2 + Q0ANN[-1]) * iw[9])
            + (Al4A[-1] * w[11])
        )
        ipQ0 = ((Q0AAA[-1] * 3 + Q0AAN[-1] * 2 + Q0ANN[-1]) * iw[9]) / (
            (Q4A[-1] * iw[6])
            + (Q2A[-1] * iw[6])
            + ((Q1AA[-1] * 2 + Q1AN[-1]) * iw[8])
            + ((Q0AAA[-1] * 3 + Q0AAN[-1] * 2 + Q0ANN[-1]) * iw[9])
            + (Al4A[-1] * w[11])
        )
        ipAl4 = (Al4A[-1] * w[11]) / (
            (Q4A[-1] * iw[6])
            + (Q2A[-1] * iw[6])
            + ((Q1AA[-1] * 2 + Q1AN[-1]) * iw[8])
            + ((Q0AAA[-1] * 3 + Q0AAN[-1] * 2 + Q0ANN[-1]) * iw[9])
            + (Al4A[-1] * w[11])
        )

        Alcoef = 1 / (1 + 3 * ipAl4)

        pQ3 = (Q3[-1] * w[0]) / (
            Q3[-1] * w[0] + Q4 * w[1] + Q2 * w[2] + Q1 * w[3] + Al5[-1] * w[4]
        )
        pQ4 = (Q4 * w[1]) / (
            Q3[-1] * w[0] + Q4 * w[1] + Q2 * w[2] + Q1 * w[3] + Al5[-1] * w[4]
        )
        pQ2 = (Q2 * w[2]) / (
            Q3[-1] * w[0] + Q4 * w[1] + Q2 * w[2] + Q1 * w[3] + Al5[-1] * w[4]
        )
        pQ1 = (Q1 * w[3]) / (
            Q3[-1] * w[0] + Q4 * w[1] + Q2 * w[2] + Q1 * w[3] + Al5[-1] * w[4]
        )
        pAl5 = (
            Alcoef
            * (Al5[-1] * w[4])
            / (Q3[-1] * w[0] + Q4 * w[1] + Q2 * w[2]
               + Q1 * w[3] + Al5[-1] * w[4])
        )

        rpQ3 = pQ3 / (pQ3 + pQ2 + pQ1 + pAl5)
        rpQ2 = pQ2 / (pQ3 + pQ2 + pQ1 + pAl5)
        rpQ1 = pQ1 / (pQ3 + pQ2 + pQ1 + pAl5)
        rpAl5 = pAl5 / (pQ3 + pQ2 + pQ1 + pAl5)

        if Q4 > 0:
            rQ4A = Q4A[-1] / (Q4A[-1] + Q4N[-1])
            rQ4N = Q4N[-1] / (Q4A[-1] + Q4N[-1])
        else:
            rQ4A = 0
            rQ4N = 0

        if Q2 > 0:
            rQ2A = Q2A[-1] / (Q2A[-1] + Q2N[-1])
            rQ2N = Q2N[-1] / (Q2A[-1] + Q2N[-1])
        else:
            rQ2A = 0
            rQ2N = 0

        if Q1 > 0:
            rQ1AA = Q1AA[-1] / (Q1AA[-1] + Q1AN[-1] + Q1NN[-1])
            rQ1AN = Q1AN[-1] / (Q1AA[-1] + Q1AN[-1] + Q1NN[-1])
            rQ1NN = Q1NN[-1] / (Q1AA[-1] + Q1AN[-1] + Q1NN[-1])
        else:
            rQ1AA = 0
            rQ1AN = 0
            rQ1NN = 0

        if Q1AA[-1] + Q1AN[-1] > 0:
            irQ1AA = Q1AA[-1] * 2 / (Q1AA[-1] * 2 + Q1AN[-1])
            irQ1AN = Q1AN[-1] / (Q1AA[-1] * 2 + Q1AN[-1])

        else:
            irQ1AA = 0
            irQ1AN = 0

        if Q0AAA[-1] + Q0AAN[-1] + Q0ANN[-1] > 0:
            irQ0AAA = (Q0AAA[-1] * 3 / (Q0AAA[-1] * 3 +
                                        Q0AAN[-1] * 2 + Q0ANN[-1]))
            irQ0AAN = (Q0AAN[-1] * 2 / (Q0AAA[-1] * 3 +
                                        Q0AAN[-1] * 2 + Q0ANN[-1]))
            irQ0ANN = Q0ANN[-1] / (Q0AAA[-1] * 3 + Q0AAN[-1] * 2 + Q0ANN[-1])

        else:
            irQ0AAA = 0
            irQ0AAN = 0
            irQ0ANN = 0

        # Draws
        if Q4 + Q2 + 2 * Q1 + 3 * Q0 < Q3[0] * w[5]:
            P = 1
        else:
            P = 0

        if (
            Q3[-1]
            - pQ3
            + 3 * pAl5 * (ipQ4 + ipQ2)
            + (-rpQ3 + 3 * rpAl5 * (ipQ4 + ipQ2)) * pQ4
            > 0
        ):
            next_Q3 = (
                Q3[-1]
                - pQ3
                + 3 * pAl5 * (ipQ4 + ipQ2)
                + (-rpQ3 + 3 * rpAl5 * (ipQ4 + ipQ2)) * pQ4
            )
        else:
            next_Q3 = 0

        if Q4N[-1] + pQ3 * P - pQ4 * rQ4N + (rpQ3 * P) * pQ4 > 0:
            next_Q4N = Q4N[-1] + pQ3 * P - pQ4 * rQ4N + (rpQ3 * P) * pQ4
        else:
            next_Q4N = 0

        if (Q4A[-1] - pQ4 * rQ4A - 3 * pAl5 * ipQ4 +
                (-3 * rpAl5 * ipQ4) * pQ4 > 0):
            next_Q4A = (
                Q4A[-1] - pQ4 * rQ4A - 3 * pAl5 * ipQ4 +
                (-3 * rpAl5 * ipQ4) * pQ4
            )
        else:
            next_Q4A = 0

        if (
            Q2A[-1]
            + pQ4 * rQ4A
            - pQ2 * rQ2A
            - 3 * pAl5 * ipQ2
            + 3 * pAl5 * ipQ1 * irQ1AA
            + (-rpQ2 * rQ2A - 3 * rpAl5 * ipQ2 + 3 *
               rpAl5 * ipQ1 * irQ1AA) * pQ4
            > 0
        ):
            next_Q2A = (
                Q2A[-1]
                + pQ4 * rQ4A
                - pQ2 * rQ2A
                - 3 * pAl5 * ipQ2
                + 3 * pAl5 * ipQ1 * irQ1AA
                + (-rpQ2 * rQ2A - 3 * rpAl5 * ipQ2 + 3 * rpAl5 *
                   ipQ1 * irQ1AA) * pQ4
            )
        else:
            next_Q2A = 0

        if (
            Q2N[-1]
            + pQ4 * rQ4N
            - pQ2 * rQ2N
            + pQ3 * (1 - P)
            + 3 * pAl5 * ipQ1 * irQ1AN
            + (-rpQ2 * rQ2N + rpQ3 * (1 - P) + 3 * rpAl5 * ipQ1 * irQ1AN) * pQ4
            > 0
        ):
            next_Q2N = (
                Q2N[-1]
                + pQ4 * rQ4N
                - pQ2 * rQ2N
                + pQ3 * (1 - P)
                + 3 * pAl5 * ipQ1 * irQ1AN
                + (-rpQ2 * rQ2N + rpQ3 * (1 - P) + 3 * rpAl5 *
                   ipQ1 * irQ1AN) * pQ4
            )
        else:
            next_Q2N = 0

        if (
            Q1AA[-1]
            - pQ1 * rQ1AA
            + 3 * pAl5 * ipQ0 * irQ0AAA
            - 3 * pAl5 * ipQ1 * irQ1AA
            + (-rpQ1 * rQ1AA + 3 * rpAl5 * ipQ0 * irQ0AAA - 3 *
               rpAl5 * ipQ1 * irQ1AA)
            * pQ4
            < 0
        ):
            next_Q1AA = 0
        else:
            next_Q1AA = (
                Q1AA[-1]
                - pQ1 * rQ1AA
                + 3 * pAl5 * ipQ0 * irQ0AAA
                - 3 * pAl5 * ipQ1 * irQ1AA
                + (
                    -rpQ1 * rQ1AA
                    + 3 * rpAl5 * ipQ0 * irQ0AAA
                    - 3 * rpAl5 * ipQ1 * irQ1AA
                )
                * pQ4
            )

        if (
            Q1AN[-1]
            + pQ2 * rQ2A
            - pQ1 * rQ1AN
            + 3 * pAl5 * ipQ0 * irQ0AAN
            - 3 * pAl5 * ipQ1 * irQ1AN
            + (
                rpQ2 * rQ2A
                - rpQ1 * rQ1AN
                + 3 * rpAl5 * ipQ0 * irQ0AAN
                - 3 * rpAl5 * ipQ1 * irQ1AN
            )
            * pQ4
            < 0
        ):
            next_Q1AN = 0
        else:
            next_Q1AN = (
                Q1AN[-1]
                + pQ2 * rQ2A
                - pQ1 * rQ1AN
                + 3 * pAl5 * ipQ0 * irQ0AAN
                - 3 * pAl5 * ipQ1 * irQ1AN
                + (
                    rpQ2 * rQ2A
                    - rpQ1 * rQ1AN
                    + 3 * rpAl5 * ipQ0 * irQ0AAN
                    - 3 * rpAl5 * ipQ1 * irQ1AN
                )
                * pQ4
            )

        if (
            Q1NN[-1]
            + pQ2 * rQ2N
            - pQ1 * rQ1NN
            + 3 * pAl5 * ipQ0 * irQ0ANN
            + (rpQ2 * rQ2N - rpQ1 * rQ1NN + 3 * rpAl5 * ipQ0 * irQ0ANN) * pQ4
            < 0
        ):
            next_Q1NN = 0
        else:
            next_Q1NN = (
                Q1NN[-1]
                + pQ2 * rQ2N
                - pQ1 * rQ1NN
                + 3 * pAl5 * ipQ0 * irQ0ANN
                + (rpQ2 * rQ2N - rpQ1 * rQ1NN + 3 *
                   rpAl5 * ipQ0 * irQ0ANN) * pQ4
            )

        if (
            Q0AAA[-1] - 3 * pAl5 * ipQ0 * irQ0AAA + (-3 * rpAl5 *
                                                     ipQ0 * irQ0AAA) * pQ4
            < 0
        ):
            next_Q0AAA = 0
        else:
            next_Q0AAA = (
                Q0AAA[-1]
                - 3 * pAl5 * ipQ0 * irQ0AAA
                + (-3 * rpAl5 * ipQ0 * irQ0AAA) * pQ4
            )

        if (
            Q0AAN[-1]
            + pQ1 * rQ1AA
            - 3 * pAl5 * ipQ0 * irQ0AAN
            + (rpQ1 * rQ1AA - 3 * rpAl5 * ipQ0 * irQ0AAN) * pQ4
            < 0
        ):
            next_Q0AAN = 0
        else:
            next_Q0AAN = (
                Q0AAN[-1]
                + pQ1 * rQ1AA
                - 3 * pAl5 * ipQ0 * irQ0AAN
                + (rpQ1 * rQ1AA - 3 * rpAl5 * ipQ0 * irQ0AAN) * pQ4
            )

        if (
            Q0ANN[-1]
            + pQ1 * rQ1AN
            - 3 * pAl5 * ipQ0 * irQ0ANN
            + (rpQ1 * rQ1AN - 3 * rpAl5 * ipQ0 * irQ0ANN) * pQ4
            < 0
        ):
            next_Q0ANN = 0
        else:
            next_Q0ANN = (
                Q0ANN[-1]
                + pQ1 * rQ1AN
                - 3 * pAl5 * ipQ0 * irQ0ANN
                + (rpQ1 * rQ1AN - 3 * rpAl5 * ipQ0 * irQ0ANN) * pQ4
            )

        if Q0NNN[-1] + pQ1 * rQ1NN + (rpQ1 * rQ1NN) * pQ4 < 0:
            next_Q0NNN = 0
        else:
            next_Q0NNN = Q0NNN[-1] + pQ1 * rQ1NN + (rpQ1 * rQ1NN) * pQ4

        if Al5[-1] - pAl5 - pQ4 * rpAl5 > 0:
            next_Al5 = Al5[-1] - pAl5 - pQ4 * rpAl5
        else:
            next_Al5 = 0

        if Al4A[-1] - 3 * pAl5 * ipAl4 + (-3 * rpAl5 * ipAl4) * pQ4 > 0:
            next_Al4A = (Al4A[-1] - 3 * pAl5 * ipAl4 +
                         (-3 * rpAl5 * ipAl4) * pQ4)
        else:
            next_Al4A = 0

        if (Al4N[-1] + pAl5 + 3 * pAl5 * ipAl4 +
                (rpAl5 + 3 * rpAl5 * ipAl4) * pQ4 > 0):
            next_Al4N = (
                Al4N[-1] + pAl5 + 3 * pAl5 * ipAl4 + (rpAl5 + 3 *
                                                      rpAl5 * ipAl4) * pQ4
            )
        else:
            next_Al4N = 0

        Q3shinanigens.append(
            -3 * pAl5 * ipAl4
            + (-3 * rpAl5 * ipAl4) * pQ4
            + 3 * pAl5 * ipAl4
            + (3 * rpAl5 * ipAl4) * pQ4
        )

        Q3.append(next_Q3)

        Q4A.append(next_Q4A)
        Q4N.append(next_Q4N)

        Q2A.append(next_Q2A)
        Q2N.append(next_Q2N)

        Q1AA.append(next_Q1AA)
        Q1AN.append(next_Q1AN)
        Q1NN.append(next_Q1NN)

        Q0AAA.append(next_Q0AAA)
        Q0AAN.append(next_Q0AAN)
        Q0ANN.append(next_Q0ANN)
        Q0NNN.append(next_Q0NNN)

        Al5.append(next_Al5)

        Al4A.append(next_Al4A)
        Al4N.append(next_Al4N)

        Q3Na.append(next_Q3)
        Al5Na.append(next_Al5)
        Al4Na.append(next_Al4A + next_Al4N)
        Q4Na.append(next_Q4A + next_Q4N)
        Q2Na.append(next_Q2A + next_Q2N)
        Q1Na.append(next_Q1AA + next_Q1AN + next_Q1NN)
        Q0Na.append(next_Q0AAA + next_Q0AAN + next_Q0ANN + next_Q0NNN)

    # Vi laver lister over teoretiske værdier udregnet fra modellen
    mod_m = min(M2O, key=lambda x: abs(x - mod))
    ind = M2O.index(mod_m)
    Q3_m = Q3Na[ind]
    Q4_m = Q4Na[ind]
    Q2_m = Q2Na[ind]
    Q1_m = Q1Na[ind]
    Q0_m = Q0Na[ind]
    Al5_m = Al5Na[ind]
    Al4_m = Al4Na[ind]

    return mod_m, Q3_m, Q4_m, Q2_m, Q1_m, Q0_m, Al5_m, Al4_m


def AlB_SSE(w, data, frac=None, s_plt=False, s_dat=False, p=False):
    mod_data = data[0]
    r_data = data[1]
    B3_data = data[2]
    B4_data = data[3]
    B2_data = data[4]
    B1_data = data[5]
    B0_data = data[6]
    Al5_data = data[7]
    Al4_data = data[8]

    B3_m = []
    B4_m = []
    B2_m = []
    B1_m = []
    B0_m = []
    Al5_m = []
    Al4_m = []

    for i in range(len(mod_data)):
        mod_m, Q3_m, Q4_m, Q2_m, Q1_m, Q0_m, Al5_mod, Al4_mod = AlB_draw(
            w, r_data[i], mod_data[i]
        )

        B3_m.append(Q3_m)
        B4_m.append(Q4_m)
        B2_m.append(Q2_m)
        B1_m.append(Q1_m)
        B0_m.append(Q0_m)
        Al5_m.append(Al5_mod)
        Al4_m.append(Al4_mod)

    if s_plt is False and p is True:
        t = list(range(100))
        plt.plot(
            B4_data,
            B4_m,
            "rd",
            B3_data,
            B3_m,
            "bd",
            B2_data,
            B2_m,
            "gd",
            B1_data,
            B1_m,
            "yd",
            Al4_data,
            Al4_m,
            "kd",
            Al5_data,
            Al5_m,
            "cd",
        )
        plt.plot(t, t, "k--")
        plt.axis([0, 100, 0, 100])
        plt.legend(["$B^4$", "$B^3$", "$B^2$", "$B^1$", "$Al^4$", "$Al^5$"])
        plt.xlabel("Experimental Data")
        plt.ylabel("Model Data")
        plt.title("Quality of fit")
        plt.show()

    if s_plt is True:
        if not os.path.exists("AlB_Structure"):
            os.mkdir("AlB_Structure")
        t = list(range(100))
        plt.plot(
            B4_data,
            B4_m,
            "rd",
            B3_data,
            B3_m,
            "bd",
            B2_data,
            B2_m,
            "gd",
            B1_data,
            B1_m,
            "yd",
            Al4_data,
            Al4_m,
            "kd",
            Al5_data,
            Al5_m,
            "cd",
        )
        plt.plot(t, t, "k--")
        plt.axis([0, 100, 0, 100])
        plt.legend(["$B^4$", "$B^3$", "$B^2$", "$B^1$", "$Al^4$", "$Al^5$"])
        plt.xlabel("Experimental Data")
        plt.ylabel("Model Data")
        plt.title("Quality of fit")
        plt.savefig(os.path.join("AlB_Structure", "Model_v_data.png"))
        plt.show()

    if s_dat is True:
        if not os.path.exists("AlB_Structure"):
            os.mkdir("AlB_Structure")
        m_data_plt_kval = np.column_stack(
            [
                B3_data,
                B3_m,
                B4_data,
                B4_m,
                B2_data,
                B2_m,
                B1_data,
                B1_m,
                Al4_data,
                Al4_m,
                Al5_data,
                Al5_m,
            ]
        )
        np.savetxt(os.path.join("AlB_Structure",
                                "Model_v_data.csv"), m_data_plt_kval)

    if p is False:

        SSE = sum(
            ((B4_data - B4_m) ** 2)
            + ((B3_data - B3_m) ** 2)
            + ((B2_data - B2_m) ** 2)
            + ((B1_data - B1_m) ** 2)
            + ((B0_data - B0_m) ** 2)
            + ((Al5_data - Al5_m) ** 2)
            + ((Al4_data - Al4_m) ** 2)
        )

        return SSE


def AlB_engine(fil, data, it=10):
    dat = data
    w0 = [10, 20, 30, 30, 30, 30, 30, 30, 30, 30, 30]

    minimizer_kwargs = {"method": "COBYLA", "args": (dat,)}
    res = scipy.optimize.basinhopping(
        AlB_SSE,
        w0,
        niter=it,
        T=2.0,
        stepsize=5,
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
