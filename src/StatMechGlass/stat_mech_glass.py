#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 14:40:28 2020

@author: M.S. Bødker

"""

import matplotlib.pyplot as plt
import csv
import numpy as np
import os
import math
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import scipy.optimize


if "StatMechGlass" in os.getcwd():
    import stat_mech_module as smm
else:
    from . import stat_mech_module as smm


def _data_load(path, file_name, col_nr):
    """
    This function will load the data required by the other functions.
    It takes a path, file name and the required column as inputs
    """
    current_dir = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    data_col = []
    with open(os.path.join(path, f"{file_name}.csv"), newline="") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=",", quotechar="|")
        for row in spamreader:
            data_col.append(row[col_nr])
    os.chdir(current_dir)
    data_col = [float(i) for i in data_col]
    data_col = np.array(data_col)
    return data_col


def _form_lookup(former, modifier=None, path_in=None):
    """
    This function is used to lookup the correct terms and formers when
    running the functions
    """
    fil = modifier
    dat = 0
    if former == "Si":
        path = "Parameters/SiO2"
        path2 = "Data/SiO2"
        engine_fun, SSE_fun, draw_fun = (
            smm.stat_mech_silicate.Si_engine,
            smm.stat_mech_silicate.Si_SSE,
            smm.stat_mech_silicate.Si_onedraw,
        )
        if modifier:
            if path_in:
                dat = (
                    _data_load(path_in, fil, 0),
                    _data_load(path_in, fil, 1),
                    _data_load(path_in, fil, 2),
                    _data_load(path_in, fil, 3),
                    _data_load(path_in, fil, 4),
                    _data_load(path_in, fil, 5),
                )
            else:
                dat = (
                    _data_load(path2, fil, 0),
                    _data_load(path2, fil, 1),
                    _data_load(path2, fil, 2),
                    _data_load(path2, fil, 3),
                    _data_load(path2, fil, 4),
                    _data_load(path2, fil, 5),
                )
        s_conc = {"Si4": 100, "Si3": 0, "Si2": 0, "Si1": 0, "Si0": 0}
        weight = ["wSi4", "wSi3", "wSi2", "wSi1"]
        atom_frac = 2
        data_q = ["Si4", "Si3", "Si2", "Si1", "Si0"]
        first_draw = None

    elif former == "B":
        path = "Parameters/B2O3"
        path2 = "Data/B2O3"
        engine_fun, SSE_fun, draw_fun = (
            smm.stat_mech_borate.B_engine,
            smm.stat_mech_borate.B_SSE,
            smm.stat_mech_borate.B_onedraw,
        )
        if modifier:
            if path_in:
                dat = (_data_load(path_in, fil, 0), _data_load(path_in, fil, 1))
            else:
                dat = (_data_load(path2, fil, 0), _data_load(path2, fil, 1))
        s_conc = {"B3": 100, "B4": 0, "B2": 0, "B1": 0, "B0": 0}
        weight = ["wb3", "wb4", "wb2", "wb1"]
        atom_frac = 1
        data_q = ["B4"]
        first_draw = None

    elif former == "Al":
        path = "Parameters/Al2O3"
        path2 = "Data/Al2O3B2O3"
        engine_fun, SSE_fun, draw_fun, first_draw = (
            smm.stat_mech_aluminoborate.AlB_engine,
            smm.stat_mech_aluminoborate.AlB_SSE,
            smm.stat_mech_aluminoborate.AlB_one_draw,
            smm.stat_mech_aluminoborate.AlB_first_draw,
        )
        if modifier:
            if path_in:
                dat = (_data_load(path_in, fil, 0), _data_load(path_in, fil, 1))
            else:
                dat = (_data_load(path2, fil, 0), _data_load(path2, fil, 1))
        s_conc = {"Al6": 100, "Al4": 0}
        weight = ["wal6"]
        atom_frac = 1
        data_q = ["Al4", "Al6"]

    elif former == "P":
        path = "Parameters/P2O5"
        path2 = "Data/P2O5"
        engine_fun, SSE_fun, draw_fun = (
            smm.stat_mech_phosphate.P_engine,
            smm.stat_mech_phosphate.P_SSE,
            smm.stat_mech_phosphate.P_onedraw,
        )
        if modifier:
            if path_in:
                dat = (
                    _data_load(path_in, fil, 0),
                    _data_load(path_in, fil, 1),
                    _data_load(path_in, fil, 2),
                    _data_load(path_in, fil, 3),
                    _data_load(path_in, fil, 4),
                )
            else:
                dat = (
                    _data_load(path2, fil, 0),
                    _data_load(path2, fil, 1),
                    _data_load(path2, fil, 2),
                    _data_load(path2, fil, 3),
                    _data_load(path2, fil, 4),
                )
        s_conc = {"P3": 100, "P2": 0, "P1": 0, "P0": 0}
        weight = ["wp3", "wp2", "wp1"]
        atom_frac = 1
        data_q = ["p3", "p2", "p1", "p0"]
        first_draw = None

    elif former == "AlB":
        path = "Data/Al2O3B2O3"
        engine_fun, SSE_fun = (
            smm.stat_mech_aluminoborate.AlB_engine,
            smm.stat_mech_aluminoborate.AlB_SSE,
        )
        if modifier:
            if path_in:
                dat = (
                    _data_load(path_in, fil, 0),
                    _data_load(path_in, fil, 3),
                    _data_load(path_in, fil, 4),
                    _data_load(path_in, fil, 5),
                    _data_load(path_in, fil, 6),
                    _data_load(path_in, fil, 7),
                    _data_load(path_in, fil, 8),
                    _data_load(path_in, fil, 9),
                    _data_load(path_in, fil, 10),
                )
            else:
                
                dat = (
                    _data_load(path, fil, 0),
                    _data_load(path, fil, 3),
                    _data_load(path, fil, 4),
                    _data_load(path, fil, 5),
                    _data_load(path, fil, 6),
                    _data_load(path, fil, 7),
                    _data_load(path, fil, 8),
                    _data_load(path, fil, 9),
                    _data_load(path, fil, 10),
                )
        first_draw = None

    return (
        path,
        engine_fun,
        SSE_fun,
        draw_fun,
        s_conc,
        weight,
        atom_frac,
        dat,
        path2,
        data_q,
        first_draw,
    )


def _tg_fit(tg_data, mod):
    """
    This function takes Tg and modifier data to fit the Tg values as 
    a function of modifier concentration.
    """
    x = np.array(tg_data[0])
    x = x.reshape(-1, 1)
    y = np.array(tg_data[1])
    polynomial_features = PolynomialFeatures(degree=3)
    x_poly = polynomial_features.fit_transform(x)
    lin = LinearRegression()
    lin.fit(x_poly, y)

    tg = []
    mod = np.array(mod)
    for m in range(len(mod)):
        next_tg = lin.predict(polynomial_features.fit_transform(
                                mod[m].reshape(1, -1)))
        tg.append(next_tg)
    tg = np.array(tg)
    return tg


def smg_structure(val, tg, p=None):
    """
       This function will calculate the structural distribution of any glass
       composition. The function requires accurate relative reaction enthalpies
       for all possible chemical interactions in the glass melt.

    =============================================================================
       smg_structure(val, tg, p = None)
    =============================================================================

       where val is the chemical composition of the desired glass.
       val should be a python dictionary in the form: {"Si":25,"B":25,"Na":50}.
       Please refer to the README file for elaboration on the naming convention

       tg is the temperature atoms in the glass-forming liquid
       stops rearranging due to the kinetic barrier.
       This is assumed to be equal to the fictive temperature of the glass.
       The higher the tg, the more disorder in the structural distribution.

       p is only used when the function is called by the
       ternary glass parameter optimization.
       This parameter should not be altered manually


       Example:

       >>> res_structures = smg_structure({"Si":25, "B": 25, "Na":50}, tg=700)
    """
    comp = {
        "formers": ["Si", "B", "P"],
        "intermediates": ["Al"],
        "modifiers": ["Na", "K", "Li", "Ca"],
    }

    formers_s = comp["formers"]
    formers = []
    f_conc = []
    intermediates_s = comp["intermediates"]
    intermediates = []
    i_conc = []
    modifiers_s = comp["modifiers"]
    modifiers = []
    m_conc = []

    for i in range(len(formers_s)):
        a_frac = _form_lookup(formers_s[i])[6]
        try:
            f_conc.append(val[formers_s[i]] / a_frac)
        except:
            pass
        else:
            formers.append(formers_s[i])

    for i in range(len(intermediates_s)):
        try:
            a_frac = _form_lookup(intermediates_s[i])[6]
            i_conc.append(val[intermediates_s[i]] / a_frac)
        except:
            pass

    if sum(i_conc) > sum(f_conc):
        print(
            "The results may be inaccurate since the concentration of"
            "intermediates is higher than the concentration of formers"
            )
    for i in range(len(intermediates_s)):
        a_frac = _form_lookup(intermediates_s[i])[6]
        try:
            f_conc.append(val[intermediates_s[i]] / a_frac)
        except:
            pass
        else:
            intermediates.append(intermediates_s[i])

    for i in range(len(modifiers_s)):
        try:
            m_conc.append(val[modifiers_s[i]])
        except:
            pass
        else:
            modifiers.append(modifiers_s[i])

    t_n_draws = (sum(m_conc) / sum(f_conc)) * 100

    n_draws = int(t_n_draws)

    if t_n_draws - n_draws > 0.5:
        n_draws += 1

    # Starting concentrations:
    structures = {}
    for i in range(len(formers)):
        start_conc = _form_lookup(formers[i])[4]
        for i2 in start_conc:
            structures[i2] = start_conc[i2] * f_conc[i] / sum(f_conc)

    for i in range(len(intermediates)):
        start_conc = _form_lookup(intermediates[i])[4]
        for i2 in start_conc:
            structures[i2] = (
                            start_conc[i2] *
                            f_conc[i + len(formers)] / sum(f_conc)
                            )

    if len(intermediates) > 0:
        structure_val = []
        for i in structures:
            structure_val.append(structures[i])

        H_int = []

        for i in formers:
            w_data = list(_data_load(_form_lookup(i)[0], intermediates[0], 0))
            for i2 in range(len(w_data)):
                H_int.append(w_data[i2])

        for i in intermediates:
            w_data = list(_data_load(_form_lookup(i)[0], intermediates[0], 0))
            for i2 in range(len(w_data)):
                H_int.append(w_data[i2])

        w_int = []

        for i in range(len(H_int)):
            w_int.append(math.exp(-H_int[i] / (tg * 0.008314)))

        structure_alb = _form_lookup(intermediates[0])[10](
            w_int, structure_val, formers[0]
        )

        indi = 0
        for i in structures.keys():
            structures[i] = structure_alb[indi]
            indi += 1
    else:
        pass

    # Defining weighting factors
    for m in range(n_draws):
        weights = {}
        draws = []

        if len(intermediates) > 0:
            for i in range(len(intermediates)):
                if m == 0:
                    formers.append(intermediates[i])
                else:
                    pass
            for i in range(len(modifiers)):
                m_weights = {}
                m_draws = []
                for i2 in range(len(formers)):
                    if i2 == 0:
                        f = 1
                    else:
                        if p:
                            f = p
                        else:
                            f_p = formers[0] + formers[i2]
                            f_path = "Parameters/MF"
                            f = _data_load(f_path, f_p, 0)[0]
                    path = _form_lookup(formers[i2])[0]
                    w_names = _form_lookup(formers[i2])[5]
                    Hi = _data_load(path, modifiers[i], 0)

                    m_weights[w_names[0]] = 1 * f

                    struc_keys = list(_form_lookup(formers[i2])[4].keys())
                    m_draws.append(m_weights[w_names[0]] *
                                   structures[struc_keys[0]])

                    if len(Hi) > 1:
                        for i3 in range(len(Hi)):
                            m_weights[w_names[i3 + 1]] = (
                                math.exp(-Hi[i3] / (tg * 0.008314))
                            ) * f
                            m_draws[i2] += (
                                m_weights[w_names[i3 + 1]]
                                * structures[struc_keys[i3 + 1]]
                            )
                    else:
                        pass

                for i4 in m_weights:
                    try:
                        weights[i4] += (
                                        m_weights[i4] *
                                        (m_conc[i] / sum(m_conc))
                                        )
                    except:
                        weights[i4] = (
                                        m_weights[i4] *
                                        (m_conc[i] / sum(m_conc))
                                        )

                for i5 in range(len(m_draws)):
                    try:
                        draws[i5] += m_draws[i5]
                    except:
                        draws.append(m_draws[i5])

            draws_norm = []
            for i in range(len(draws)):
                draws_norm.append((draws[i] / sum(draws)))

            for i in range(len(formers)):
                step_conc_ind = list(_form_lookup(formers[i])[4].keys())
                step_conc = []

                for i2 in step_conc_ind:
                    step_conc.append(structures[i2])
                step_w_ind = _form_lookup(formers[i])[5]
                step_w = []
                for i2 in step_w_ind:
                    step_w.append(weights[i2])

                onedraw_fun = _form_lookup(formers[i])[3]

                new_conc = list(onedraw_fun(step_w, step_conc, draws_norm[i]))

                for i3 in range(len(step_conc_ind)):
                    structures[step_conc_ind[i3]] = new_conc[i3]

                if formers[i] in intermediates:
                    back_draw = -draws_norm[i] * 3

                    for i in range(len(formers)):
                        if formers[i] not in intermediates:

                            step_conc_ind = list(_form_lookup(
                                            formers[i])[4].keys()
                                                            )
                            step_conc = []

                            for i2 in step_conc_ind:
                                step_conc.append(structures[i2])

                            step_w_ind = _form_lookup(formers[i])[5]
                            step_w = []
                            for i2 in step_w_ind:
                                step_w.append(weights[i2])

                            onedraw_fun = _form_lookup(formers[i])[3]
                            new_conc = list(
                                onedraw_fun(step_w, step_conc,
                                            back_draw, back=True)
                            )

                            for i3 in range(len(step_conc_ind)):
                                structures[step_conc_ind[i3]] = new_conc[i3]

        else:

            for i in range(len(modifiers)):
                m_weights = {}
                m_draws = []
                for i2 in range(len(formers)):
                    if i2 == 0:
                        f = 1
                    else:
                        if p:
                            f = p
                        else:
                            f_p = formers[0] + formers[i2]
                            f_path = "Parameters/MF"
                            f = _data_load(f_path, f_p, 0)[0]
                    path = _form_lookup(formers[i2])[0]
                    w_names = _form_lookup(formers[i2])[5]
                    Hi = _data_load(path, modifiers[i], 0)

                    m_weights[w_names[0]] = 1 * f

                    struc_keys = list(_form_lookup(formers[i2])[4].keys())
                    m_draws.append(m_weights[w_names[0]] *
                                   structures[struc_keys[0]])

                    for i3 in range(len(Hi)):
                        m_weights[w_names[i3 + 1]] = (
                            math.exp(-Hi[i3] / (tg * 0.008314))
                        ) * f
                        m_draws[i2] += (
                            m_weights[w_names[i3 + 1]] *
                            structures[struc_keys[i3 + 1]]
                        )

                for i4 in m_weights:
                    try:
                        weights[i4] += (m_weights[i4] *
                                        (m_conc[i] / sum(m_conc))
                                        )
                    except:
                        weights[i4] = m_weights[i4] * (m_conc[i] / sum(m_conc))

                for i5 in range(len(m_draws)):
                    try:
                        draws[i5] += m_draws[i5]
                    except:
                        draws.append(m_draws[i5])

            draws_norm = []
            for i in range(len(draws)):
                draws_norm.append((draws[i] / sum(draws)))

            for i in range(len(formers)):
                step_conc_ind = list(_form_lookup(formers[i])[4].keys())
                step_conc = []
                for i2 in step_conc_ind:
                    step_conc.append(structures[i2])

                step_w_ind = _form_lookup(formers[i])[5]
                step_w = []
                for i2 in step_w_ind:
                    step_w.append(weights[i2])

                onedraw_fun = _form_lookup(formers[i])[3]

                new_conc = list(onedraw_fun(step_w, step_conc, draws_norm[i]))

                for i3 in range(len(step_conc_ind)):
                    structures[step_conc_ind[i3]] = new_conc[i3]

    return structures


def smg_basin_binary(former, modifier, it=10, path_in=None):
    """
       This function will calculate interaction
       enthalpies for binary oxide glasses.
       If you wish to automatically save the parameter to your
       /Parameter directory, please refer to the smg_binary_par function
       This function requires quantitative experimental data for
       structural distribution and fictive temperature.
       Refer to README for details of data file formatting

    =============================================================================
       smg_basin_binary(former, modifier, it=10, path_in=None)
    =============================================================================

       where former and modifier are string parameters such as "Si" and "Na".
       Please refer to the README file for elaboration on the naming convention

       it is the number of iterations the basinhopping parameter optimazation
       function will run. Please refer to the manuscript for elaboration
       
       path_in can be set to any desired path, where the files will be found
       instead of the default path

       The function requires structural data in the /Data directory
       under the directory with the same name as the desired former.
       In the sodium silicate example, a Na.csv file should be placed in
       the Data/SiO2 directory. Additionally, tg data should
       be provided in a seperate Na_tg.csv file in the same directory

       Example:

       >>> optimal_parameters = smg_basin_binary("Si", "Na", it=500)
    """
    draw_nr = list(range(400))
    draw_ar = np.array(draw_nr)

    M2O = []

    a_frac = _form_lookup(former)[6]

    for i in draw_ar:
        next_mod = (
                    (draw_ar[i] / a_frac) / (100 +
                                             ((draw_ar[i] / a_frac))) * 100
                    )
        M2O.append(next_mod)
    
    if path_in:
        path = path_in
    else:
        path = _form_lookup(former)[8]
    tg_des = modifier + "_Tg"
    tg_data = (_data_load(path, tg_des, 0), _data_load(path, tg_des, 1))

    tg = _tg_fit(tg_data, M2O)

    fil = modifier
    dat, engine_fun, SSE_fun = (
        _form_lookup(former, fil, path_in)[7],
        _form_lookup(former, fil, path_in)[1],
        _form_lookup(former, fil, path_in)[2],
    )

    par = engine_fun(fil, dat, tg, it)
    SSE_fun(par, dat, tg, frac=None, s_plt=False, s_dat=False, p=True)

    return par


def smg_binary_par(former, modifier, it=10, path_in=None):
    """
       This function will calculate and save interaction enthalpies for binary
       oxide glasses. If you don't wish to automatically save the parameter to
       your /Parameter directory, please refer to the smg_basin_binary function
       This function requires quantitative experimental data for structural
       distribution and fictive temperature. Refer to README for details of
       data file formatting

    =============================================================================
       smg_binary_par(former, modifier, it=10)
    =============================================================================

       where former and modifier are string parameters such as "Si" and "Na".
       Please refer to the README file for elaboration on the naming convention

       it is the number of iterations the basinhopping parameter optimazation
       function will run. Please refer to the manuscript for elaboration
       
       path_in can be set to any desired path, where the files will be found
       instead of the default path

       The function requires structural data in the /Data directory under the
       directory with the same name as the desired former. In the sodium
       silicate example, a Na.csv file should be placed in the Data/SiO2
       directory. Additionally, tg data should be provided in a seperate
       Na_tg.csv file in the same directory

       Example:

       >>> smg_binary_par("Si", "Na", it=500)
    """

    par = smg_basin_binary(former, modifier, it)

    path = _form_lookup(former, path_in)[0]

    current_dir = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    np.savetxt(os.path.join(path, "{}.csv".format(modifier)), par)

    os.chdir(current_dir)

    return print("Parameters {} saved to {} in {}".format(par, modifier, path))


def _smg_ternary_SSE(p, formers, modifier):
    """
    This function will return the SSE when
    fitting ternary oxide glass parameters.
    """

    data_path = "Data/" + formers[0] + formers[1]

    m_data, f1_data, f2_data, tg_data = (
        _data_load(data_path, modifier, 0),
        _data_load(data_path, modifier, 1),
        _data_load(data_path, modifier, 2),
        _data_load(data_path, modifier, 3),
    )

    SSE = 0

    for i in range(len(m_data)):
        values = {formers[0]: f1_data[i],
                  formers[1]: f2_data[i], modifier: m_data[i]}
        tg = tg_data[i]
        res_struc = smg_structure(values, tg, p)
        data_ind = 4
        for i2 in range(len(formers)):
            data_q = _form_lookup(formers[i2])[9]
            for i3 in range(len(data_q)):
                data = _data_load(data_path, modifier, data_ind)[i3]
                SSE += (data - float(res_struc[data_q[i3]])) ** 2
                data_ind += 1

    return SSE


def smg_ternary_p_opt(formers, modifier, it=10):
    """
       This function will fit former/former interactions for ternary oxide
       glasses. If you wish to automatically save the parameter to your
       /Parameter directory, please refer to the smg_ternary_par function
       This function requires quantitative experimental data for structural
       distribution and fictive temperature.
       Refer to README for details of data file formatting

    =============================================================================
       smg_ternary_p_opt(formers, modifier, it=10)
    =============================================================================

       where formers is a list of strings such as ["B", "Si"] and modifier
       such as "Na". Please refer to the README file for elaboration on the
       naming convention.

       it is the number of iterations the basinhopping parameter optimazation
       function will run. Please refer to the manuscript for elaboration

       The function requires structural data in the /Data directory under the
       directory with the same name as the desired formers. In the sodium
       borosilicate example, a Na.csv file should be placed in the Data/BSi
       directory. Additionally, tg data should be provided in the same file.

       Example:

       >>> optimal_parameters = smg_ternary_p_opt(["B", "Si"], "Na", it=500)
    """

    w0 = 1
    minimizer_kwargs = {
        "method": "COBYLA",
        "args": (
            formers,
            modifier,
        ),
    }
    res = scipy.optimize.basinhopping(
        _smg_ternary_SSE,
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


def smg_ternary_par(formers, modifier, it=10):
    """
       This function will fit former/former interactions for ternary oxide
       glasses. If you don't wish to automatically save the parameter to your
       /Parameter directory, please refer to the smg_ternary_p_opt function.
       This function requires quantitative experimental data for structural
       distribution and fictive temperature. Refer to README for details of
       data file formatting

    =============================================================================
       smg_ternary_par(formers, modifier, it=10)
    =============================================================================

       where formers is a list of strings such as ["B", "Si"] and modifier
       such as "Na". Please refer to the README file for elaboration on the
       naming convention.

       it is the number of iterations the basinhopping parameter optimazation
       function will run. Please refer to the manuscript for elaboration

       The function requires structural data in the /Data directory under the
       directory with the same name as the desired formers. In the sodium
       borosilicate example, a Na.csv file should be placed in the Data/BSi
       directory. Additionally, tg data should be provided in the same file.

       Example:

       >>> smg_ternary_par(["B", "Si"], "Na", it=500)
    """

    par = float(smg_ternary_p_opt(formers, modifier, it))
    par_in = 1 / par
    path = "Parameters/MF/"

    name1 = formers[0] + formers[1]
    name2 = formers[1] + formers[0]

    current_dir = os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    with open("{}{}.csv".format(path, name1), "w") as f:
        f.write(str(par))
        f.close()

    with open("{}{}.csv".format(path, name2), "w") as f:
        f.write(str(par_in))
        f.close()

    os.chdir(current_dir)

    print("Parameter {} saved to {} in {}".format(par, name1, path))
    print("Parameter {} saved to {} in {}".format(par_in, name2, path))


def smg_plot(comps, free_comp, tg, plt_save=False):
    """
       This function will plot the structures obtained by the smg_structure
       function. The function takes the composition of the glass, the variable
       component and a glass transition temperature.

    =============================================================================
       smg_plot(comps, free_comp, tg, plt_save = False):
    =============================================================================

       where comps is the chemical composition of the desired glass.
       comps should be a python dictionary in the
       form: {"Si":50, "B":50, "Na":0}.
       Please refer to the README file for elaboration on the naming convention

       free_comp is free component which will be plotted on the x-axis.
       It should be a string such as "Na"

       tg is the glass transition temperature used for calculating the
       structures for plotting.

       plt_save will determine whether or not the plot will
       be saved in png format


       Example:

       >>> smg_plot({"Si":50, "B":50, "Na":0}, "Na", 800, plt_save = True)
    """

    for i in range(101):
        comps[free_comp] = i
        structures = smg_structure(comps, tg)
        if i == 0:
            structures_end = structures
            for key in structures_end:
                structures_end[key] = [structures_end[key]]
        else:
            for key in structures:
                structures_end[key].append(structures[key])

    plt_legend = []
    for key in structures_end:
        plt.plot(range(101), structures_end[key])
        plt_legend.append(key)
    plt.legend(plt_legend)
    plt.xlabel("x")
    plt.ylabel("Structure species concentration")
    if plt_save:
        plt_name = ""
        for key in comps:
            plt_name += key
        plt.savefig(("{}_plot.png".format(plt_name)))
    plt.show()

