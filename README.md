1. Introduction

This python package uses statistical mechanics-based modeling
to accurately predict the short-range order structural distribution
in oxide glasses. To make the predictions, the model requires interaction
enthalpies obtained by fitting to experimental data (typically obtained with
NMR spectroscopy techniques).

Despite some already obtained enthalpies, the package is designed for you to
build your own library of enthalpies by providing relevant data to the package.

For a detailed guide, please refer to "the article" (will be filled when published)
Any bugs or questions, do not hesitate to contact msb@bio.aau.dk or mos@bio.aau.dk

2. Package format

The basic format of the package follows:

StatMechGlass
  .
  |-Data
  ||-SiO2
  |||-Na.csv
  |||-Na_Tg.csv
  ||-SiB
  |||-Na.csv
  |-Parameter
  ||-SiO2
  |||-Na.csv
  ||-MF
  |||-SiB.csv
  |-stat_mech_module
  ||-__init__.py
  ||-stat_mech_silicate.py
  ||-stat_mech_borate.py
  |-LICENCE.txt
  |-CITATION.txt
  |-README.md
  |-__init__.py
  |-stat_mech_glass.py

Here, the /Data directory is where you want to place your experimentally
obtained data. The /parameter directory is where the package will automatically
store the enthalpies obtained by fitting the provided data.

3. Usage

When using the package, you may type your commands in the stat_mech_glass.py
file directly. It is advised to write your code within the "if __name__ == '__main__'".
Alternatively, you may place the package in the working directory and import it with

        import StatMechGlass.stat_mech_glass as smg

The four main functions, you will use:

  smg.smg_binary_par(former, modifier, it=10)
  smg.smg_ternary_par(formers, modifier, it=10)
  smg.smg_structure(val, tg, p = None)
  smg.smg_plot(comps, free_comp, tg, plt_save = False)

3.1. Fitting enthalpy parameters on binary oxide glasses

When building the enthalpy database, use smg.smg_binary_par:

  1.  Place your data in the Data/"Former" directory, where "Former" corresponds
      to the network forming specie of the glass. Currently supported formers:
          "SiO2", "B2O3", "P2O5"
      The data file should be named appropriately such as "Na.csv" or "K.csv"

  2.  Place Tg data for the same glass system in the same directory with the name:
          "modifier"_Tg.csv
      Here, "modifier" should be the same as in 1. such as "Na_Tg.csv" or "K_Tg.csv"
      The Tg data does not need to be for the same glass compositions as the
      structural data

  3.  Execute the function
          smg.smg_binary_par(former, modifier, it=10)
      Example:
          smg.smg_binary_par("Si", "Na", it=500)
      500 or more iterations are advised for accurate enthalpies (refer to the
      manuscript for more details)

3.2. Fitting interaction parameters on ternary oxide glasses

When building the former/former interaction database, use smg.smg_ternary_par:

  1.  Place your data in the Data/"Former""Former" directory, where "Former"
      corresponds to the network forming or intermediate specie of the glass.
      Currently supported formers and intermediates:
          SiO2, B2O3, P2O5, Al2O5
      Here, the folder must be named according to the naming convention (section 4):
          "SiB", "BP", "PSi", "AlB"
      The data file should be named appropriately such as "Na.csv" or "K.csv"

  2.  Tg data should be provided for each glass in the data file.
      Refer to section 4 for clarification on datafile content

  3.  Execute the function
          smg.smg_ternary_par(formers, modifier, it=10)
      Example:
          smg.smg_binary_par(["Si", "B"], "Na", it=100)
      100 or more iterations are advised for accurate parameter (refer to the
      manuscript for more details)
