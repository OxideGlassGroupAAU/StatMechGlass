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
  ├── Data  
  │   ├── SiB  
  │   │   ├── Na.csv  
  │   └── SiO2  
  │       ├── Na.csv  
  │       ├── Na_Tg.csv  
  ├── Parameters  
  │   ├── MF  
  │   │   ├── BAl.csv  
  │   │   └── SiB.csv  
  │   └── SiO2  
  │       ├── K.csv  
  │       └── Na.csv  
  ├── stat_mech_module  
  │   ├── __init__.py  
  │   ├── stat_mech_borate.py  
  │   ├── stat_mech_phosphate.py  
  │   └── stat_mech_silicate.py  
  ├── __init__.py  
  ├── LICENCE  
  ├── README.md  
  ├── stat_mech_glass.py  

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

3.3. Predicting the structural distribution in a given composition

The number of compositions to be predicted by the model increases exponentially
when building the enthalpy database. This is due to the tranferability of the
interaction enthalpies. That is, enthalpies established from Na2O-SiO2 and Na2O-B2O3
glasses may also be used to predict any Na2O-B2O3-SiO2 composition.  

When using the model to predict structural distributions, use smg.smg_structure:

  1.  Define the glass composition using python directory  
      Example:  
          glass_comp = {"Si": 25, "B": 25, "Na": 25, "K":25}  
  2.  Run the function with a defined tg  
      Example:  
          results = smg.smg_structure(glass_comp, 700)

This way, users may easily build a database of structures from a large set of
glass compositions

3.4. Simple 2D plotting

As glass compositions may consist of many different elements, smg.smg_structure
can be used to return structures which the user may plot themselves. For simple
visualization, the smg.smg_plot function may be used:

  1.  Define the glass composition using python directory  
      Example:  
          glass_comp = {"Si": 25, "B": 25, "Na": 0}  
      Alternatively, leave out the free component:  
          glass_comp = {"Si": 25, "B": 25}  
  2.  Run the function with a defined tg and free component  
      Example:  
          smg.smg_plot(glass_comp, "Na", 800, plt_save = True)  
      Set plt_save to True for saving the plot as .png file

4. Naming convention

For the script to locate files, parameters and functions, a certain naming
convention was introduced which must be followed:

4.1. Network formers and modifiers

When calling network formers in any function these are abbreviated to their
basic atom:  
"Si", "B", "Al", "P"  
Example:  
25SiO2-25B2O3-50Na2O should be {"Si": 25, "B": 25, "Na": 50}. Note that the
glass contains 2 boron atoms for each silicon atom but in the naming convention
they seem to contain the same number of atoms.

When using modifiers in the functions, these should be named according to the data
files provided by the user. If the datafile is named "Na.csv", "Na" should be used
in the functions.

4.2 Datafiles

4.2.1 Binary oxide glass data
