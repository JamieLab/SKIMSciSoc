# Shelf break flows calculation



# Instructions
To run the analysis: the surface current data needs to be downloaded and reformatted, shelf break identified and the direction on to shelf determined, currents calculated with respect to the direction onto shelf and the mean data extracted.
Analysis to be run as follows:
1. Download GEBCO bathymetry data from: https://www.gebco.net/data_and_products/gridded_bathymetry_data/
2. Run DATA_retrieve_inputs.py which downloads all the current data - FTP username and password need to be added to downloadekmansurf(), downloadekmandepth() and downloadgeostrophic().
3. Run DATA_calculate mean_monthly_currents.py
4. Run analysis.py

These steps will output a csv file with mean annual, seaonal currents for 14 shelf seas in Laruelle et al. (2021; https://doi.org/10.1038/s41467-017-02738-z) for further analysis (example found in output_xlsx/table2data_global_500m.csv). 
This file must be manually split into a annual and 4 seasonal xlsx files (examples in output_xlsx) and the dpCO2 trends from Laruelle et al. (2021; https://doi.org/10.1038/s41467-017-02738-z) added manually.

5. Run PLOT_means_data.py to produce final results.
6. Run plot_current_along_shelf.py to produce plots of different currents for shelf regions
7. Run plot_current_component_dominance.py to produce global maps of the dominance of each current component
8. Run plot_method_figs.py to produce the method example figure
9. Run PLOT_supplementary_S2.py to produce the supplementary figure with different shelf break depths and the resulting total across shelf break currents.

# Supporting manuscript
J.D. Shutler, D.J. Ford, T. Holding, C. Ubelmann, L. Gaultier, F. Collard, B. Chapron, M. Rio, C. Donlon and C. Roberts. Shelf sea carbon accumulation rates are consitent with variability in wind driven shelf break exchange. (in prep)

# Funding
This code and analysis was funded by European Space Agency (ESA) Sea surface KInematics Multiscale monitoring (SKIM) Mission Science study, the ESA SKIM Scientific Performance Evaluation study and the Convex Seascape Survey (https://convexseascapesurvey.com/)
