# GEO242_Final_Project
Final project codes, data files, and plots.

## Contents
1. Contreras_final_project_latek (directory containing .tex files, figures, and .bib file, and PDF of final report)
2. HF_empty (directory containing .csv files that list non-data lines in the RE catalogs)
3. Hayward_RE_catalogs (directory of RE_catalog txt files of FARESearch outputs)
4. code_outputs (directory of outputs produced by the codes contianed in this repository)
5. figures (directory of figures produced)
6. NAC_final_project.m (main code; loads data, calculates slip, creep, and produces plots)
7. beta_calib.m (matlab function that finds the least squares fit for the beta constant)
8. creep_models.csv (USGS measurments on the Hayward fault polygon)
9. find_empty_lines.sh (shell script that counts non-data lines int he RE catalogs)
10. grids.txt (txt files of fault polygon grids)
11. plot_fault.gmt6 (gmt script that plots REs, USGS measurement locations, fault traces)
12. slip_creep.m (matlab function that calculates slip and creep rates for a given grid)

## Workflow
1. 

## Directory organization
In order for the INSERTFILENAME.m script to run correctly, ensure that the current working directory mus have the following organization:
- Parent directory
  - Subdirectory: Hayward RE catalogs
    - RE catalog text files (.txt)
  - Subdirectory: Grid workspaces
  - Subdirectory: Final figures
  - 
