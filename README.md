# GEO242_Final_Project
Final project codes, data files, and plots.

## Contents
1. Contreras_final_project_latek 
    - (directory containing .tex files, figures, and .bib file, and PDF of final report)
2. HF_empty 
    - (directory containing .csv files that list non-data lines in the RE catalogs)
3. Hayward_RE_catalogs 
    -(directory of RE_catalog txt files of FARESearch outputs)
4. code_outputs 
    - (directory of outputs produced by the codes contianed in this repository)
5. figures 
    - (directory of figures produced)
6.  HF_creep_model.txt 
    - (txt files of USGS measurements on fault polygon)
7. NAC_final_project.m 
    - (main code; loads data, calculates slip, creep, and produces plots)
8. beta_calib.m 
    - (matlab function that finds the least squares fit for the beta constant)
9. creep_models.csv 
    - (USGS surface creep rate measurments)
10. find_empty_lines.sh 
    - (shell script that counts non-data lines int he RE catalogs)
11. grids.txt 
    - (txt files of fault polygon grids)
12. plot_fault.gmt6 
    - (gmt script that plots REs, USGS measurement locations, fault traces)
13. q_faults.tar.gz 
    - (tar file containing kml files of relevant quaternary faults)
14. slip_creep.m 
   - (matlab function that calculates slip and creep rates for a given grid)

## Workflow
1. run 'find_empty_lines.sh'
    - must be in the same directory as the RE catalogs
    - produces empty lines csv files used to separate families (must copy into HF_empty directory or change NAC_final_project.m file load path)
2 run NAC_final_project.m
    - slip_creep.m, beta_calib.m , Hayward_RE_catalogs, grids.txt, empty lines .csv files, and creep_models.csv ust be in the same directory
    - produces figures for slip, creep rates, constant caibrations, HF_creep_model.txt (USGS measurements on the fault polygon), short- and long-term creep rates txt files, and .mat files used to load data ack in for each grid
    - produces .mat files used to load data back in for each grid
3. Run plot_fault.gmt
    - grids.txt, RE catalog text files (.txt), HF_creep_model.txt, untarred q_faults.tar.gz contents must be in the same directory
    - produces plot of fault (creep_locs.ps & creep_locs.jpg)
