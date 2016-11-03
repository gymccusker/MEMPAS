# MEMPAS
Manchester Electron Microscope Particle Analysis Software

Authors: Gillian Young (University of Manchester)

Contact gillian.young@manchester.ac.uk with any queries.

-----------------------------------------------------------------------------

MATLAB software used to classify elemental composition from scanning electron
microscopes into atmospherically-relevant aerosol particle species.

-----------------------------------------------------------------------------

1. Installation instructions
2. Instructions for use
3. Version control
4. Planned updates
5. Licensing information

-----------------------------------------------------------------------------

1. Installation instructions

  All functions necessary for use are packaged within the MEMPAS_V* directory 
  (V* depending on version number). Place this directory in a location with 
  write access, and add this location to the MATLAB path to enable use. 

  NOTE: MEMPAS was written in MATLAB R2013a, if any issues arise from differing 
  MATLAB versions, please let the authors know.

2. Instructions for use

  **FOR OLD SEM**

  Run sem_load from the MATLAB command line and follow the instructions in the 
  MATLAB workspace. Make sure the data to be analysed is in the current working 
  directory.

  **FOR NEW SEM**

  Load in function is a work in progress due to format of new SEM output files.
  For now, the following instructions allow this software to be used with the
  new SEM files:

  - Save output file (<filename>.xlsx) as a .csv file **of the same name**.
    This can be done in Microsoft Excel, File -> Save As -> <filename>.csv.

  - **In the .csv version**, copy columns B (Field No), F (Area), G (Average 
    Diameter), and L (Volume) to the end of the file, after the last column 
    filled. See example_file.csv for the required format. **DO NOT ALTER THE
    .XLSX VERSION OF THE FILE**

  - Run sem_load from the MATLAB command line. **Choose the .xlsx version of the
    file to load in.** Follow the instructions in the MATLAB workspace. Make sure
    the data to be analysed is in the current working directory.

3. Version control

  Version control is in place with MEMPAS. Updates and identified bugs must be 
  submitted to the branch on GitHub. Authors of successful updates will be added
  to the author list.

4. Planned updates

  - sem_load will be updated to avoid issues with loading in new SEM data.
  - Switches will be updated to a GUI in a subsequent version. 
  - Possible addition of cluster analysis.

5. Licensing information

  This software is free to use for any interested parties, with the requirement 
  that it, and the authors, are acknowledged in any resulting publications and 
  the journal articles used (Young et al., 2016; Kandler et al., 2011) are 
  acknowledged where appropriate.

  - Kandler, K., Lieke, K., Benker, N., Emmel, C., Küpper, M., Müller-Ebert, 
  D., Ebert, M., Scheuvens, D., Schladitz, A., Schütz, L., and Weinbruch, S.:
  Electron microscopy of particles collected at Praia, Cape Verde, during the 
  Saharan Mineral Dust Experiment: particle chemistry, shape, mixing 	state and 
  complex refractive index, Tellus Series B, 63, 475–496, 	
  doi:10.1111/j.1600-0889.2011.00550.x, 2011.

  - Young, G., Jones, H. M., Darbyshire, E., Baustian, K. J., McQuaid, J. B., 
  Bower, K. N., Connolly, P. J., Gallagher, M. W., and Choularton, T. W.: 
  Size-segregated compositional analysis of aerosol particles collected in the 
  European Arctic during the ACCACIA campaign, Atmospheric Chemistry and Physics, 
  16, 4063–4079, doi:10.5194/acp-16-4063-2016, 2016.
