# Generalised-Multistep-Dynamic-TOPMODEL version 2 (v2)
The main difference with v1 is that v2 uses a new vertical hydrulic conductivity profile and therefore requires one less parameter, i.e., Hmax (max subsurface depth), which was a calibrated parameter in v1, is no longer necessary. Also, some bugs have been fixed and some functions have been made more efficient. Aside from these changes everything else is kept the same as v1.


Steps to run the model:
-First make sure you download the latet release: click on the green "<> code" button and press "Download ZIP".
-Download the "Generalised-Multistep-Dynamic_TOPMODEL-main.zip" and unzip it. 
-Open "Generalised-Multistep-Dynamic_TOPMODEL-main" and unzip "DATA.zip".
-Note that this will create a path with double DATA folders: ...\DATA\DATA\...
-You will need to copy-paste the contents of the 2nd DATA folder into the 1st DATA folder and then remove the 2nd DATA folder.
-Your working directory path (pwd) should look something like this:
-C\...\Generalised-Multistep-Dynamic_TOPMODEL-main\DATA
-Once you are done unziping, open and run "GMD_TOPMODEL_run.m".
