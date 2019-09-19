# Teledyne_ADCP_Tilt_Correction_and_Tecplot_Data_Creator
This set of functions can be used to correct tilting errors for a teledyne RDI ADCP. 
How to use the correction functions-

1. From WINRIVER II create a classic ASCII output and generic ASCII outputs using template 'beam2details_template.ttf' and 'backscatter.ttf'. Note the total ensemble number.
2. Open the Classic ASCII in excel using text loader, use 'Space' as delimiter. Copy and save the 4 columns (E,F,G,H) with velocity data as .csv files
3. In matlab, Open program 'Sample Function Call', give the inputs 

ASCIIname=beam2details generic ascii;
csv_vel=velocity csv file created at step 2;
backscatterASC=backscatter generic ascii;
nens=number of ensembles;

Change the limits after checking the plots from the createfit section.
4. Copy the 4 columns from 'Corr_earth_vel_NaN_added' array. Open the excel file of classic ASCII, go to E3 and paste.
5. Find and replace all 'NaN' as ' ' (blank, write nothing).
6. Select column A,B and from right click> select format cells>select Number> put 8 in Decimal places> click 'OK'
7. Select column N and right click> select format cells>select Number> put 0 in Decimal places> click 'OK'
8. Save the file as ' filename_ASC.prn ' .
9. Change the extension from '.prn' to '.txt'
10. Open in VMT.
11. Also in step 3, change the filename in line 13, it is the velocity output with navigation information for tecplot. Rename the file as .dat and open in tecplot.    
