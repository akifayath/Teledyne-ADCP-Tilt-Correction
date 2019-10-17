**** Only for downward looking moving boat ADCPs. 
This package can 
a. Read ADCP velocity in beam coordinates, and backscatter from WinRiver II ASCII files 
b. Apply tilt, roll and pitch correction to the data and give output in earth coordinate.
c. Apply 3 beam solution where necessary (One of the beams failed) 
d. Interpolates data for a missing bin. 
e. Output data as a .dat file for Tecplot.

How to use the correction functions-

1. From WINRIVER II create a classic ASCII output and generic ASCII outputs using template 'beam2details_template.ttf' and 'backscatter.ttf'. Note the total ensemble numbers for each transect from the discharge summary.
2. Open the Classic ASCII in excel using text loader, use 'Space' as delimiter. Copy and save the 4 columns (E,F,G,H) with velocity data as .csv files.
3. In matlab, Open program 'Sample Function Call', give the inputs 

ASCIIname=beam2details generic ascii;
csv_vel=velocity csv file created at step 2;
backscatterASC=backscatter generic ascii;
nens=number of ensembles;
****Change the filename in line 13, it is the velocity output with navigation information for tecplot. Rename the file as .dat and open in tecplot. The dat file contains - Easting, Northing, Depth coordinates, Backscatter Intensity (dB), East Velocity (cm), North Velocity (cm), Vertical Velocity (cm) and error velocity (cm). The tecplot_data matrix can be used for other processing too.

**** There's a threshold limit given in the code to reject any outliers. Run the code with this threshold and check the plots from section "Data observation". If the data seem to be out of the threshold, modify that and run the code again. Also check the plots from the createfit section to see the distribution.
4. Copy the 4 columns from 'Corr_earth_vel_NaN_added' array. Open the excel file of classic ASCII (Step 2), go to E3 and paste.
5. Find and replace all 'NaN' as ' ' (blank, write nothing).
6. Select column A,B and from right click> select format cells>select Number> put 8 in Decimal places> click 'OK'
7. Select column N and right click> select format cells>select Number> put 0 in Decimal places> click 'OK'
8. Save the file as ' filename_ASC.prn '. It is mandatory to name the file as '*****_ASC' for the file to be read by VMT.
9. Change the extension from '.prn' to '.txt'
10. Open in VMT.

****For upward looking ADCP please update the code as following.
i. Go to the function HeadTilt (find it in the bottom) in tilt_correction. 
ii. Change roll as below,
    roll=double(roll)/180*pi+pi; ( It is rotating the ADCP 180 degrees)
This should solve the issue. 

**** If the beam is other than 20 degree, go to tilt_correction function and in line 32 change the input.alpha     
