% Only for downward looking moving boat ADCPs. 
% This code can 
% a. Read ADCP velocity in beam coordinates, and backscatter from WinRiver II ASCII files. 
% b. Apply tilt, roll and pitch correction to the data and give output in earth coordinate.
% c. Apply 3 beam solution where necessary (One of the beams failed) 
% d. Interpolates data for a missing bin. 
% e. Output data as a .dat file for Tecplot.
%Output Details:

%A.Corr_earth_vel_cell= Corrected Earth Velocity for all the ensembles in Cell Structure in Centi meters per sec
%A.Corr_earth_vel_cm= Corrected Earth Velocity for all the ensembles in Centi meters per sec
%A.Corr_earth_vel_NaN_added= Corrected Earth Velocity for all the ensembles in Centimeters per sec
                         % with NaN added equalizng it to the WinRiver II default ASCII format.
%backscatter= Backscatter intensity for the transect
%tecplot_data= Corrected data ready to plot in tecplot
                         
%Adapted from Code by Bart Vermeulen.
%Transformation matrix reference: ADCP Coordinate Transformation, Teledyne RD Instruments.

%Written by Mohammad Kifayath Chowdhury, Louisiana State University 
%Last Edited: 7/9/2019

clear all

ASCIIname='ALT_LAT_2_J_0_004_19-06-10_085523_beam2details_template_ASC.TXT';
csv_vel='ALT_LAT_6_J_0_000.csv';
backscatterASC='ALT_LAT_2_J_0_004_19-06-10_085523_backscatter_ASC.TXT';
nens=446;
[A,original_vel]=corr_data_gen(ASCIIname,csv_vel,nens,170,-170,150,-150,40,-40);
backscatter = backscatter(ASCIIname,backscatterASC, nens);
tecplot_data=horzcat(backscatter,A.Corr_earth_vel_cm);
tecplot_data(isnan(tecplot_data))=-32768;
save('ALT_LAT_6_June','tecplot_data','-ascii')
%%
figure
createFit(A.w_original,A.w_corr)

figure
createFit(A.e_original,A.e_corr)


figure
createFit(A.n_original,A.n_corr)
%% Data observation
e=A.e_original;
createFit1(e,'Easting')
n=A.n_original;
title('Easting')
createFit1(n,'Northing')
title('Northing')
w=A.w_original;
createFit1(w,'Vertical')
title('Vertical')