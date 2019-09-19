function [A,original_vel]=corr_data_gen(ASCIIname,filename2,nens,emax,emin,nmax,nmin,wmax,wmin)

%This function can output the tilt corrected data for a transect that can be copied into the classic ASCII from WINRIVER II. 
%Output Details:

%Corr_earth_vel_cell= Corrected Earth Velocity for all the ensembles in Cell Structure in Centi meters per sec
%Corr_earth_vel_cm= Corrected Earth Velocity for all the ensembles in Centi meters per sec
%Corr_earth_vel_NaN_added= Corrected Earth Velocity for all the ensembles in Centi meters per sec
                         % with NaN added for using in the WinRiver II default ASCII format.

                        
%filename1= filename of the ASCII txt file from modified generic ASCII template
%filename2= filename of the 4 column csv


%filename 2 fact-----
%Import the csv file of the 4 columns that contain the earth velocity (Read
%the ASCII as text in Excel, give space as delimiter, change nothing else.
%Copy the 4 velocity columns and save as .csv).
%First two blank rows of the ASCII do not come to matlab. 

%nens is the total number of ensembles from the transect

%It needs the following functions to run- 
%tilt_correction, import_csv_ASCII, createFit( for data comparison)
%Example: 
% [Corr_earth_vel_cell,Corr_earth_vel_cm,Corr_earth_vel_NaN_added,w_corr,w_original,original_vel]=corr_data_gen('A2_J_0_001_19-06-09_133106_beam2details_template_ASC.TXT','A2_J_0_001_19-06-09_133106_ASC.csv',659);

%Prepared by
%Mohammad Kifayath Chowdhury, Louisiana State University
%Last Edit: 7/9/2019


[Corr_earth_vel_cell]=tilt_correction(ASCIIname,nens);
for n=1:length(Corr_earth_vel_cell)
    Corr_earth_vel_cell{n,1}=Corr_earth_vel_cell{n,1}.*100; % Converting to cm per sec same as classic ASCII
end

%emax=170;emin=-170;nmax=150;nmin=-150;wmax=40;wmin=-40;
%% Applying Threshold
for q=1:length(Corr_earth_vel_cell)
    e=Corr_earth_vel_cell{q,1};
    if ~any(~isnan(e(:)))
        e=nan(size(e,1),size(e,2));% If all the beams failed, no need to apply the thresholds
    else    
        e_corr=e(:,1);
        n_corr=e(:,2);
        w_corr=e(:,3);
        err_corr=e(:,4);

        e_corr(e_corr<emin)=NaN;
        e_corr(e_corr>emax)=NaN;
        n_corr(n_corr<nmin)=NaN;
        n_corr(n_corr>nmax)=NaN;
        w_corr(w_corr<wmin)=NaN; 
        w_corr(w_corr>wmax)=NaN; 
        e=horzcat(e_corr,n_corr,w_corr,err_corr);
    errordata=find (all(isnan(e),1)); %if any of the components is way over the distribution nd threshold then all the columns might end up NaNs.
                                      %it means error level is very high
                                      %and the ensemble is to be ignored
          if isempty(errordata) 
              x=1; %No big error
          elseif (length(errordata))==2||length(errordata)==3||length(errordata)==4||length(errordata)==1%if all the beams failed together
              x=111; %Data is out of distribution
          end 
          
      if x==1
            s=zeros(size(e,1),size(e,2));
              for  t=size(e,1):-1:1
                  if isnan(e(t,:))
                       s(t,:)=1; 
                  else
                       s(t,:)=[];
                  end   
               end 
            length_for_resample=size(e,1)-size(s,1);
            if length_for_resample>1
                tx=1:length_for_resample; % Resampling all the columns for possible missing or bad value in any of the bins
                e(1:length_for_resample,1:4)=resample(e(1:length_for_resample,1:4),tx);
            end
      elseif x==111
             e=nan(size(e,1),size(e,2));
      end
      
    end
    up=e(:,3);
    up(up<-40)=NaN;
    e(:,3)=up;    
    Corr_earth_vel_cell{q}=e;
    
end
%%
for n=1:length(Corr_earth_vel_cell)
       Corr_earth_vel_cell{n,1}(isnan(Corr_earth_vel_cell{n,1}))=-32768; % Converting all  the NaN value as -32768, same as Classic ASCII
end


Corr_earth_vel_cm=cell2mat(Corr_earth_vel_cell); %Combining all the matrices from the cell structure
%Create a NaN Matrix
Rep.nanmat_ini=NaN(7,4); % NaN matrix for the initial 7 lines
Rep.nanmat=NaN(6,4); %NaN matrix for the intermediate 6 rows in between each ensemble
%Rep.vertcat=vertcat(Corr_earth_vel_cell{1,1},Rep.nanmat); %Trial to see if it works for one cell
Corr_earth_vel_cell_space_added=Corr_earth_vel_cell; % Creating a new matrix for entering the NaN matrices in between the transects in the for loop
    for x=1:length(Corr_earth_vel_cell_space_added)-1
        %Loop to vertcat a 6 by 4 NaN matrix to all the cells
        Corr_earth_vel_cell_space_added{x}=vertcat(Corr_earth_vel_cell_space_added{x,1},Rep.nanmat);
    end
Corr_earth_vel_cell_space_added=vertcat(Rep.nanmat_ini,Corr_earth_vel_cell_space_added);% vertcat the initial 7 rows. Now a matrix with the same dimension of the classic ascii is made for the 4 velocity columns
Corr_earth_vel_NaN_added=cell2mat(Corr_earth_vel_cell_space_added);

%%
format short
Corr_earth_vel_NaN_added=round(Corr_earth_vel_NaN_added,1);
endrow=length(Corr_earth_vel_cm)+6*(nens-1)+7; %Calculate endrow from existing data
original_vel= import_csv_ASCII(filename2, 1, endrow); %read 4 columns .csv directly
for t=1:length(Corr_earth_vel_NaN_added)
    if isnan(Corr_earth_vel_NaN_added(t,1)) 
        Corr_earth_vel_NaN_added(t,1)=original_vel(t,1);
        original_vel(t,1)=NaN;
    end
    if isnan(Corr_earth_vel_NaN_added(t,2)) 
        Corr_earth_vel_NaN_added(t,2)=original_vel(t,2); 
        original_vel(t,2)=NaN;
    end
    if isnan(Corr_earth_vel_NaN_added(t,3)) 
        Corr_earth_vel_NaN_added(t,3)=original_vel(t,3);
        original_vel(t,3)=NaN;
    end
    if isnan(Corr_earth_vel_NaN_added(t,4)) 
        Corr_earth_vel_NaN_added(t,4)=original_vel(t,4);   
    end
end
A.Corr_earth_vel_cell=Corr_earth_vel_cell;
A.Corr_earth_vel_cm=Corr_earth_vel_cm;
A.Corr_earth_vel_NaN_added=Corr_earth_vel_NaN_added;
    
A.w_original=original_vel(:,3);
A.w_original(isnan(A.w_original))=[];
A.w_original(A.w_original==-32768)=NaN;

A.w_corr=Corr_earth_vel_cm(:,3);
A.w_corr(A.w_corr==-32768)=NaN;

A.e_original=original_vel(:,1);
A.e_original(isnan(A.e_original))=[];
A.e_original(A.e_original==-32768)=NaN;

A.e_corr=Corr_earth_vel_cm(:,1);
A.e_corr(A.e_corr==-32768)=NaN;

A.n_original=original_vel(:,2);
A.n_original(isnan(A.n_original))=[];
A.n_original(A.n_original==-32768)=NaN;

A.n_corr=Corr_earth_vel_cm(:,2);
A.n_corr(A.n_corr==-32768)=NaN;
end
% createFit(w_original,w_corr)

%Now Copy 'Corr_earth_vel_NaN_added' matrix to excel, replace the data in the
%designated location,keep first two rows blank Find and replace the NaN as blank and save the edited ASCII as .prn
%Rename the extension as .txt
%The text file is now readable in VMT.

%% Functions
function A = import_csv_ASCII(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   A = import_csv_ASCII(FILENAME) Reads data from text file
%   FILENAME for the default selection.
%
%   A = import_csv_ASCII(filename, startRow, endRow)Reads
%   data from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   A = import_csv_ASCII('A2_J_0_000_19-06-09_132346_ASC.csv', 1, 14640);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2019/07/07 16:27:22

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r','n','UTF-8');
% Skip the BOM (Byte Order Mark).
fseek(fileID, 3, 'bof');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
A = [dataArray{1:end-1}];
end