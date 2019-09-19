%function [Corr_earth_vel_cell,Corr_earth_vel_cm,Corr_earth_vel_NaN_added]=tilt_correction(filename,nens)

%This function can apply tilting correction to all the ensembles of a transect. Data source can be any but you have to use the following generic ASCII template for data import. 
%Output Details:
%Corr_earth_vel_cell= Corrected Earth Velocity for all the ensembles in Cell Structure in Centimeters per sec
%Corr_earth_vel_cm= Corrected Earth Velocity array for all the ensembles in Centimeters per sec
%Corr_earth_vel_NaN_added= Corrected Earth Velocity for all the ensembles in Centimeters per sec
                         % with NaN added for using in the WinRiver II default ASCII format in corr_data_gen function.
%w_corr= Corrected vertical velocity
%w_original= Vertical velocity with no tilt correction

%First use the template beam2details_template.ttf in WinRiver II Generic ASCII To export the data of a single transect.
%Then enter the ASCII name and total ensemble number.
%It needs the following functions to run: import_beam,import_ebhpr11,beam2instr,Headtilt and MatMult.
%Example: [Corr_earth_vel_cell,Corr_earth_vel_cm,Corr_earth_vel_NaN_added]=tilt_correction('A2_J_0_000_19-06-09_132346_beam2details_template_ASC.TXT',606);


%Adapted from Code by Bart Vermeulen.
%Transformation matrix reference: ADCP Coordinate Transformation, Teledyne RD Instruments.

%Written by Mohammad Kifayath Chowdhury, Louisiana State University 
%Last Edited: 7/9/2019
filename='ALT_LAT_2_J_0_000_19-06-10_082619_beam2details_template_ASC.TXT';
nens=431;


[beam.Ensemble,beam.nbins,beam.binsize,beam.Heading,beam.Pitch,beam.Roll,beam.Lat,beam.Long]=import_ebhprll(filename, 1, nens);
nCols=max(beam.nbins)*4+3; %finding the size for formatSpec
formatSpec =['%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s' repmat('%f', [1 nCols])];
formatSpec =[formatSpec '%[^\n\r]'];
beam.beamvel=import_beam(filename, 1, nens,formatSpec);
% load('A3_001.mat')
% mat=A3_001;
% beam.beamvel=mat(:,17:683);  %Change 111 according to number of columns in the data 1, 659); 
%Use these 3 lines if column number is more than 119 then convert the file to mat and read it this way.
t_ens=length(beam.Ensemble);
input.conv=1; % 1 for convex set up
input.alpha=20; % Beam Angle in degrees
%% Sorting out all the ensemble beam velocities as matrices

reshaped.beam_vel = cell(t_ens, 1) ; % Cell structure to save the matrices output from the loop
for n=1:t_ens % Desired Ensemble Number, starting from 1
p=beam.beamvel(n,:); % Extract the ensemble from the matrix
input.p1=find(~isnan(p)); % Finding the NaN values (Spaces are converted to NaN during export)
input.e1=reshape(p(input.p1),[],4); % Converting the ensemble to 4 beam columns
input.e1(input.e1==-32768)=NaN; %Replace all bad data with NaN 
reshaped.beam_vel{n}=input.e1; %reshape.beam_vel contains all the beam velocities for all the ensembles of this particular transect
end

%% Conversion of all the beam velocities with the correction to earth velocities
reshaped.cor_earth_vel = cell(t_ens, 1) ; % Cell structure to save the matrices output from the loop
%reshaped.original_earth_vel=cell(t_ens, 1) ;
for n=1:t_ens

input.nbins=beam.nbins(n,1);% Number of bins 
input.nens=1; % Dealing with 1 ensemble only in each iteration
heading=beam.Heading(n,1); % Heading angle of that ensemble
pitch=beam.Pitch(n,1); % Pitch angle of that ensemble
roll=beam.Roll(n,1); % Roll angle of that ensemble
e=reshaped.beam_vel{n,1};
% for t=1:size(e,1)
%          if isnan(e(t,:))
%              e(t,:)=-32768;
%          end
%      for z=1:size(e,2)
%          if isnan(e(t,z))
%          e(t,z)=0;
%          end    
%      end
%          if e(t,:)==-32768
%              e(t,:)=NaN;
%          end
% end
s=zeros(size(e,1),size(e,2));

for t=size(e,1):-1:1
        if isnan(e(t,:))
           s(t,:)=1; 
        else
           s(t,:)=[];
        end   
end 
length_for_resample=beam.nbins(n,1)-size(s,1);
tx=1:length_for_resample;
e(1:length_for_resample,1:4)=resample(e(1:length_for_resample,1:4),tx);%Taking the beam values for the ensembles created in the previous loop
%beam to earth corrected
 TM=beam2instr(input.conv,input.alpha); % transformation from beam 2 instrument coordinates
 TM=MatMult(HeadTilt(heading,pitch,roll,input.nens),TM); % Transformation for rotation including tilt
 tempVEL= cat(1,reshape(e(:,1),1,1,[]),reshape(e(:,2),1,1,[]),reshape(e(:,3),1,1,[]),reshape(e(:,4),1,1,[])); % Reshaping the velocity matrix for multiplication with the transformation matrix
 tempVEL2= MatMult(TM,tempVEL);
 reshaped.cor_earth_vel{n}=cat(2,reshape(tempVEL2(1,1,:),input.nbins,input.nens,1),reshape(tempVEL2(2,1,:),input.nbins,input.nens,1),reshape(tempVEL2(3,1,:),input.nbins,input.nens,1),reshape(tempVEL2(4,1,:),input.nbins,input.nens,1)); %Reshaping the multiplied data into a  matrix format and writing them into a cell structure.

end

Corr_earth_vel_cell=reshaped.cor_earth_vel;
%Corr_earth_vel_cell=round(Corr_earth_vel_cell,0);

for n=1:length(Corr_earth_vel_cell)
    Corr_earth_vel_cell{n,1}=Corr_earth_vel_cell{n,1}.*100; % Converting to cm per sec same as classic ASCII
    Corr_earth_vel_cell{n,1}(isnan(Corr_earth_vel_cell{n,1}))=-32768; % Converting all  the NaN value as -32768, same as Classic ASCII
end

Corr_earth_vel_cm=cell2mat(Corr_earth_vel_cell); %Combining all the matrices from the cell structure

% for n=length(Corr_earth_vel_cm):-1:1
%       if isnan(Corr_earth_vel_cm(n,1))
%         Corr_earth_vel_cm(n,:)=-32768;
%       end
%     %(isnan(Corr_earth_vel_cm(n,1)))=-32768; % Converting all  the NaN value as -32768, same as Classic ASCII
% end
%Additional NaN values are being added to equal the dimension of the
%matrix output from here and the one from excel .csv.

%Create a NaN Matrix
Rep.nanmat_ini=NaN(7,4); % NaN matrix for the initial 7 lines
Rep.nanmat=NaN(6,4); %NaN matrix for the intermediate 6 rows in between each ensemble
%Rep.vertcat=vertcat(Corr_earth_vel_cell{1,1},Rep.nanmat); %Trial to see if it works for one cell
Corr_earth_vel_cell_space_added=Corr_earth_vel_cell; % Creating a new matrix for entering the NaN matrices in between the transects in the for loop
for x=1:length(Corr_earth_vel_cell_space_added)-1
    %Loop to vertcat a 6 by 4 NaN matrix to all the cells
    Corr_earth_vel_cell_space_added{x}=vertcat(Corr_earth_vel_cell_space_added{x,1},Rep.nanmat);
end
Corr_earth_vel_NaN_added=cell2mat(Corr_earth_vel_cell_space_added);
Corr_earth_vel_NaN_added=vertcat(Rep.nanmat_ini,Corr_earth_vel_NaN_added); % vertcat the initial 7 rows. Now a matrix with the same dimension of the classic ascii is made for the 4 velocity columns





%%               **********Functions*********

function import_vel = import_beam(filename, startRow, endRow,formatSpec)
%% Initialize variables.
delimiter = {',',' '};
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Open the text file.
fileID = fopen(filename,'r');

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
import_vel = [dataArray{1:end-1}];
end


function [Ensemble,nbins,binsize,Heading,Pitch,Roll,Lat,Long] = import_ebhprll(filename, startRow, endRow)
%% Initialize variables.
delimiter = {',',' '};
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column3: double (%f)
%   column5: double (%f)
%	column7: double (%f)
%   column9: double (%f)
%	column11: double (%f)
%   column13: double (%f)
%	column15: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%*s%f%*s%f%*s%f%*s%f%*s%f%*s%f%*s%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
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

%% Allocate imported array to column variable names
Ensemble = dataArray{:, 1};
nbins = dataArray{:, 2};
binsize = dataArray{:, 3};
Heading = dataArray{:, 4};
Pitch = dataArray{:, 5};
Roll = dataArray{:, 6};
Lat = dataArray{:, 7};
Long = dataArray{:, 8};
end

function biM=beam2instr(conv,alpha)
%FUNCTION TO COMPUTE TRANSFORMATION MATRIX Beam to instrument
%Determine beam 2 instruments matrix 4 beams adcp's, Ref: page 11 ADCP
%Coordinate Transformation
        a=1/(2*sind(alpha));
        b=1/(4*cosd(alpha));
        d=a/sqrt(2);
        biM=[conv*a,   -conv*a,      0,     0;
            0,      0,   -conv*a,   conv*a;
            b,      b,      b,     b;
            d,      d,     -d,    -d];
end

function htM=HeadTilt(heading, pitch, roll,nens)
% FUNCTION TO COMPUTE ROTATION MATRIX, Ref: Page 18-19 ADCP Coordinate Transformation
% % all angles should be given in degrees 

    pitch=double(pitch)/180*pi;
    roll=double(roll)/180*pi;
pitch=atan(tan(pitch).*cos(roll));    


 heading=double(heading)/180*pi;

sh=reshape(sin(heading),1,1,[]);
ch=reshape(cos(heading),1,1,[]);
sp=reshape(sin(pitch),1,1,[]);
cp=reshape(cos(pitch),1,1,[]);
sr=reshape(sin(roll),1,1,[]);
cr=reshape(cos(roll),1,1,[]);
% heading=double(heading)/180*pi;
% sh=sin(heading);
% ch=cos(heading);
% sp=sin(pitch);
% cp=cos(pitch);
% sr=sin(roll);
% cr=cos(roll);
htM = [ch.*cr+sh.*sp.*sr,  sh.*cp,          ch.*sr-sh.*sp.*cr,  zeros(1,1,nens);
       -sh.*cr+ch.*sp.*sr, ch.*cp,          -sh.*sr-ch.*sp.*cr, zeros(1,1,nens);
       -cp.*sr,            sp,              cp.*cr,             zeros(1,1,nens);
       zeros(1,1,nens),    zeros(1,1,nens), zeros(1,1,nens),    ones(1,1,nens)];

end


function MM=MatMult(A,B)  % Multiply matrices along third dimension
nrows=size(A,1);
ncols=size(B,2);

indim=size(A,2);
if indim~=size(B,1)
    error('MatMult:InDim','Inner matrix dimensions should agree')
end
n3=max(size(A,3),size(B,3));
if (size(A,3)~=size(B,3)) && size(A,3)~=1 && size(B,3)~=1
    error('MatMult:Wrong3dim','Third dimension size must be equal or \n at least one variable must have third dimension size equal to one')
end
MM=zeros(nrows,ncols,n3);
for cntrow=1:nrows
    for cntcol=1:ncols
        for cntdim=1:indim
               MM(cntrow,cntcol,:)=MM(cntrow,cntcol,:)+A(cntrow,cntdim,:).*B(cntdim,cntcol,:);
        end
    end
end
end

