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
