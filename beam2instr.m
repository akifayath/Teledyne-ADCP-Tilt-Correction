function biM=beam2instr(conv,alpha)
%FUNCTION TO COMPUTE TRANSFORMATION MATRIX Beam to instrument
%Determine beam 2 instruments matrix 4 beams adcp's, Ref: page 11 ADCP
% Function created by Bart Vermeulen.

%Coordinate Transformation
        a=1/(2*sind(alpha));
        b=1/(4*cosd(alpha));
        d=a/sqrt(2);
        biM=[conv*a,   -conv*a,      0,     0;
            0,      0,   -conv*a,   conv*a;
            b,      b,      b,     b;
            d,      d,     -d,    -d];
end
