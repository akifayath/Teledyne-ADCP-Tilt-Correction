function [fitresult, gof] = createFit1(w,name)
%CREATEFIT1(W)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      Y Output: w
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 13-Aug-2019 11:52:02


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( [], w );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name', name );
h = plot( fitresult, xData, yData );
legend( h, name, 'fit 1', 'Location', 'NorthEast' );
% Label axes
ylabel (name)
grid on


