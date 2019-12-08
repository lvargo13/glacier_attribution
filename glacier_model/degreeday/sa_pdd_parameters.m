
global CONFIG

% temp and radiation factors
%CONFIG.DegreeDay.RadiationFactor = 0.13 ;
%CONFIG.DegreeDay.DDF = 1.5 ;

% climate offsets (based on brewster data- p x1.3, t -1.25)  
% CONFIG.DegreeDay.PptnFactor = 0.8 ;   
% CONFIG.DegreeDay.TempOffset = -0.6; 

CONFIG.DegreeDay.AlbedoIce = 0.35; % cuffey and patterson
CONFIG.DegreeDay.AlbedoSnow = 0.85;
CONFIG.DegreeDay.AlbedoCoef = 0.112;  % Brock 2000
CONFIG.DegreeDay.SnowTempThreshold = 1; 

% snow th
CONFIG.DegreeDay.SnowMaxThickness = 0;   % max snow thickness 
CONFIG.DegreeDay.SnowLineElevation = 1940;   % (m a.s.l.) elevation of snowline at start of model
CONFIG.DegreeDay.SnowGradient = 1; % (mm/m) (gradient of snow to use above initial snowline
CONFIG.DegreeDay.SLThreshMin = 15; % minimum snow depth (mm) to calc ELA
CONFIG.DegreeDay.SLThreshMax = 150; % max snow depth (mm) to calc ELA
