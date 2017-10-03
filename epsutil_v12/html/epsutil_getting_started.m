%% Upslope Area Toolbox: Getting Started
%
%% Acknowledgements
% The Milford (milford_ma_dem) and Natick (natick_ned) sample data sets are
% courtesy of the <http://www.usgs.gov U.S. Geological Survey (USGS)>.
%
%% What are Digital Elevation Models?
% A _digital elevation model_ (or DEM) is a computer representation of surface
% topography.  A _raster DEM_ is a rectilinear grid of values, each of which
% represents the height of a surface at the corresponding grid location.
%
% High-quality, high-resolution DEMs are now widely available and being used for
% a wide variety of terrain analysis.  DEM data for the United States can be
% obtained through the U.S. Geological Survey (USGS) and its data providers.  
%
% The file milford_ma_dem.mat contains an example DEM covering a portion of
% Massachusetts in the United States.  You can load this MAT-file and display
% the elevation matrix |Z| as follows: 

load milford_ma_dem
imshow(Z, [])  % imshow is in the Image Processing Toolbox

%% What Does "Upslope Area" Mean?
% Imagine you are standing on the side of hill somewhere in the rain. Some of
% the water that falls uphill from your position will flow directly toward and
% then past your shoes.  Some of the water, though, will flow downhill in a
% different direction, away from you.  The area of land above you that drains
% directly through where you are standing is called the _upslope area_ of your
% position.
%
% If you were standing at the very top of the hill, the upslope area there would
% be 0; no water flows to you from anywhere else.  On the other hand, if you
% stood at the deepest point in a crater with high rims all the way around, the
% upslope area would be the entire area of the crater.
%
% Upslope area is an important hydrology measurement used to study water
% drainage networks, the motion of sediments and contaminants, erosion,
% landslides.
%
%% What Can Upslope Area Toolbox Do?
% The Upslope Area Toolbox can:
%
% * Compute the upslope area for every location in a DEM.
% * Compute the downhill water drainage region (influence map) for a given
% location in a DEM.
% * Compute the uphill water drainage region (dependence map) for a given
% location in a DEM.
% * Fill interior sinks in a DEM.
% * Visualize upslope area, influence maps, and dependence maps superimposed on
% the original DEM.
% * Identify fill values connected to the edges of a DEM.
%
%% What Data Can I Use?
% DEM data for the United States can be obtained from the U.S. Geological Survey
% (USGS).  For example, DEM data and other datasets can be obtained from the
% <http://seamless.usgs.gov/index.php National Map Seamless Server>. You can
% download this data in BIL format, which can be read using the MATLAB function
% |multibandread|.  There are several web sites that offer information and
% tutorials on getting data from the Seamless Server, including
% <http://www.yale.edu/ceo/Documentation/dem_import.pdf this one at Yale
% University>.
%
% For locating DEM data covering other regions, you might try the listing of
% <http://www.terrainmap.com/rm39.html "Free Digital Elevation Model (DEM) and
% Free Satellite Imagery Download Links"> at <http://www.terrainmap.com/
% terrainmap.com>.
%
% If you have a recent version (R2009b or later) of Mapping Toolbox, you can use
% the Web Map Service (WMS) features to obtain DEM data.  The 
% <http://www.mathworks.com/access/helpdesk/help/toolbox/map/ref/wmsread.html
% |wmsread| reference page> has an example showing how to obtain DEM data
% directly from the JPL WMS server.
%
% For sources of other geospatial data that can be read using functions in the
% Mapping Toolbox, see <http://www.mathworks.com/support/tech-notes/2100/2101.html 
% MathWorks Technical Note 2101 -  Accessing Geospatial Data on the Internet for
% the Mapping Toolbox>. 
%
%% Example: Where Does Runoff Go when it Leaves my Property? 
% The sample data file natick_ned.bil is a 1/3-Arc Second National Elevation
% Dataset downloaded using the <http://seamless.usgs.gov/ USGS Seamless Server>.
% The MathWorks headquarters in Natick, Massachusetts is located at the center
% of this data set. In this example you will determine where water runoff goes
% when it leaves The MathWorks. 
%
% The file natick_ned_output_parameters.txt, included in the download, contains
% enough information to read the data file using the MATLAB function
% |multibandread|.

type natick_ned_output_parameters.txt

%%
% Read the DEM data using |multibandread| and display it using |imshow|.

Z = multibandread('natick_ned.bil', [360 403 1], 'int16', 0, 'bil', 'ieee-le');
imshow(Z, [])

%%
% Calculate the flow matrix for the DEM.

R = demFlow(Z);
T = flowMatrix(Z, R);

%%
% Calculate the influence map for the location in the center of the DEM.
% Visualize the influence map using |visMap|.

[M, N] = size(Z);
r = round(M/2);
c = round(N/2);
I = influenceMap(Z, T, r, c);
visMap(I, Z, r, c)

%%
% Zoom in.

axis([195 295 110 200])

%%
% You can see that water draining from The MathWorks headquarters flows
% generally to the northeast.

%%
% Copyright 2009 The MathWorks, Inc.
