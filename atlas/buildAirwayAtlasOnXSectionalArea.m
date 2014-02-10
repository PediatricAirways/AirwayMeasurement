% 
% build pediatric airway atlas on cross sectional area curves
% Author: Yi Hong
% Email: yihong@cs.unc.edu
%


clear all
close all

% input parameters
subject_info = readtext( '/playpen/Project/airwayAtlas/Pediatric_Airway_Atlas_Code/data/selected_subjects_CRL_SGS_Carina.csv' );     % data sheet for measurement
data_prefix = '/playpen/Project/airwayAtlas/Pediatric_Airway_Atlas_Code/results_perimeter';                  % processing results
results_prefix = '/playpen/Project/airwayAtlas/Pediatric_Airway_Atlas_Code/results_perimeter/Atlas24_68_19'; % output location
gaussian_sigma = 24;    % band-width for gaussian function, corresponding to sigma
pos_SGS = 69;           % SGS-pre starting position in the data sheet 
pos_SGS_Post = 75;      % SGS-post starting position in the data sheet, used for classification


% read the subject information from the csv file of data sheet
[nRow, nCol] = size( subject_info );
for iI = 1:nCol
    if strcmp( subject_info{1, iI}, 'PatientId' ) == 1
        idPatient = iI;
    end
    if strcmp( subject_info{1, iI}, 'Age' ) == 1
        idAge = iI;
    end
    if strcmp( subject_info{1, iI}, 'Weight' ) == 1
        idWeight = iI;
    end
end

cases = {};
ages = {};
weights = {};
for icase = 1:nRow-1
    cases{icase} = subject_info{icase+1, idPatient};    % subject id
    ages{icase} = subject_info{icase+1, idAge};         % subject age, used for age-based atlas
    weights{icase} = subject_info{icase+1, idWeight};   % subject weight, used for weight-based atlas
end
ages = cell2mat( ages );
weights = cell2mat( weights );


% plot ages distribution
figure, hold on
tmp = gaussmf( ages, [gaussian_sigma, floor(max(ages)/2)] );
[tmpAges, tmpId] = sort( ages );
tmp = tmp( tmpId );
hist( ages, ceil(max(ages)/10.0)*10 );
title_text = sprintf( 'Histogram of %d subjects'' ages', length(ages) );
title( title_text, 'FontSize', 20, 'FontWeight', 'Bold' );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
xlabel( 'Age: month(s)' );
ylabel( 'Subject number' );

% generate the filenames for centerline, landmarks, and so on
for icase = 1:length(cases)
    tmp = cases{icase};
    if isnumeric(tmp)
        cases{icase} = num2str(tmp);
    end
    tmpStr = sprintf( '%s/IsoSurface%s/%s_MeanAndNormal.txt', data_prefix, cases{icase}, cases{icase} );
    meanNormFile{icase} = tmpStr;
    tmpStr = sprintf( '%s/Contour%s/%s_Area.txt', data_prefix, cases{icase}, cases{icase} );
    areaFile{icase} = tmpStr;
    tmpStr = sprintf( '%s/IsoSurface%s/%s_LandmarksIdOnCenterline.txt', data_prefix, cases{icase}, cases{icase} );
    landmarksFile{icase} = tmpStr;
end

% compute the size of the uniform window for functional boxplot
% make it compariable with the Gaussian window size
% only used control subjects' age because the atlas built on control subjects
ages_CRL = ages( 1:pos_SGS-1 );

% calculate the compariable uniform window 
x = 1:max(ages_CRL);
fx_gaussian = zeros( 1, length(x));
% gaussian density estimation
for iI = 1:length(ages_CRL)
    tmp = gaussmf( x, [gaussian_sigma, ages_CRL(iI)] );
    fx_gaussian = fx_gaussian + tmp;
end
fx_gaussian = fx_gaussian ./ length(ages_CRL);
figure, hist( ages_CRL, ceil(max(ages_CRL)/10.0)*10 );
hold on
plot( x, fx_gaussian, 'r--', 'LineWidth', 2);

% matched with gaussian window
winSizeMatched = 0;
minValueMatched = -1;
for winSize = gaussian_sigma:2*gaussian_sigma
    fx_uniform = zeros( 1, length(x) );
    for iI = 1:length(ages_CRL)
        pickedCurves = ( (x >= ages_CRL(iI) - winSize) .* (x <= ages_CRL(iI) + winSize) );
        fx_uniform = fx_uniform + pickedCurves;
    end
    fx_uniform = fx_uniform ./ length(ages_CRL);
    tmpValue = norm( fx_gaussian-fx_uniform, 2 );
    if minValueMatched < 0 || tmpValue < minValueMatched
        winSizeMatched = winSize;
        minValueMatched = tmpValue;
    end
end

% uniform window density estimation
fx_uniform = zeros( 1, length(x) );
for iI = 1:length(ages_CRL)
    pickedCurves = ( (x >= ages_CRL(iI) - winSizeMatched) .* (x <= ages_CRL(iI) + winSizeMatched) );
    fx_uniform = fx_uniform + pickedCurves;
end
fx_uniform = fx_uniform ./ length(ages_CRL);
plot( x, fx_uniform, 'b-.', 'LineWidth', 2 );
hold off
	

curveAndLandmarksFile = sprintf('%s/multCurvesAndLandmarks.fig', results_prefix);

% atlas building using pointwise, functional, and weighted functional boxplots
params.meanNormFile = meanNormFile;      % centerline 
params.areaFile = areaFile;              % cross-sectional area 
params.landmarksFile = landmarksFile;    % landmarksFile, landmark's id on centerline
params.outputPrefix = results_prefix;    % output location
params.cases = cases;                    % subject id
params.ages = ages;                      % subject age
params.weights = weights;                % subject weight
params.gaussian_sigma = gaussian_sigma;  % gaussian window size
params.uniform_winSize = winSizeMatched; % uniform window size
params.pos_SGS = pos_SGS;                % starting position for pre-SGS
params.pos_SGS_Post = pos_SGS_Post;      % starting position for post-SGS
params.curveAndLandmarksFile = curveAndLandmarksFile;

buildAirwayAtlasForEachSubjectXS( params );


