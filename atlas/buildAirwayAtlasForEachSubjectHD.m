% 
% build atlas for each subject based on cross-sectional area curves

function buildAirwayAtlasForEachSubjectHD( params )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: read files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read ceneterline and compute dist in 1D, from 0 to 1
numf = max( size( params.meanNormFile ) );
numP = zeros(numf, 1);
for n = 1:numf
    disp(params.meanNormFile{n});
    fid = fopen( params.meanNormFile{n}, 'r');
    tline = fgets(fid);
    count = 0;
    while ~feof(fid)
        tline = fgets(fid);
        count = count + 1;
        tmp = sscanf(tline, '%f');
        A(:, count, n) = tmp(1:3);
        normCases(:, count, n) = tmp(4:6);
    end
    fclose(fid);
    numP(n) = count;
    dist(1,n) = 0;
    for i=2:count
        p1_x = A(1,i,n); p1_y = A(2,i,n); p1_z = A(3,i,n);
        p2_x = A(1,i-1,n); p2_y = A(2,i-1,n); p2_z = A(3,i-1,n);
        dist(i,n) = sqrt( (p1_x-p2_x)^2 + (p1_y-p2_y)^2 + (p1_z-p2_z)^2 ) + dist(i-1,n);
    end
    %dist_copy( 1:count, n ) = dist( 1:count, n );
    dist(1:count,n) = dist(1:count,n) ./ dist(count,n);
end

% read landmarks id on centerline
for n = 1:numf
    fidLandmarksIdFile = fopen( params.landmarksFile{n}, 'r' );
    tline = fgets( fidLandmarksIdFile );
    nLandmarksId = sscanf( tline, '%d' );
    nInc = 0;
    for iI = 1:nLandmarksId
        tline = fgets( fidLandmarksIdFile );
        landmarksIdOnCenterline = sscanf( tline, '%d' );
        % remove the subglottic landmark
        if nLandmarksId > 5 && iI == 5 
            nInc = nInc + 1;
            continue;
        end
        LandmarksId( n, iI-nInc ) = landmarksIdOnCenterline;
        Landmarks(n, iI-nInc) = dist( landmarksIdOnCenterline, n);
    end
end

% read area file
for n = 1:numf
    fid = fopen(params.areaFile{n}, 'r');
    tline = fgets(fid);
    nNumArea = sscanf(tline, '%d');
    for i = 1:nNumArea
        tline = fgets(fid);
        if(tline == -1) break; end
        area(i,n) = sscanf(tline, '%f');
    end
    fclose(fid);
end

% read perimeter file, use the area variable to store the hydraulicDiameter, 
% so that we don't need to change the following code
% a little ugly here
for n = 1:numf
    fid = fopen(params.perimeterFile{n}, 'r');
    tline = fgets(fid);
    nNumArea = sscanf(tline, '%d');
    for i = 1:nNumArea
        tline = fgets(fid);
        if(tline == -1) break; end
        perimeter(i,n) = sscanf(tline, '%f');
		% compute hydraulic diameter
		if perimeter(i,n) < 1e-7
			area(i,n) = 0;
			disp(params.perimeterFile{n});
		else
			area(i,n) = 4*area(i,n)/perimeter(i,n);
		end
    end
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curve registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ('./Matlabfunctions/fdaM')
addpath ('./FunctionalBoxplot')
%  set up a basis for the functions W(t) that define the warping functions
%  create curves to fit the area 
for n = 1:numf
    rng      = [0,1];
    knots    = dist(1:numP(n),n)';
    norder   = 6;
    nbasis   = length(knots) + norder - 2;
    areabasis = create_bspline_basis(rng, nbasis, norder, knots);
    Lfdobj   = int2Lfd(2);
    lambda   = 1e-5;   % adjust this number to make the curve fit the scattered point                                   
    areafdPar = fdPar(areabasis, Lfdobj, lambda);

    % curve fitting function for the area
    % try to make the first and last points on the curve.
    xValue = dist(1:numP(n), n);
    xValue(1) = xValue(1) + 100 * 1e-7;
    xValue(end) = xValue(end) - 100 * 1e-7;
    yValue = area(1:numP(n), n);
    for iTmp = 1:100
        xValue = [ xValue(1) - 1e-7; xValue ; xValue(end) + 1e-7 ];
        yValue = [ yValue(1); yValue; yValue(end) ];
    end
    areafd_tmp = smooth_basis( xValue, yValue, areafdPar );
    % to make the nasal spine to TVC smooth
    % and TVC to trachea carina fit points well
    areaSmooth(1:LandmarksId(n, 4)-1, n) = eval_fd(dist(1:LandmarksId(n, 4)-1, n), areafd_tmp );
    areaSmooth(LandmarksId(n, 4):numP(n), n) = area(LandmarksId(n, 4):numP(n), n);
end

% make all curves represented by the same number of points
agefine = linspace(0,1,501)';
for n = 1:numf
    rng      = [0,1];
    knots    = dist(1:numP(n),n)';
    norder   = 6;
    nbasis   = length(knots) + norder - 2;
    areabasis = create_bspline_basis(rng, nbasis, norder, knots);
    Lfdobj   = int2Lfd(2);
    lambda   = 0.5 * 1e-6;  % adjust this number to make the curve fit the scattered point
    areafdPar = fdPar(areabasis, Lfdobj, lambda);

    % curve fitting function for the area
    % try to make the first and last points on the curve.
    xValue = dist(1:numP(n), n);
    xValue(1) = xValue(1) + 100 * 1e-7;
    xValue(end) = xValue(end) - 100 * 1e-7;
    yValue = areaSmooth(1:numP(n), n);
    for iTmp = 1:100
        xValue = [ xValue(1) - 1e-7; xValue ; xValue(end) + 1e-7 ];
        yValue = [ yValue(1); yValue; yValue(end) ];
    end
    areafd = smooth_basis( xValue, yValue, areafdPar );
    
    clear yValue   
    areafine(:,n) = eval_fd(agefine, areafd);
end

% curve presentation
rng_new = [0,1];
knots_new = agefine;
norder_new = 6;
nbasis_new = length(knots_new) + norder_new - 2;
basis_new = create_bspline_basis(rng_new, nbasis_new, norder_new, knots_new);
Lfdobj   = int2Lfd(2);
lambda   = 0.5 * 1e-6;
fdPar_new = fdPar(basis_new, Lfdobj, lambda);
areafd_new = smooth_basis(knots_new, areafine, fdPar_new);

% find the value for landmarks on curves
for icase = 1:numf
    valueLM(icase, :) = eval_fd( Landmarks(icase, :), areafd_new(icase));
end

% draw curve and landmarks
style = {'yx', 'rx', 'gx', 'bx', 'cx'};
areafine_new = eval_fd(agefine, areafd_new);

% the color to display different ages
[cmap, nBins] = generateColormap( params.ages );
colormap(cmap);

% normalized weight for coloring
minWeight = 0;
maxWeight = ceil(max( params.weights ));
colorPosW = round( ( params.weights - minWeight ) / ( maxWeight - minWeight + 1e-10 ) * nBins + 1 );
colorPosW = min( colorPosW, nBins );
colorPosW = max( colorPosW, 1 );

% plot all the unregistered curves
figure('Colormap', cmap), hold on
for icase = 1:numf
    if icase < params.pos_SGS
    	plot( agefine, areafine_new( :, icase ), 'Color', cmap(round(params.ages(icase)), :), 'LineWidth', 1 );
    else
        plot( agefine, areafine_new( :, icase ), 'Color', cmap(round(params.ages(icase)), :), 'LineWidth', 2, 'LineStyle', '--' );
    end
end
% display landmarks on curves
for icase = 1:size(Landmarks, 2)
     plot( Landmarks(:, icase), valueLM(:, icase), char(style(mod(icase, 5)+1)), 'LineWidth', 3);
end
hold off
caxis( [0 nBins] );
h = colorbar( 'peer', gca );
set( get( h, 'ylabel' ), 'String', 'Age: month(s)', 'FontSize', 20, 'FontWeight', 'Bold' );
set( h, 'ylim', [0, nBins] );
set( gca, 'XTick', 0:0.1:1 );
xlabel( 'Depth along the centerline', 'FontSize', 20, 'FontWeight', 'Bold' );
ylabel( 'Cross-sectional area (mm^2)', 'FontSize', 20, 'FontWeight', 'Bold' );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Unregistered Curves', 'FontSize', 20, 'FontWeight', 'Bold' );
saveas( gca, params.curveAndLandmarksFile );

% do the rigistration with landmarks
Landmarks( :, 1 ) = Landmarks( :, 1 ) + 1e-7;
Landmarks( :, end ) = Landmarks( :, end ) - 1e-7;
nLandmarks = size( Landmarks, 2 );
 
wbasisLM = create_bspline_basis([0,1], max(nLandmarks+3-2,4), 3);
WfdLM    = fd(zeros(max(nLandmarks+3-2,4),1),wbasisLM);
WfdParLM = fdPar(WfdLM,1,1e-12);
 
%  carry out the landmark registration
LandmarksMean = mean(Landmarks);
[areafdLM, areawarpfdLM, WfdLM] = landmarkreg(areafd_new, Landmarks, LandmarksMean, WfdParLM, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot registeration results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
areamatUR = eval_fd(agefine, areafd_new);
areamatLM = eval_fd(agefine, areafdLM);
areawarpmatLM  = eval_fd(agefine, areawarpfdLM);
areawarpmatLM(1,:) = 0; areawarpmatLM(length(agefine),:) = 1;

% plot the warping functions
figure('Colormap', cmap), hold on
for icase = 1:numf
    plot( agefine, areawarpmatLM( :, icase ), 'Color', cmap(round(params.ages(icase)), :), 'LineWidth', 1 );
end
for icase = 1:numf
    for iLM = 1:size( Landmarks, 2 )
	plot(mean(Landmarks(:,iLM)), Landmarks(icase, iLM), char(style(mod(iLM, 5)+1)), 'LineWidth', 3);
    end
end
ylabel( 'Physical position of the landmark', 'FontSize', 20, 'FontWeight', 'Bold' );
xlabel( 'Mean position of the landmark', 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Warping function', 'FontSize', 20, 'FontWeight', 'Bold' );
set( gca, 'XTick', 0:0.1:1 );
set( gca, 'YTick', 0:0.1:1 );
caxis( [0 nBins] );
h = colorbar;
set( get( h, 'ylabel' ), 'String', 'Age: month(s)', 'FontSize', 20, 'FontWeight', 'Bold' );
set( h, 'ylim', [0, nBins] );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
filename = sprintf( '%s/warpFunctionLandmarks.fig', params.outputPrefix );
saveas( gca, filename );
filename = sprintf( '%s/warpFunctionLandmarks.png', params.outputPrefix );
saveas( gca, filename );


% plot all the registered curves, colored by age
figure('Colormap', cmap); hold on
for icase = params.pos_SGS:size( areamatLM, 2 )
	plot( agefine, areamatLM( :, icase ), 'Color', cmap(round(params.ages(icase)), :), 'LineWidth', 2, 'LineStyle', '-.' );
	hold on
end
for icase = 1:numf
    if icase < params.pos_SGS
        plot( agefine, areamatLM( :, icase ), 'Color', cmap(round(params.ages(icase)), :), 'LineWidth', 1 );
    else
        plot( agefine, areamatLM( :, icase ), 'Color', cmap(round(params.ages(icase)), :), 'LineWidth', 2, 'LineStyle', '-.' );
    end
end
for icase = 1:numf
    warpInvfd = smooth_basis(areawarpmatLM(:,icase), agefine, fdPar_new);
    warpedLM(icase,:) = eval_fd(Landmarks(icase,:), warpInvfd);
    for iLM = 1:size(warpedLM, 2)
        warpedLM(icase, iLM) = max( warpedLM(icase, iLM), 0 );
        warpedLM(icase, iLM) = min( warpedLM(icase, iLM), 1 );
    end
    valueLM(icase, :) = eval_fd(warpedLM(icase,:), areafdLM(icase));
    for iLM = 1:size(warpedLM, 2)
    	plot(warpedLM(icase, iLM), valueLM(icase, iLM), char(style(mod(iLM, 5)+1)), 'LineWidth', 3);
    end
end
hold off
caxis( [0 nBins] );
h = colorbar;
set( get( h, 'ylabel' ), 'String', 'Age: month(s)', 'FontSize', 20, 'FontWeight', 'Bold' );
set( h, 'ylim', [0 nBins] );
xlabel( 'Depth along the centerline', 'FontSize', 20, 'FontWeight', 'Bold' );
ylabel( 'Cross-sectional area (mm^2)', 'FontSize', 20, 'FontWeight', 'Bold' );
if params.pos_SGS <= size( areamatLM, 2 )
	legend( params.cases(params.pos_SGS:end) );
end
set( gca, 'XTick', 0:0.1:1 );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Registered Curves', 'FontSize', 20, 'FontWeight', 'Bold' );
filename = sprintf( '%s/landmark_registration_age.fig', params.outputPrefix );
saveas( gca, filename );
filename = sprintf( '%s/landmark_registration_age.png', params.outputPrefix );
saveas( gca, filename );

% plot all the registered curves, colored by weight
minWeight = min( params.weights );
maxWeight = max( params.weights );
figure('Colormap', cmap), hold on
for icase = 1:numf
    if icase < params.pos_SGS
    	plot( agefine, areamatLM( :, icase ), 'Color', cmap(colorPosW(icase), :), 'LineWidth', 1 );
    else
        plot( agefine, areamatLM( :, icase ), 'Color', cmap(colorPosW(icase), :), 'LineWidth', 2, 'LineStyle', '--' );
    end
end
% plot landmarks
for icase = 1:numf
    for iLM = 1:size(warpedLM, 2)
        plot(warpedLM(icase, iLM), valueLM(icase, iLM), char(style(mod(iLM, 5)+1)), 'LineWidth', 3);
    end
end
hold off
caxis( [minWeight maxWeight] );
h = colorbar
set( h, 'ylim', [minWeight, maxWeight] );
xlabel( 'Depth along the centerline', 'FontSize', 20, 'FontWeight', 'Bold' );
ylabel( 'Cross-sectional area', 'FontSize', 20, 'FontWeight', 'Bold' );
set( gca, 'XTick', 0:0.1:1 );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Registered Curves (purple-red : light-heavy, kg)', 'FontSize', 20, 'FontWeight', 'Bold' );
filename = sprintf( '%s/landmark_registration_weight.png', params.outputPrefix );
set( gcf, 'PaperUnits', 'inches' );
set( gcf, 'PaperSize', [11 8.5] );
set( gcf, 'PaperPositionMode', 'manual' );
set( gcf, 'PaperPosition', [0 0 11 8.5] );
set( gcf, 'renderer', 'painters' );
print( gcf, '-dpng', filename );

filename = sprintf( '%s/landmark_registration_weight.fig', params.outputPrefix );
saveas( gca, filename );
filename = sprintf( '%s/landmark_registration_weight.png', params.outputPrefix );
saveas( gca, filename );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistical analysis: functional boxplot, pca, point-wise boxplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default settings
factor = 1.5;
fullout = false;
barcol = 'b';
outliercol = 'r';
color = 'm';
prob = 0.5;
show = true;
method = 'MBD';
depth = [];

% plot the functional boxplot for all the curves
figure;
for icase = params.pos_SGS : size( areamatLM, 2 )
	plot( agefine, areamatLM(:, icase), 'LineWidth', 2, 'Color', cmap(round(params.ages(icase)), :), 'LineStyle', '-.' );
 	hold on
end
if params.pos_SGS <= size( areamatLM, 2 )
	legend( params.cases(params.pos_SGS:end), 'FontSize', 20, 'FontWeight', 'Bold' );
end
fbplot( areamatLM(:, 1:params.pos_SGS-1), agefine, depth, method, show, prob, color, outliercol, barcol, fullout, factor );
hold on
for icase = params.pos_SGS : size( areamatLM, 2 )
	plot( agefine, areamatLM(:, icase), 'LineWidth', 2, 'Color', cmap(round(params.ages(icase)), :), 'LineStyle', '-.' );
end
hold off
filename = sprintf( '%s/boxplot/functional_boxplot_allCurves_%s.png', params.outputPrefix, method );
set( gca, 'ylim', [0, 1000] );
title( 'Functional boxplot for all curves' );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
set( gcf, 'PaperUnits', 'inches' );
set( gcf, 'PaperSize', [11 8.5] );
set( gcf, 'PaperPositionMode', 'manual' );
set( gcf, 'PaperPosition', [0 0 11 8.5] );
set( gcf, 'renderer', 'painters' );
print( gcf, '-dpng', filename );

% plot the point-wise boxplot for all the curves
figure; 
for icase = params.pos_SGS : size( areamatLM, 2 )
	plot( 1:size(areamatLM, 1), areamatLM(:, icase), 'LineWidth', 2, 'Color', cmap(round(params.ages(icase)), :), 'LineStyle', '-.' );
	hold on
end
if params.pos_SGS <= size( areamatLM, 2 )
	legend( params.cases(params.pos_SGS:end), 'FontSize', 20, 'FontWeight', 'Bold');
end
boxplot( (areamatLM(:, 1:params.pos_SGS-1))' );
hold on
for icase = params.pos_SGS : size( areamatLM, 2 )
	plot( 1:size(areamatLM, 1), areamatLM(:, icase), 'LineWidth', 2, 'Color', cmap(round(params.ages(icase)), :), 'LineStyle', '-.' );
end
hold off
set( gca, 'ylim', [0, 1000] );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Pointwise boxplot for all curves' );
filename = sprintf( '%s/boxplot/pointwise_boxplot_allCurves_allPoints.png', params.outputPrefix );
set( gcf, 'PaperUnits', 'inches' );
set( gcf, 'PaperSize', [11 8.5] );
set( gcf, 'PaperPositionMode', 'manual' );
set( gcf, 'PaperPosition', [0 0 11 8.5] );
set( gcf, 'renderer', 'painters' );
print( gcf, '-dpng', filename );


% PCA and boxplot performed on parts of the curves

% find the id for each landmark after registration
landmarksNewId = zeros(size(warpedLM,2), 1);
for iSample = 1:size(warpedLM,2)
    [nValue, nId] = min( abs( agefine - mean( warpedLM( :, iSample ) ) ) );
    landmarksNewId( iSample ) = nId;
end

% subsection of curves: from the current landmark to the end
part_sample = [ size( warpedLM, 2 ) - 1 ];    % from TVC to trachea carina
for iSample = 1:length( part_sample )
	nStartLM = part_sample(iSample);
	[nStartValue, nStartId] = min( abs( agefine - mean( warpedLM( :, nStartLM ) ) ) );
	[nEndValue, nEndId] = min( abs( agefine - mean( warpedLM( :, nStartLM+1 ) ) ) );
	area_part = areamatLM( nStartId:end, : );
	agefine_part = agefine( nStartId:end ); 
	yLim = ceil( max( max(area_part) ) / 100.0 ) * 100;
	rng = [ min(agefine_part), max(agefine_part) ];
	knots = agefine_part;
	norder = 6;
	nbasis = length(knots) + norder - 2;
	areabasis_part = create_bspline_basis( rng, nbasis, norder, knots );
	Lfdobj = int2Lfd(2);
	lambda = 0.5 * 1e-6;                      % adjust for fitting
	if part_sample(iSample) == 1
		lambda = 0.5 * 1e-5;
	end
	areafdPar_part = fdPar( areabasis_part, Lfdobj, lambda );
	areafd_part = smooth_basis( agefine_part, area_part, areafdPar_part );
	areamatLM_part = eval_fd( agefine_part, areafd_part );

	% plot the subsection of the curves
	figure, hold on
	for icase = params.pos_SGS : size( areamatLM_part, 2 )
		plot( agefine_part, areamatLM_part(:, icase), 'Color', cmap( round(params.ages(icase)), : ), 'LineWidth', 2, 'LineStyle', '-.' );
	end
	for icase = 1:numf
	    if icase < params.pos_SGS
	    	plot( agefine_part, areamatLM_part(:, icase), 'Color', cmap( round(params.ages(icase)), : ), 'LineWidth', 1 );
	    else
		plot( agefine_part, areamatLM_part(:, icase), 'Color', cmap( round(params.ages(icase)), : ), 'LineWidth', 1, 'LineStyle', '-.' );
	    end
	end
	hold off
	set( gca, 'xlim', [min(agefine_part), max(agefine_part)] );
	set( gca, 'ylim', [0, yLim] );
	if params.pos_SGS <= size( areamatLM_part, 2 )
		legend( params.cases{params.pos_SGS:end} );
	end
	set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
	filename = sprintf( '%s/PCA/Curve_Part_StartFromLandmark%d.png', params.outputPrefix, part_sample(iSample) );
	set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	set( gcf, 'PaperUnits', 'inches' );
	set( gcf, 'PaperSize', [11 8.5] );
	set( gcf, 'PaperPositionMode', 'manual' );
	set( gcf, 'PaperPosition', [0 0 11 8.5] );
	set( gcf, 'renderer', 'painters' );
	print( gcf, '-dpng', filename );

	% principle component analysis
	nharm = 5;
	areafd_part_CRL = areafd_part(1:params.pos_SGS-1);
	areafdLM_PCA_part = pca_fd( areafd_part_CRL-mean(areafd_part_CRL), nharm );
	areamatLM_PCA_part = eval_fd( agefine_part, areafdLM_PCA_part.harmfd );
	for iharm = 1:nharm
	    newmat1_part = areamatLM_PCA_part( :, iharm ) * sqrt( areafdLM_PCA_part.values( iharm ) ) + mean( areamatLM_part(:, 1:params.pos_SGS-1), 2 );
	    newmat2_part = areamatLM_PCA_part( :, iharm ) * -sqrt( areafdLM_PCA_part.values( iharm ) ) + mean( areamatLM_part(:, 1:params.pos_SGS-1), 2 );
	    figure, plot( agefine_part, mean( areamatLM_part(:, 1:params.pos_SGS-1), 2 ), 'r-', 'LineWidth', 2 );
	    hold on
	    plot( agefine_part, newmat1_part, 'g--', 'LineWidth', 2 );
	    plot( agefine_part, newmat2_part, 'b-.', 'LineWidth', 2 );
	    hold off
	    title(['Cross-sectional Area: Part PC', int2str(iharm), ' (', num2str(areafdLM_PCA_part.varprop(iharm)*100, '%10.2f'), '%)'], ...
		'FontSize', 20, 'FontWeight', 'Bold');
	    filename = sprintf( '%s/PCA/PC%02d_part.png', params.outputPrefix, iharm );
	    set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
	    set( gca, 'xlim', [min(agefine_part), max(agefine_part)] );
	    set( gca, 'ylim', [0, yLim] );
	    legend( 'Mean', 'Mean+\lambda*v', 'Mean-\lambda*v' );
	    set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	    set( gcf, 'PaperUnits', 'inches' );
	    set( gcf, 'PaperSize', [11 8.5] );
	    set( gcf, 'PaperPositionMode', 'manual' );
	    set( gcf, 'PaperPosition', [0 0 11 8.5] );
	    set( gcf, 'renderer', 'painters' );
	    print( gcf, '-dpng', filename );

	    % plot all the curves with just one component
	    figure,
	    for icase = 1:size(areafdLM_PCA_part.harmscr, 1)
            curveComponent_part = areamatLM_PCA_part( :, iharm ) * areafdLM_PCA_part.harmscr( icase, iharm ) + mean( areamatLM_part(:, 1:params.pos_SGS-1), 2 );
            plot( agefine_part, curveComponent_part, 'Color', cmap(round(params.ages( icase )), :), 'LineWidth', 1 );
            hold on
	    end
	    hold off
	    title( ['Curve component (Part) ', int2str(iharm)], 'FontSize', 20, 'FontWeight', 'Bold' );
	    caxis( [0 nBins] );
	    h = colorbar;
	    set( h, 'ylim', [0, nBins] );
	    set( gca, 'xlim', [min(agefine_part), max(agefine_part)] );
	    set( gca, 'ylim', [0, yLim] );
	    set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
	    filename = sprintf( '%s/PCA/CurveComponentPart%02d_StartFromLandmark%d.png', params.outputPrefix, iharm, part_sample(iSample) );
	    set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	    set( gcf, 'PaperUnits', 'inches' );
	    set( gcf, 'PaperSize', [11 8.5] );
	    set( gcf, 'PaperPositionMode', 'manual' );
	    set( gcf, 'PaperPosition', [0 0 11 8.5] );
	    set( gcf, 'renderer', 'painters' );
	    print( gcf, '-dpng', filename );

	    % display scores
	    figure, hold on
	    for icase = 1:size(areafdLM_PCA_part.harmscr, 1 )
            plot( areafdLM_PCA_part.harmscr( icase, iharm ), rand(1)/100, 'Color', cmap( round(params.ages( icase )), : ), 'Marker', 'o', 'MarkerSize', 6, 'LineWidth', 2 );
	    end
	    hold off
	    set( gca, 'ylim', [-0.05, 0.05] );
	    set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
	    title( 'Scores of part of the curves' );
	    filename = sprintf( '%s/PCA/Score_part%02d_StartFromLandmark%d.png', params.outputPrefix, iharm, part_sample(iSample) );
	    set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	    set( gcf, 'PaperUnits', 'inches' );
	    set( gcf, 'PaperSize', [11 8.5] );
	    set( gcf, 'PaperPositionMode', 'manual' );
	    set( gcf, 'PaperPosition', [0 0 11 8.5] );
	    set( gcf, 'renderer', 'painters' );
	    print( gcf, '-dpng', filename );
	end
	figure, hold on
	for icase = 1:size(areafdLM_PCA_part.harmscr, 1 )
	    plot( areafdLM_PCA_part.harmscr( icase, 1 ), areafdLM_PCA_part.harmscr( icase, 2 ), 'Color', cmap( round(params.ages( icase )), : ), 'Marker', 'o', 'MarkerSize', 6, 'LineWidth', 2 );
	end
	hold off
	filename = sprintf( '%s/PCA/Score1_Score2_part_StartFromLandmark%d.png', params.outputPrefix, part_sample(iSample) );
 	set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
	set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	set( gcf, 'PaperUnits', 'inches' );
	set( gcf, 'PaperSize', [11 8.5] );
	set( gcf, 'PaperPositionMode', 'manual' );
	set( gcf, 'PaperPosition', [0 0 11 8.5] );
	set( gcf, 'renderer', 'painters' );
	print( gcf, '-dpng', filename );

    % display PCA results in one figure
	displayPCA(agefine_part, areafdLM_PCA_part, areamatLM_PCA_part, areamatLM_part, areafd_part, cmap, params.ages, params.outputPrefix, part_sample(iSample), params.pos_SGS);

    % age-adapted atlas 
	outputPrefix_part = sprintf( '%s/Age/boxplotPart%d', params.outputPrefix, part_sample(iSample) );
    displayBoxplotForSubjects( agefine_part, areamatLM(landmarksNewId(part_sample(iSample)):end, :), params.ages, params.cases, outputPrefix_part, yLim, params.gaussian_sigma, params.uniform_winSize, cmap, params.pos_SGS, params.pos_SGS_Post, areamatLM, landmarksNewId, nStartLM );
	
    % weight-adapted atlas
    %outputPrefix_part = sprintf( '%s/Weight/boxplotPart%d', outputPrefix, part_sample(iSample) );
	%displayBoxplotForSubjects( agefine_part, areamatLM(landmarksNewId(part_sample(iSample)):end, :), params.weights, params.cases, outputPrefix_part, yLim, params.gaussian_sigma/3, params.uniform_winSize/3, cmap, params.pos_SGS, areamatLM, landmarksNewId, nStartLM );
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction: plot results of pca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayPCA(agefine_part, areafdLM_PCA_part, areamatLM_PCA_part, areamatLM_part, areafd_part, cmap, ages, outputPrefix, partId, pos_SGS)
	close all
	if pos_SGS <= size(areamatLM_part, 2)
		areaCoef = getcoef( (areafd_part - mean(areafd_part(1:pos_SGS-1))) );
		sgsScores = (areaCoef(:,pos_SGS:end))' * (areafdLM_PCA_part.Lmat)' * (areafdLM_PCA_part.eigvecs);
	end
	fprintf( 'size part ....' );
	size( areafdLM_PCA_part.harmscr, 1 )
	figure;
    for iSub = 1:4
    	ix = floor( (iSub-1)/2 ) + 1;
        iy = mod( iSub-1, 2 ) + 1;
        subplot( 2, 2, iSub ),
        if ix == iy
        	for icase = 1:size(areafdLM_PCA_part.harmscr, 1)
		        curveComponent_part = areamatLM_PCA_part( :, ix ) * areafdLM_PCA_part.harmscr( icase, ix ) + mean( areamatLM_part(:, 1:pos_SGS-1), 2 );
				plot( agefine_part, curveComponent_part, 'Color', cmap(round(ages( icase )), :), 'LineWidth', 1 );
                hold on
            end
			for icase = pos_SGS:size(areamatLM_part, 2)
				curveComponent_part = areamatLM_PCA_part( :, ix ) * sgsScores(icase-pos_SGS+1, ix) + mean( areamatLM_part( :, 1:pos_SGS-1 ), 2 );
				plot( agefine_part, curveComponent_part, 'Color', cmap( round(ages( icase )), : ), 'LineWidth', 2, 'LineStyle', '-.' );
			end
            hold off
			axis tight
			xlim( [min(agefine_part) max(agefine_part)] );
			ylim( [min(min(areamatLM_part)) max(max(areamatLM_part))] );
			set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
			xlabel( 'Depth along the centerline', 'FontSize', 20, 'FontWeight', 'Bold' );
			ylabel( 'X-sectional area (mm^2)', 'FontSize', 20, 'FontWeight', 'Bold'  );
			title_name = sprintf( 'PCA Component %d, %3.2f%%', ix, areafdLM_PCA_part.varprop(ix)*100 );
			title( title_name, 'FontSize', 20, 'FontWeight', 'Bold' );
        else
            for icase = 1:size(areafdLM_PCA_part.harmscr,1)
				plot( areafdLM_PCA_part.harmscr( icase, ix ), areafdLM_PCA_part.harmscr( icase, iy ), 'Color', cmap( round(ages( icase )), : ), 'Marker', 'o', 'MarkerSize', 8, 'LineWidth', 2 );
				hold on
            end
			for icase = pos_SGS:size(areamatLM_part, 2)
				plot( sgsScores(icase-pos_SGS+1, ix), sgsScores(icase-pos_SGS+1, iy), 'Color', cmap( round(ages(icase)), : ), 'Marker', '*', 'MarkerSize', 9, 'LineWidth', 2 ); 
			end
            hold off
			axis tight
			xlabel_name = sprintf( 'PCA %d', ix );
			ylabel_name = sprintf( 'PCA %d', iy );
			xlabel( xlabel_name, 'FontSize', 20, 'FontWeight', 'Bold'  );
			ylabel( ylabel_name, 'FontSize', 20, 'FontWeight', 'Bold'  );
			set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
			title( 'Scores', 'FontSize', 20, 'FontWeight', 'Bold'  );
        end
        hold off
    end
    
    filename = sprintf( '%s/PCA/Part_StartFromLandmark%d.fig', outputPrefix, partId );
    saveas( gca, filename );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction: plot boxplot for both control and SGS subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function displayBoxplotForSubjects( agefine, areamatLM, ages, cases, outputPrefix, yLim, gaussian_sigma, uniform_winSize, cmap, pos_SGS, pos_SGS_Post, areamatLM_whole, landmarksId, iSample )
close all
bflag = true;
% default settings
factor = 1.5;
fullout = false;
barcol = 'b';
outliercol = 'r';
color = 'm';
prob = 0.5;
show = true;
method = 'MBD';
depth = [];

nBins = ceil(max(ages)/10.0) * 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% functional boxplot for all curves %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fbplot of all the curves
figure;
for icase = pos_SGS : size( areamatLM, 2 )
	plot( agefine, areamatLM(:, icase), 'Color', cmap(round(ages(icase)), :), 'LineWidth', 2, 'LineStyle', '-.' );
	hold on
end
fbplot( areamatLM(:, 1:pos_SGS-1), agefine, depth, method, show, prob, color, outliercol, barcol, fullout, factor );
hold on
for icase = pos_SGS : size( areamatLM, 2 )
	plot( agefine, areamatLM(:, icase), 'Color', cmap(round(ages(icase)), :), 'LineWidth', 2, 'LineStyle', '-.' );
end
if pos_SGS <= size( areamatLM, 2 )
	legend( cases( pos_SGS:end ), 'FontSize', 20, 'FontWeight', 'Bold' );
end
hold off
filename = sprintf( '%s/functional_boxplot_allCurvesWithSGS_%s.png', outputPrefix, method );
set( gca, 'xlim', [min(agefine), max(agefine)] );
set( gca, 'ylim', [0, yLim] );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Functional boxplot for all curves' );
if bflag 
	set( gcf, 'PaperUnits', 'inches' );
	set( gcf, 'PaperSize', [11 8.5] );
	set( gcf, 'PaperPositionMode', 'manual' );
	set( gcf, 'PaperPosition', [0 0 11 8.5] );
	set( gcf, 'renderer', 'painters' );
	print( gcf, '-dpng', filename );
else
	saveas( gca, filename );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% point-wise boxplot for all curves %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for icase = pos_SGS : size( areamatLM, 2 )
	plot( 1:size(areamatLM, 1), areamatLM(:, icase), 'Color', cmap(round(ages(icase)), :), 'LineWidth', 2, 'LineStyle', '-.' );
	hold on
end
boxplot( areamatLM(:, 1:pos_SGS-1)' );
hold on
for icase = pos_SGS : size( areamatLM, 2 )
	plot( 1:size(areamatLM, 1), areamatLM(:, icase), 'Color', cmap(round(ages(icase)), :), 'LineWidth', 2, 'LineStyle', '-.' );
end
if pos_SGS <= size( areamatLM, 2 )
	legend( cases( pos_SGS:end ), 'FontSize', 20, 'FontWeight', 'Bold' );
end
hold off
set( gca, 'ylim', [0, yLim] );
set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
title( 'Pointwise boxplot for all curves' );
filename = sprintf( '%s/pointwise_boxplot_allCurves_allPointsWithSGS.png', outputPrefix );
if bflag 
	%set( gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1] );
	set( gcf, 'PaperUnits', 'inches' );

	set( gcf, 'PaperSize', [11 8.5] );
	set( gcf, 'PaperPositionMode', 'manual' );
	set( gcf, 'PaperPosition', [0 0 11 8.5] );
	set( gcf, 'renderer', 'painters' );
	print( gcf, '-dpng', filename );
else
	saveas( gca, filename );
end

% score and age of a curve for point-wise, functional, and weighted functional boxplots
posInAtlas = zeros( size( areamatLM, 2 ), 6 );

% window sliding
start_id_tmp = 1;  % build atlas for all the cases, including both SGS and controls
for icase = start_id_tmp : size( areamatLM, 2 )
  ageId = ages( icase );
  pickedCurves = (ages(1:pos_SGS-1) >= ageId - uniform_winSize) .* (ages(1:pos_SGS-1) <= ageId + uniform_winSize);
 
  close all;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  figure 1, square window, functional boxplot, left %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % square window, functional boxplot
  figure('Position', [100, 100, 1024, 786] ), subplot(2,2,1), 
  plot( agefine, areamatLM(:, icase ), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap(round(ageId), :) );
  hold on
  [depthTmp, outliers, medianCurveId, posPercentile] = fbplot( areamatLM( :, logical( pickedCurves ) ), ...
                  agefine, depth, method, show, prob, color, outliercol, barcol, fullout, factor, areamatLM( :, icase ) );
  hold on
  plot( agefine, areamatLM(:, icase ), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap(round(ageId), :) );
  hold off
  legend( cases{icase});
  iCountId = 0;
  originalMedianCurveId = -1;
  for iCurveId = 1:length(pickedCurves)
      if pickedCurves( iCurveId ) > 0
          iCountId = iCountId + 1;
      	  if iCountId == medianCurveId
	      originalMedianCurveId = iCurveId;
              break;
          end
      end
  end
  title( 'Functional boxplot', 'FontSize', 20, 'FontWeight', 'Bold' );
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'xlim', [min(agefine), max(agefine)] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
   
  subplot(2,2,3), 
  hist(ages(1:pos_SGS-1), nBins);
  hold on
  %bar( ages(originalMedianCurveId), length(find(ages==ages(originalMedianCurveId))), 'm' );
  plot( [ages(originalMedianCurveId), ages(originalMedianCurveId)], [0, length(find(ages(1:pos_SGS-1)==ages(originalMedianCurveId)))], 'm', 'LineWidth', 2 );
  textString = sprintf( '\\downarrow %d (median)', ages(originalMedianCurveId) );
  text( ages(originalMedianCurveId)-2, 1.35, textString, 'Color', 'm', 'FontSize', 20, 'FontWeight', 'Bold');
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-2, -1.0, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  rectangle('Position', [max(ageId-uniform_winSize,0.1), 0, min(ageId+uniform_winSize, nBins)-max(ageId-uniform_winSize, 0), 1.1], 'Curvature', [0, 0], 'LineWidth', 2, 'LineStyle', '-.', 'EdgeColor', 'r');
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
  hold off


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  figure 1, square window, point-wise boxplot, right %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % square window, point-wise boxplot
  %figure, boxplot( areamatLM( tmpLinePoint( logical( mod( tmpLinePoint, 30 ) == 0 ) ), logical( (ages >= ageId - sizeWindow) .* (ages <= ageId + sizeWindow) ) )' );
  subplot(2,2,2),
  plot( 1:size(areamatLM, 1), areamatLM( :, icase ), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap( round(ageId), : ) ); 
  hold on
  boxplot( areamatLM( 1:size(areamatLM, 1), logical( pickedCurves ) )', 'boxstyle', 'filled' );
  hold on
  plot( 1:size(areamatLM, 1), areamatLM( :, icase ), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap( round(ageId), : ) ); 
  hold off
  legend( cases{icase} );
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
  title( 'Pointwise boxplot', 'FontSize', 20, 'FontWeight', 'Bold' );
  
  subplot(2,2,4), 
  hist(ages(1:pos_SGS-1), nBins);
  hold on
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-2, -1.0, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  rectangle('Position', [max(ageId-uniform_winSize,0.1), 0, min(ageId+uniform_winSize, nBins)-max(ageId-uniform_winSize, 0), 1.1], 'Curvature', [0, 0], 'LineWidth', 2, 'LineStyle', '-.', 'EdgeColor', 'r');
  hold off
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
  filename = sprintf( '%s/functional_pointwise_boxplot_%s_age%0.2f.fig', outputPrefix, cases{icase}, ageId );
  saveas( gca, filename );
  filename = sprintf( '%s/functional_pointwise_boxplot_%s_age%0.2f.png', outputPrefix, cases{icase}, ageId );
  saveas( gca, filename );

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  figure 2, gaussian window, weighted functional boxplot, left %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % weighted gaussian window
  %w_j_unNormalized = gaussmf( ages(1:pos_SGS-1), [ gaussian_sigma, ageId ] );
  w_j_unNormalized = gaussianWeighting( ages(1:pos_SGS-1), 0, 192, gaussian_sigma, ageId );
  w_j = w_j_unNormalized / sum( w_j_unNormalized ) ;
  figure('Position', [100, 100, 1024, 786]), subplot(2,2,1), 
  plot( agefine, areamatLM( :, icase ), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap( round(ageId), : ) );
  hold on
  [depthTmp, outliers, medianCurveId, posPercentile, centerId] = wfbplot( areamatLM(:, 1:pos_SGS-1), agefine, w_j, depth, method, show, prob, color, outliercol, barcol, fullout, factor, areamatLM( :, icase ) );
  textString = sprintf( 'Percentile: %f, Age: %d', posPercentile, round(ageId) );
  text( mean(agefine), yLim*0.8, textString );
  hold on
  plot( agefine, areamatLM( :, icase ), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap( round(ageId), : ) );
  hold off
  legend( cases{icase} );
  titlename = sprintf( 'Weighted functional boxplot' );
  title( titlename, 'FontSize', 20, 'FontWeight', 'Bold' );
  set( gca, 'xlim', [min(agefine), max(agefine)] );
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );

  subplot(2,2,3), 
  hist(ages(1:pos_SGS-1), nBins);
  hold on
  %bar( ages(medianCurveId), length(find(ages==ages(medianCurveId))), 'm' );
  plot( [ages(medianCurveId), ages(medianCurveId)], [0, length(find(ages(1:pos_SGS-1)==ages(medianCurveId)))], 'm', 'LineWidth', 2);
  textString = sprintf( '\\downarrow %d (median)', ages(medianCurveId) );
  text( ages(medianCurveId)-2, 1.35, textString, 'Color', 'm', 'FontSize', 20, 'FontWeight', 'Bold');
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-2, -1.0, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  [sort_age, sort_Id] = sort(ages(1:pos_SGS-1));
  plot(sort_age, w_j_unNormalized(sort_Id), 'r-.*', 'LineWidth', 2);
  hold off
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  figure 2, gaussian window, wighted point-wise boxplot, right %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % the quartile_5_curve and quartile_95_curve are actually the minimal and 
  % maximal curves, should change the name later
  [median_curve, quartile_25_curve, quartile_75_curve, quartile_5_curve, quartile_95_curve, scorePW] = ...
      wpbplot( areamatLM(:, 1:pos_SGS-1), w_j, factor, areamatLM( :, icase ));

  % store the score based on the weighted point-wise boxplot
  posInAtlas( icase, 1 ) = scorePW;
  posInAtlas( icase, 2 ) = ageId;

  subplot(2,2,2), 
  plot(agefine, quartile_95_curve, 'r-', 'LineWidth', 2 );
  hold on
  plot(agefine, quartile_75_curve, 'g-.', 'LineWidth', 2 );
  plot(agefine, median_curve, 'm-', 'LineWidth', 2);
  plot(agefine, quartile_25_curve, 'c-.', 'LineWidth', 2 );
  plot(agefine, quartile_5_curve, 'b-', 'LineWidth', 2);
  legend( 'Maximum', '75%', 'Median', '25%', 'Minimum' );
  plot( agefine, areamatLM(:, icase), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap(round(ageId), :) );
  hold off
  set( gca, 'xlim', [min(agefine), max(agefine)] );
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
  title( 'Weighted percentile', 'FontSize', 20, 'FontWeight', 'Bold' );

  subplot(2,2,4), 
  hist(ages(1:pos_SGS-1), nBins);
  hold on
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-0.7, -1.0, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  plot(sort_age, w_j_unNormalized(sort_Id), 'r-.*', 'LineWidth', 2);
  hold off
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
  filename = sprintf( '%s/weighted_functional_pointwise_boxplot_%s_age%0.2f.fig', outputPrefix, cases{icase}, ageId );
  saveas( gca, filename );
  filename = sprintf( '%s/weighted_functional_pointwise_boxplot_%s_age%0.2f.png', outputPrefix, cases{icase}, ageId );
  saveas( gca, filename );

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  figure 3, square window, functional boxplot, left %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  h=figure('Position', [100, 100, 1024, 786] );
  subplot(2,2,1), plot( agefine, areamatLM(:, icase), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap(round(ageId), :) );
  hold on
  [depthTmp, outliers, medianCurveId, posPercentile] = fbplot( areamatLM( :, logical( pickedCurves ) ), ...
                  agefine, depth, method, show, prob, color, outliercol, barcol, fullout, factor, areamatLM( :, icase ) );
  hold on
  plot( agefine, areamatLM(:, icase), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap(round(ageId), :) );
  hold off
  iCountId = 0;
  originalMedianCurveId = -1;
  for iCurveId = 1:length(pickedCurves)
      if pickedCurves( iCurveId ) > 0
          iCountId = iCountId + 1;
          if iCountId == medianCurveId
              originalMedianCurveId = iCurveId;
              break;
          end
      end
  end
  legend( cases{icase} );
  title( 'Functional boxplot', 'FontSize', 20, 'FontWeight', 'Bold' );
  set( gca, 'xlim', [min(agefine), max(agefine)] );
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );

  subplot(2,2,3), 
  hist(ages(1:pos_SGS-1), nBins);
  hold on
  plot( [ages(originalMedianCurveId), ages(originalMedianCurveId)], [0, length(find(ages(1:pos_SGS-1)==ages(originalMedianCurveId)))], 'm', 'LineWidth', 2 );
  textString = sprintf( '\\downarrow %d (median)', ages(originalMedianCurveId) );
  text( ages(originalMedianCurveId)-2, 1.35, textString, 'Color', 'm', 'FontSize', 20, 'FontWeight', 'Bold');
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-2, -1.0, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  rectangle('Position', [max(ageId-uniform_winSize,0.1), 0, min(ageId+uniform_winSize, nBins)-max(ageId-uniform_winSize, 0), 1.1], 'Curvature', [0, 0], 'LineWidth', 2, 'LineStyle', '-.', 'EdgeColor', 'r');
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );
  hold off

  % store the score based on the functional boxplot
  posInAtlas( icase, 3 ) = posPercentile;
  posInAtlas( icase, 4 ) = ageId;
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%  figure 3, gaussian window, weighted functional boxplot, right %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,2), 
  plot( agefine, areamatLM(:, icase), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap(round(ageId), :) );
  hold on
  [depthTmp, outliers, medianCurveId, posPercentile] = wfbplot( areamatLM(:, 1:pos_SGS-1), agefine, w_j, depth, method, show, prob, color, outliercol, barcol, fullout, factor, areamatLM(:, icase) );
  textString = sprintf( 'Percentile: %f, Age: %d', posPercentile, round(ageId) );
  text( mean(agefine), mean(areamatLM(:, icase)), textString );
  hold on
  plot( agefine, areamatLM(:, icase), 'LineWidth', 2, 'LineStyle', '-.', 'Color', cmap(round(ageId), :) );
  hold off
  legend( cases{icase} );
  titlename = sprintf( 'Weighted functional boxplot' );
  title( titlename, 'FontSize', 20, 'FontWeight', 'Bold' );
  set( gca, 'xlim', [min(agefine), max(agefine)] );
  set( gca, 'ylim', [0, yLim] );
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );

  subplot(2,2,4), 
  %hist(ages, max(ages));
  hist(ages(1:pos_SGS-1), nBins);
  hold on
  %bar( ages(medianCurveId), length(find(ages==ages(medianCurveId))), 'm' );
  plot( [ages(medianCurveId), ages(medianCurveId)], [0, length(find(ages(1:pos_SGS-1)==ages(medianCurveId)))], 'm', 'LineWidth', 2);
  textString = sprintf( '\\downarrow %d (median)', ages(medianCurveId) );
  text( ages(medianCurveId)-2, 1.35, textString, 'Color', 'm', 'FontSize', 20, 'FontWeight', 'Bold');
  plot([ageId, ageId], [0, 1], 'r-', 'LineWidth', 2);
  textString = sprintf( '\\uparrow %d months', ageId );
  text( ageId-2, -1.0, textString, 'Color', 'r', 'FontSize', 20, 'FontWeight', 'Bold');
  plot(sort_age, w_j_unNormalized(sort_Id), 'r-.*', 'LineWidth', 2);
  hold off 
  set( gca, 'FontSize', 20, 'FontWeight', 'Bold' );

  % store the score based on the weighted functional boxplot
  posInAtlas( icase, 5 ) = posPercentile;
  posInAtlas( icase, 6 ) = ageId;

  filename = sprintf( '%s/weighted_unweighted_functional_boxplot_%s_age%0.2f.fig', outputPrefix, cases{icase}, ageId );
  saveas( gca, filename );
  filename = sprintf( '%s/weighted_unweighted_functional_boxplot_%s_age%0.2f.png', outputPrefix, cases{icase}, ageId );
  saveas( gca, filename );
end

% plot the scores for subjects (weighted pointwise boxplots)
index_pw_tmp = 1;
figure, hold on
for icase = 1:size(posInAtlas, 1)
	if icase < pos_SGS
		plot( posInAtlas( icase, index_pw_tmp), posInAtlas( icase, index_pw_tmp+1), 'Marker', 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 8 ); 
	elseif icase < pos_SGS_Post   % the number of the pre-surgery subjects
		plot( posInAtlas( icase, index_pw_tmp), posInAtlas( icase, index_pw_tmp+1), 'Marker', 'p', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 12 );
	else
		plot( posInAtlas( icase, index_pw_tmp), posInAtlas( icase, index_pw_tmp+1), 'Marker', 's', 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm', 'MarkerSize', 10 );
	end
end
hold off
title('Score based on functional boxplots');
filename = sprintf( '%s/posInAtlas_pw.fig', outputPrefix );
saveas( gca, filename );
filename = sprintf( '%s/posInAtlas_pw.png', outputPrefix );
saveas( gca, filename );

% plot the scores for subjects (functional boxplots)
index_fb_tmp = 3;
figure, hold on
for icase = 1:size(posInAtlas, 1)
	if icase < pos_SGS
		plot( posInAtlas( icase, index_fb_tmp), posInAtlas( icase, index_fb_tmp+1), 'Marker', 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 8 ); 
	elseif icase < pos_SGS_Post   % the number of the pre-surgery subjects
		plot( posInAtlas( icase, index_fb_tmp), posInAtlas( icase, index_fb_tmp+1), 'Marker', 'p', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 12 );
	else
		plot( posInAtlas( icase, index_fb_tmp), posInAtlas( icase, index_fb_tmp+1), 'Marker', 's', 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm', 'MarkerSize', 10 );
	end
end
hold off
title('Score based on functional boxplots');
filename = sprintf( '%s/posInAtlas_fb.fig', outputPrefix );
saveas( gca, filename );
filename = sprintf( '%s/posInAtlas_fb.png', outputPrefix );
saveas( gca, filename );

% plot the scores for subjects (weighted functional boxplots)
index_wfb_tmp = 5;
figure, hold on
for icase = 1:size(posInAtlas, 1)
	if icase < pos_SGS
		plot( posInAtlas( icase, index_wfb_tmp), posInAtlas( icase, index_wfb_tmp+1), 'Marker', 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'MarkerSize', 8 ); 
	elseif icase < pos_SGS_Post   % the number of the pre-surgery subjects
		plot( posInAtlas( icase, index_wfb_tmp), posInAtlas( icase, index_wfb_tmp+1), 'Marker', 'p', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 12 );
	else
		plot( posInAtlas( icase, index_wfb_tmp), posInAtlas( icase, index_wfb_tmp+1), 'Marker', 's', 'MarkerEdgeColor', 'm', 'MarkerFaceColor', 'm', 'MarkerSize', 10 );
	end
end
hold off
title('Score based on weighted functional boxplots');
filename = sprintf( '%s/posInAtlas_wfb.fig', outputPrefix );
saveas( gca, filename );
filename = sprintf( '%s/posInAtlas_wfb.png', outputPrefix );
saveas( gca, filename );

filename = sprintf( '%s/diagnosis_data/posInAtlas.mat', outputPrefix );
save( filename, 'posInAtlas' );

filename = sprintf( '%s/diagnosis_data/area.mat', outputPrefix );
save( filename, 'areamatLM_whole' );

filename = sprintf( '%s/diagnosis_data/cases.mat', outputPrefix );
save( filename, 'cases' );

end

%% generate a colormap
function [cmap, nBins] = generateColormap( ages )
    nBins = ceil( max(ages)/100.0 ) * 100;
    cmap = zeros( nBins, 3 );
    stepColor = 1.0/(nBins/6.0);
    for iI = 1:nBins
        if iI < (nBins/12.0)
            cmap( iI, 1 ) = 1;
            cmap( iI, 2 ) = 0.5 - iI * stepColor;
            cmap( iI, 3 ) = 1;
        elseif iI < (nBins/12.0*3)
            cmap( iI, 1 ) = 1 - (iI-nBins/12.0) * stepColor;
            cmap( iI, 2 ) = 0;
            cmap( iI, 3 ) = 1;
        elseif iI < (nBins/12.0*6)
            cmap( iI, 1 ) = 0;
            cmap( iI, 2 ) = (iI - nBins/12.0*3) * stepColor * 2/3;
            cmap( iI, 3 ) = 1;
        elseif iI < (nBins/12.0*7)
            cmap( iI, 1 ) = 0;
            cmap( iI, 2 ) = 1;
            cmap( iI, 3 ) = 1 - (iI - nBins/12.0*6.0) * stepColor * 2;
        elseif iI < (nBins/12.0*8)
            cmap( iI, 1 ) = (iI - nBins/12.0*7.0) * stepColor * 2;
            cmap( iI, 2 ) = 1;
            cmap( iI, 3 ) = 0;
        elseif iI < (nBins/12.0*11)
            cmap( iI, 1 ) = 1;
            cmap( iI, 2 ) = 1 - (iI - nBins/12.0*8.0) * stepColor * 2 / 3;
            cmap( iI, 3 ) = 0;
        else
            cmap( iI, 1 ) = 1 - (iI - nBins/12.0*11.0) * stepColor;
            cmap( iI, 2 ) = 0;
            cmap( iI, 3 ) = 0;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end ----
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
