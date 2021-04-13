% EXAMPLE6 
% Reads data from file and introduce demonstrates reconstruction with a 
% dictionary consisting of rotations of a basic motif
close all; clear
initpkg

%% ===== Read line data ===== $$
% Assign dataset

% %%%% 1 triangle_1
% d = DataSpec('030220',1,7,1);    
% d.set_scan_distance(6.31e-3, 4.25, 240, 10);
% angles = [0,20,40,60,80,100,120];

% %%%% 1 triangle_2
% d = DataSpec('030220',1,8,2);    
% d.set_scan_distance(6.31e-3, 4.25, 240, 10);
% angles = [0,20,40,60,80,100,120,140];

% %%%% 1 triangle_2 with 3 scans(0,60,120)
% d = DataSpec('030220',1,3,2);    
% d.set_scan_distance(6.31e-3, 4.25, 240, 10);
% angles = [0,60,120];

% %%%% 1 triangle_2 with 3 scans(0,20,40)
% d = DataSpec('030220',1,3,3);    
% d.set_scan_distance(6.31e-3, 4.25, 240, 10);
% angles = [0,20,40];

% %%%% 1 triangle_2 with 3 scans(100,120,140)
% d = DataSpec('030220',1,3,4);    
% d.set_scan_distance(6.31e-3, 4.25, 240, 10);
% angles = [100,120,140];

% %%%% 1 triangle_2 with 4 scans(0,60,120,140)
% d = DataSpec('030220',1,4,2);    
% d.set_scan_distance(6.31e-3, 4.25, 240, 10);
% angles = [0,60,120,140];

% %%%% 1 triangle_2 with 4 scans(0,20,40,60)
% d = DataSpec('030220',1,4,3);    
% d.set_scan_distance(6.31e-3, 4.25, 240, 10);
% angles = [0,20,40,60];

% %%%% 1 triangle_2 with 5 scans
% d = DataSpec('030220',1,5,2);    
% d.set_scan_distance(6.31e-3, 4.25, 240, 10);
% angles = [0,20,60,120,140];

% %%%% 1 triangle_2 with 5 scans(0,20,40,60,80)
% d = DataSpec('030220',1,5,3);    
% d.set_scan_distance(6.31e-3, 4.25, 240, 10);
% angles = [0,20,40,60,80];

% %%%% 1 triangle_2 with 6 scans
% d = DataSpec('030220',1,6,2);    
% d.set_scan_distance(6.31e-3, 4.25, 240, 10);
% angles = [0,20,40,60,120,140];

% %%%% 1 triangle_2 with 6 scans
% d = DataSpec('030220',1,6,3);    
% d.set_scan_distance(6.31e-3, 4.25, 240, 10);
% angles = [0,20,40,60,80,100];

% %%%% 1 triangle_2 with 6 scans
% d = DataSpec('030220',1,6,4);    
% d.set_scan_distance(6.31e-3, 4.25, 240, 10);
% angles = [0,20,40,60,80,140];

% %%%% 1 triangle_2 with 7 scans
% d = DataSpec('030220',1,7,2);    
% d.set_scan_distance(6.31e-3, 4.25, 240, 10);
% angles = [0,20,40,60,80,120,140];

% %%%% 1 triangle_2 with 7 scans
% d = DataSpec('030220',1,7,3);    
% d.set_scan_distance(6.31e-3, 4.25, 240, 10);
% angles = [0,20,40,60,80,100,120];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 4 triangles_1
d = DataSpec('030320',4,8,1);    
d.set_scan_distance(6.31e-3, 4.25, 240, 10);
angles = [0,20,40,60,80,100,120,140];

% %% 4 triangles_2
% d = DataSpec('030320',4,8,2);    
% d.set_scan_distance(6.31e-3, 4.25, 240, 10);
% angles = [0,20,40,60,80,100,120,140];

% %%% 8 triangles_1
% d = DataSpec('030620',8,8,1);    
% d.set_scan_distance(6.31e-3, 4.25, 240, 10);
% angles = [0,20,40,60,80,100,120,140];

% %%% 8 triangles_2
% d = DataSpec('030620',8,8,2);    
% d.set_scan_distance(6.31e-3, 4.25, 240, 10);
% angles = [0,20,40,60,80,100,120,140];

% %%% dics
% d = DataSpec('071019',4,8,2);    
% d.set_scan_distance(41.67e-3, 4.25, 240, 10);
% angles = [0,20,40,60,80,100,120,140];
 
%% ===== Scan line object ====== %%
% Read clpconfig for CLP parameter setting.
lines = d.get_clpsecm_data();
%lines.downsample(4);
lines.zeroing();
lines.params.angles = ProbeParam(angles);
 
%% ===== Plot data ===== %%
figure();
bplines   = ScanLines(lines);
p0        = ProbeParams(0);
p0.angles = ProbeParam(angles);
bplines.params = p0;
bpimage = bplines.back_project();
subplot(131); lines.plot_lines();   title('Lines');
subplot(132); bpimage.draw_image(); title('Back project image');
subplot(133); lines.plot_psf();     title('Point-spread function');

%% ===== Setup Signal Assumptions and Probe Parameters ===== %%
% -Setup disc D
% p        = ProbeParams(0);     % Initialize 0 parameters.
% p.angles = ProbeParam(angles); % Set constant parameter angles.

% Generate equilateral triangle
L = 0.45; % triangle side length %0.25 for 1 triangle data
deltaTheta = (2*pi/3) / 40; % discretization in angle

trianglefunc = @(x,y) (y > -L/(2*sqrt(3))) & (sqrt(3)*x/2 + y/2 < L/(2*sqrt(3))) & (-sqrt(3)*x/2 + y/2 < L/(2*sqrt(3))); 
theta = 0:deltaTheta:(2*pi/3);
k = length(theta);
rotatedTriangleFunc = cell(1,k);

for i = 1:length(theta)
    rotatedTriangleFunc{i} = @(x,y) trianglefunc( x * cos(theta(i)) + y * sin(theta(i)), y * cos(theta(i)) - x * sin(theta(i)) );
end

D = DictProfileArray(lines.ticks, rotatedTriangleFunc);

%% ====== Plot Dictionaries ===== %%
figure();
clf;
for i = 1:length(theta)
    imagesc(D.images{i});
    axis equal;
    pause(.1);
end

%% ====== Signal Reconstruction ===== %%
% Generate problem with formulation 'Calibrated Lasso'
p = deepcopy(lines.params);
% p.intensity = ProbeParam(ones(lines.nlines,1),ones(lines.nlines,1)*[0.8,1.2],...
%               @(v)ones(lines.nmeasures,1)*v');
 
prb = CalibLassoMultiMotif(lines,D,p);

alg = ReweightIPalmSecmRealdataMultiMotif(prb,D,lines);
%alg = IPalmSecmRealdata(prb,D,lines);
alg.set_maxiter(50);
figure; 
alg.solve();
