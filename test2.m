clear;close all;clc
initpkg

%% ===== Read line data ===== $$
% Assign dataset
% d = DataSpec('081819',3,7,1);    
% d.set_scan_distance(6.3e-3, 3.3, 3, 240);
% angles = [0,20,40,60,80,100,120];

% d = DataSpec('071019',4,6,6);    
% d.set_scan_distance(6.31e-3, 4.25, 1.54, 240);
% angles = [0,20,40,60,80,100];

% d = DataSpec('071019',4,8,5);    
% d.set_scan_distance(425e-3, 4.25, 1.54, 240);
% angles = [0,20,40,60,80,100,120,140];

% d = DataSpec('071019',4,8,2);    
% d.set_scan_distance(41.67e-3, 4.25, 1.54, 240);
% angles = [0,20,40,60,80,100,120,140];

% %%%% 1 triangle
% d = DataSpec('030220',1,7,1);    
% d.set_scan_distance(6.31e-3, 4.25, 1.6889, 240);
% angles = [0,20,40,60,80,100,120];

%%% 4 triangles
d = DataSpec('030320',4,8,1);    
d.set_scan_distance(6.31e-3, 4.25, 240, 1.481);
angles = [0,20,40,60,80,100,120,140];

%% ===== Scan line object ====== %%
% [Note] When input w/o parmas, read clpconfig for CLP parameter setting.
lines = d.get_clpsecm_data();
lines.downsample(4);
lines.zeroing();
lines.params.angles = ProbeParam(angles); 

%% ===== Setup Signal Assumptions and Probe Parameters ===== %%
% Setup disc D
% Generate equilateral triangle
L = 0.1; % triangle side length 
deltaTheta = (2*pi/3) / 40; % discretization in angle

trianglefunc = @(x,y) (y > -L/(2*sqrt(3))) & (sqrt(3)*x/2 + y/2 < L/(2*sqrt(3))) & (-sqrt(3)*x/2 + y/2 < L/(2*sqrt(3))); 
theta = 0:deltaTheta:(2*pi/3);
k = length(theta);
rotatedTriangleFunc = cell(1,k);

for i = 1:length(theta)
    rotatedTriangleFunc{i} = @(x,y) trianglefunc( x * cos(theta(i)) + y * sin(theta(i)), y * cos(theta(i)) - x * sin(theta(i)) );
end

D = DictProfileArray(lines.ticks, rotatedTriangleFunc);
%% ===== Plot data ===== %%
figure();
bplines   = ScanLines(lines);
bplines.currents = circshift(bplines.currents, 0, 1);
p0        = ProbeParams(0);
p0.angles = ProbeParam(angles);
% p0.centerbias.set_value(lines.params.centerbias.value)
bplines.params = p0;
bpimage = bplines.back_project();
subplot(131); lines.plot_lines();   title('Lines');
subplot(132); bpimage.draw_image(); title('Back project image');
subplot(133); lines.plot_psf();     title('Point-spread function');

pause();

%% ====== Signal Reconstruction ===== %%
% Generate algorithm 'IPalm' and solve the simulated prob lem.
p = deepcopy(lines.params);
lda = 0.2*max(pos(D*back_project(lines))); % 0.05*
prb = CalibLasso(lines,D,p,lda);
alg = ReweightIPalmSecmRealdata(prb,D,lines);
% alg = IPalmSecmRealdata(prb,D,lines);
alg.set_maxiter(500);
figure; 
alg.solve(); 
