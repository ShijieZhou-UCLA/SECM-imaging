initpkg

%% ===== Read line data ===== $$
% Assign dataset

% d = DataSpec('071019',4,8,1);   % 6.31  
% d.set_scan_distance(6.31e-3, 4.25, 240);
% angles = [0,20,40,60,80,100,120,140];
% psf: [.6, .1, 5, 2, 1.5, 8, -5]
% 
d = DataSpec('071019',4,8,2);   % 41.67 
%d.set_scan_distance(41.67e-3, 4.25, 240);
d.set_scan_distance(41.67e-3, 4.25, 240, 1.54);
angles = [0,20,40,60,80,100,120,140];
% psf: [.6, .1, 5, 2, 1.5, 8, 0]
% 
% d = DataSpec('071019',4,8,3);   % 121.36
% d.set_scan_distance(121.36e-3, 4.25, 240, 2);
% angles = [0,20,40,60,80,100,120,140];
% psf: [.6, .1, 5, 2, 1.5, 8, 10]
% 
% d = DataSpec('071019',4,8,4);   % 223
% d.set_scan_distance(223e-3, 4.25, 240);
% angles = [0,20,40,60,80,100,120,140]; 
% psf: [.6, .1, 5, 2, 1.5, 8, 13] 
% 
% d = DataSpec('071019',4,8,5);   % 425
% d.set_scan_distance(425e-3, 4.25, 240, 1);
% angles = [0,20,40,60,80,100,120,140] + 180;
% psf: [.6, .1, 5, 2, 1.5, 5, 27]
% 
% d = DataSpec('071019',4,6,6);   % 6.31 new    
% d.set_scan_distance(6.31e-3, 4.25, 240);
% angles = [0,20,40,60,80,100];
% psf: [.6, .1, 5, 2, 1.5, 8, -20]
%
% d = DataSpec('071019',4,7,7);   % 700
% d.set_scan_distance(700e-3, 7, 240);
% angles = [0,20,40,60,80,100,120];
% 
% d = DataSpec('071019',4,7,8);   % 1000   
% d.set_scan_distance(1000e-3, 4.25, 240);
% angles = [0,20,40,60,80,100,120,140];
% 
% d = DataSpec('21119',0,7,1);   % 6.31   
% d.set_scan_distance(9e-3, 4.25, 240);
% angles = [0,30,60,90,120,150,180];

% d = DataSpec('022820',8,8,1);   % 6.31
% d.set_scan_distance(6.31e-3, 5.97, 200, 10);
% angles = [0,20,40,60,80,100,120,140];

% d = DataSpec('022820',8,8,2);   % 224
% d.set_scan_distance(224e-3, 4.25, 240, 1);
% angles = [0,20,40,60,80,100,120,140];

%% ===== Scan line object ====== %%
lines = d.get_clpsecm_data();
lines.zeroing();
lines.params.angles = ProbeParam(angles); 

%% ===== Setup Signal Assumptions and Probe Parameters ===== %%
% Setup disc D
disc_radius = 0.075;
func = @(x,y) x.^2 + y.^2 < disc_radius^2;
D = DictProfile(lines.ticks, func);

%% ===== Plot data ===== %%
figure();
bplines   = ScanLines(lines);
bplines.currents = circshift(bplines.currents, 0, 1);
p0        = ProbeParams(0);
p0.angles = ProbeParam(angles);
bplines.params = p0;
bpimage = bplines.back_project();
subplot(131); lines.plot_lines();   title('Lines');
subplot(132); bpimage.draw_image(); title('Back project image');
subplot(133); lines.plot_psf();     title('Point-spread function');

%pause();

%% ===== Ground Truth Image ===== %%
x_truth = [   -0.5130   -0.1390
    0.0250    0.3620
    0.1960   -0.3540
    0.2920    0.1310
];
theta = 1.15 * pi;
x_truth = x_truth * [cos(theta),-sin(theta);sin(theta),cos(theta)];

img_truth = SecmImage(lines.ticks);
id_truth = floor((x_truth - lines.ticks(1)) / lines.resolution);
for k = 1:size(x_truth,1)
    img_truth.image(id_truth(k,1), id_truth(k,2)) = 1;
end
img_truth = img_truth * D;
figure
img_truth.draw_image()

%% ===== Reset psf parameters ===== %%
psf_params = [.6, .1, 5, 2, 1.5, 8, 10]';
lines.params.psf.set_value(psf_params);
lines.params.psf.set_bound(psf_params);


%% ====== Signal Reconstruction ===== %%
% Generate algorithm 'IPalm' and solve the simulated problem.
p = deepcopy(lines.params);
lda = 0.02*max(pos(D*back_project(lines))); % 0.05*
% lda = 0;
prb = CalibLasso(lines,D,p,lda);
alg = ReweightIPalmSecmRealdata(prb,D,lines);
% alg = IPalmSecmRealdata(prb,D,lines);
alg.set_maxiter(8);
figure; 
alg.solve(); 

%% ====== Calculate Square Loss ===== %%
img_recovered = alg.vars{1} * D;
img_recovered.normalize();
loss = img_recovered.min_square_loss(x_truth,D);

