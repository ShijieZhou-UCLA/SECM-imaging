% Visulization of angle experiments

%%%
% 3scans: 0,20,40
% 4scans: 0,20,40,60
% 5scans: 0,20,40,60,80
% 6scans: 0,20,40,60,80,100
% 7scans: 0,20,40,60,80,100,120
% 8scans: 0,20,40,60,80,100,120,140
%%%
scans = [3,4,5,6,7,8];
angles1 = [53.33, 58.69, 69.04, 75.37, 84.53, 94.06];
angles2 = [90, 96.02, 95.76, 90, 95.06, 94.06];
figure;
plot(scans,angles1)
hold on
plot(scans,angles2)
xticks(3:1:8)
title('Orientation angles vs number of scans')
xlabel('number of scans')
ylabel('orientation angles')
legend('Group1 (0,20,40)','Group2 (0,60,120)')


% %%%
% % 3scans: 0,60,120
% % 4scans: 0,60,120,140
% % 5scans: 0,20,60,120,140
% % 6scans: 0,20,40,60,120,140
% % 7scans: 0,20,40,60,80,120,140
% % 8scans: 0,20,40,60,80,100,120,140
% %%%
% s = [3,4,5,6,7,8];
% a = [90, 96.02, 95.76, 90, 95.06, 94.06];
% figure;
% plot(s,a)
% yticks(50:5:100)
% xlabel('number of scans')
% ylabel('angles')
