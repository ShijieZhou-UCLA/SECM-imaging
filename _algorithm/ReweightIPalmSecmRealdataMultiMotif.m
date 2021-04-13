classdef ReweightIPalmSecmRealdataMultiMotif < ReweightIPalmMultiMotif
    properties
        D      % true dictionary of kernels
        lines  % observed lines
    end
    
    methods
        function obj = ReweightIPalmSecmRealdataMultiMotif(problem,dict,lines)
            obj = obj@ReweightIPalmMultiMotif(problem);
            obj.D = dict;
            obj.lines = lines;
        end

        function display_result(obj)
            subplot(221);
            display_result@ReweightIPalmMultiMotif(obj);
            title('Objective value');
            
            subplot(222);
            Yhat = obj.D*obj.vars{1};
            Yhat.draw_image();
            title('Current image Y');
            
            subplot(223); 
            obj.lines.plot_psf();
            title('Point spread function');
            legend off;
            
            subplot(224);
            obj.vars{1}.draw_sum_images; %%% Shijie: display the sum sparse map.
            %obj.vars{1}.draw_images; %%% Shijie: display the separate sparse maps, please uncomment line 51 & 52 when using this.
            set(gca,'YDir','normal');
            title('Current map X0');
            
%             angles = obj.vars{2}.angles.value;
%             if iscolumn(angles); angles = angles'; end
%             disp(['angles(deg): ', num2str(angles,'%.2f  ')]);
%             
%             shifts = obj.vars{2}.shifts.value;
%             if iscolumn(shifts); shifts = shifts'; end
%             disp(['shifts(mm) : ', num2str(shifts,'%.2f  ')]);
%             
%             intensity = obj.vars{2}.intensity.value;
%             if iscolumn(intensity); intensity = intensity'; end
%             disp(['intensity  : ', num2str(intensity,'%.2f  ') ]);
%             
%             psf = obj.vars{2}.psf.value;
%             if iscolumn(psf); psf = psf'; end
%             disp(['psf        : ', num2str(psf,'%.2f  ') ]);

             %pause(20); % For observe the separate sparse maps
             %close all;
        end
    end
end