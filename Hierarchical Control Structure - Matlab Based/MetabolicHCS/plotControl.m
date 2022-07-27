%% Figure for Publication
% Author: Carles Ciurans
% Date: 21/11/2021, Bus Lyon-Clermont Ferrand

%% Data Reference preparation
index0 = 1*24*3600/PARA_C4b.sample_time;     % 1 Day
fontSize = 15;
lgdfontSize = 20;
LineWidth = 2;


figControl = figure;
  subplot(3,1,1)
        hold on
        %box on
        p2 = plot(C4b_matrices.time(2:index0+1)/3600,U(1:index0,4),'r','LineWidth',LineWidth);
        p3 = plot(C4b_matrices.time(2:index0+1)/3600,Ys(1:index0)*100*PARA_C4b.Air_Con,'-r','LineWidth',LineWidth, 'LineStyle', '--');
        p4 = plot(C4b_matrices.time(2:index0+1)/3600,Y_MPC(2,1:index0)*100*PARA_C4b.Air_Con,':r','LineWidth',LineWidth);
        p5 = plot(C4b_matrices.time(2:index0+1)/3600,Y(2:index0+1,37)*100*PARA_C4b.Air_Con,'b','LineWidth',LineWidth); 
        hold off
        ylabel('Concentration(%)','FontSize', fontSize)
        legend('Ref O_2', 'SSTO Ref O_2','Prediction MPC O_2','Measurement O_2', 'Location', 'eastoutside', 'FontSize', fontSize);
        title('O_{2} Control: Controlled Variable', 'FontSize', fontSize)
        ax = gca;
        ax.YColor = 'k';
        set(gca,'FontSize', fontSize)
        set(gca,'fontname','times')
        xlim([0 24])
        xticks([0 2 4 6 8 10 12 14 16 18 20 22 24])
        legend boxoff
        
        subplot(3,1,2)
        %box on
        p1 = plot(C4b_matrices.time(2:index0+1)/3600,U_MPC(1,1:index0)*1000/60,'g','LineWidth',LineWidth);
        ylabel('Flow (mL/min)', 'FontSize', fontSize);
        legend('Gas Flow', 'Location', 'eastoutside', 'FontSize', fontSize)
        title('O_{2} Control: Manipulated Variable', 'FontSize', fontSize)
        ax = gca;
        ax.YColor = 'k';
        set(gca,'FontSize', fontSize)
        set(gca,'fontname','times')
        xlim([0 24])
        xticks([0 2 4 6 8 10 12 14 16 18 20 22 24])
        legend boxoff

        subplot(3,1,3)
        p7 =  plot(C4b_matrices.time(2:index0+1)/3600, Y(2:index0+1,38)*100*PARA_C4b.Air_Con , 'b','LineWidth',LineWidth);
        ylabel('Concentration(%)','FontSize', fontSize)
        hold on 
        ylim([0 0.2])
        yyaxis right
        p6 = plot(C4b_matrices.time(2:index0+1)/3600, MV(2:index0+1,1) ,'g', 'LineWidth',LineWidth);
        ylim([0 2e-4])
        title('CO_{2} Control: Manipulated Variable', 'FontSize', fontSize)
        ylabel('Flow (L/s)')
        xlabel('Time (h)')
        xlim([0 24])
        xticks([0 2 4 6 8 10 12 14 16 18 20 22 24])
        legend('Measurement CO_{2}', 'CO_{2} Flow','Prediction MPC O_2','Measurement O_2', 'Location', 'eastoutside', 'FontSize', fontSize);
        legend boxoff
        ax = gca;
        ax.YColor = 'k';
        set(gca,'FontSize', fontSize)
        set(gca,'fontname','times')
        
        %         ylabel('Rate (g/h');
%         yyaxis right
%         hold on
%         p7 = plot(C4b_matrices.time(2:index0+1)/3600,C4b_matrices.Cint(2:index0+1,1),'-r');
%         p8 = plot(C4b_matrices.time(2:index0+1)/3600,C4b_matrices.Cint(2:index0+1,2),'-r');
%         p9 = plot(C4b_matrices.time(2:index0+1)/3600,C4b_matrices.Cint(2:index0+1,3),'-r');
%         p10 = plot(C4b_matrices.time(2:index0+1)/3600,C4b_matrices.Cint(2:index0+1,4),'-r');
%         hold off
%         ylabel('Concentration(mol m^{-3})')
%         legend('Net Photosynthesis Rate', 'Group 1: Internal CO_2', ...
%             'Group 2: Internal CO_2','Group 3: Internal CO_2', ...
%             'Group 4: Internal CO_2','Location', 'eastoutside');
%         xlabel('Time(h)')
%         title('Photosynthesis')