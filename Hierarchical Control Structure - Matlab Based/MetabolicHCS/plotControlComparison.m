%% Figure for Publication
% Author: Carles Ciurans
% Date: 22/11/2021, UCA, Clermont-Ferrand, France

%% Data Reference preparation
index0 = 1*24*3600/PARA_C4b.sample_time;     % 1 Day
fontSize = 15;
lgdfontSize = 20;
LineWidth = 2;
load ('24hDMCAPerturbationInput')
U1 = U;
Ys1 = Ys;
Y_MPC1 = Y_MPC;
Y1 = Y;
U_MPC1 = U_MPC;

load ('24hDMCAPerturbationInputOffset')
U2 = U;
Ys2 = Ys;
Y_MPC2 = Y_MPC;
Y2 = Y;
U_MPC2 = U_MPC;

figControl = figure;
        hold on
        %box on
        p2 = plot(C4b_matrices.time(2:index0+1)/3600,U1(1:index0,4),'r','LineWidth',LineWidth);
        p3 = plot(C4b_matrices.time(2:index0+1)/3600,Y1(2:index0+1,37)*100*PARA_C4b.Air_Con,'b','LineWidth',LineWidth,'LineStyle', '--'); 
        p4 = plot(C4b_matrices.time(2:index0+1)/3600,Y2(2:index0+1,37)*100*PARA_C4b.Air_Con,'b','LineWidth',LineWidth,'LineStyle', ':'); 
        hold off
        ylabel('Concentration(%)','FontSize', fontSize)
        xlabel('Time (hr)','FontSize', fontSize)
        legend('Ref O_2', 'Measurement O_2 with disturbance', 'Measurement O_2 without disturbance','Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', fontSize);
        title('Control performance comparison', 'FontSize', fontSize)
        ax = gca;
        ax.YColor = 'k';
        set(gca,'FontSize', fontSize)
        set(gca, 'fontname', 'times')
        xlim([0 24])
        xticks([0 2 4 6 8 10 12 14 16 18 20 22 24])
        legend boxoff
    