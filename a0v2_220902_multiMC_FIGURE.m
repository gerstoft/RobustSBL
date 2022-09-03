%%
clear; clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% File name
filename = 'SNRn3SimN2.mat';
load(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varList = who('-file',filename);

%% Figure
figure;
set(gcf,'position',[750,200,560,560*1.26]);
set(gca,'position',[0.157,0.122,0.748,0.803])
hold on;

% output load and plot
for n_output=1:length(varList)-2

    dataLoadchar = ['dataLoad = ',char(varList(n_output+2))];
    eval(dataLoadchar);

    for ind=1:length(SNRs)
        totET = [];

        for index=1:NmonteCarlo
            totET = [totET;dataLoad((ind-1)*NmonteCarlo+index).error];
        end

        Nout = 0.0; % Portion of Outliers, (ignore)
        totET = sort(abs(totET));

        rmseSNR(ind) = sqrt(mean(power(totET(1:length(totET)-floor(length(totET)*Nout)),2)));
    end

    plotColor = turbo(length(varList)-2);
    figConChar = ['h',num2str(n_output),'= plot(SNRs,rmseSNR,''color'',plotColor(',num2str(n_output),',:),''linewidth'',1.8,''displayname'',''',char(varList(n_output+2)),''');'];
    eval(figConChar)

end
hold off;

box on; grid on;

% axis([-20 20 .002 13])
% set(gca,'fontsize',24,'TickLabelInterpreter','latex','yscale','log')

xlabel('SNR~[dB]','interpreter','latex')
ylabel('RMSE~[$^\circ$]','interpreter','latex')
% legend((figH),'location','best','interpreter','latex',...
%     'position',[0.16822,0.135136904489426,0.45207211630685,0.194964286259242])
