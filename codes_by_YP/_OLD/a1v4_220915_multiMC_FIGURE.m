%%
clear; clc;
% close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% File name
filename = 'Gmode_s3MC250SNRn9.mat';
% filename = 'emode_s3MC100SNRn15.mat';
% filename = 'Cmode_s2MC100SNRn15.mat';
load(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveCharVar = who('outputs*');
saveCharVar(6:8) = [];

%% Figure
figure;
set(gcf,'position',[750,200,560,560*1.26]);
set(gca,'position',[0.157,0.122,0.748,0.803])
hold on;

% output load and plot
for n_output=1:length(saveCharVar)

    dataLoadchar = ['dataLoad = ',char(saveCharVar(n_output))];
    eval(dataLoadchar);

    for ind=1:length(SNRs)
        if n_output == 1 || n_output == 2
            rmseSNR(ind) = sqrt( mean(dataLoad(1+(ind-1)*NmonteCarlo:NmonteCarlo+(ind-1)*NmonteCarlo))*180/pi*180/pi );
        else
            totET = [];
            for index=1:NmonteCarlo
                totET = [totET;dataLoad((ind-1)*NmonteCarlo+index).error];
            end
            Nout = 0.0; % Portion of Outliers, (ignore)
            totET = sort(abs(totET));
            rmseSNR(ind) = sqrt(mean(power(totET(1:length(totET)-floor(length(totET)*Nout)),2)));
        end
    end

    plotColor = turbo(length(saveCharVar));
    figConChar = ['h',num2str(n_output),'= plot(SNRs,rmseSNR,''color'',plotColor(',num2str(n_output),',:),''linewidth'',1.8,''displayname'',''',char(saveCharVar(n_output)),''');'];
    eval(figConChar)

end
hold off;

box on; grid on;

axis([-6 36 7e-3 1e1])
set(gca,'fontsize',24,'TickLabelInterpreter','latex','yscale','log')

xlabel('SNR~[dB]','interpreter','latex')
ylabel('RMSE~[$^\circ$]','interpreter','latex')
% legend((figH),'location','best','interpreter','latex',...
%     'position',[0.16822,0.135136904489426,0.45207211630685,0.194964286259242])
