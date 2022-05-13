function [OutputPara, OutputFval] = GA(GenN, PopS,startPara, startFval,ver)

warning('error','MATLAB:ode15s:IntegrationTolNotMet');
%% 
load_para(ver,'WT');


%% ub and lb
ub = startPara.*1.5; 

lb = startPara.*0.4; 

% Multi-GA
tic
options = optimoptions('gamultiobj','MaxGenerations',GenN,'PopulationSize', PopS, ...
                       'InitialPopulation', startPara,'UseParallel',true, 'PlotFcn',@gaplotpareto);



    [val, fval] = gamultiobj(@fitness2, 42,[],[],[],[],lb,ub,options)%1/20

toc

index = 0;
direxist=1;
while direxist == 1
    filename = strcat('MultiGA_Output',num2str(index),'.mat');
    if ~exist(filename,'file')
        direxist=0;
    else
        index=index+1;
    end
end
    filename=['MultiGA_Output',num2str(index),'.mat'];
% save("MultiGA_Outputs.mat", "val", "fval");
save(filename, "val", "fval");
figurename = ['map',num2str(index),'.fig'];
saveas(gcf,figurename);
close all;
NN=length(fval(:,1));
CurrentPara = startPara;
CurrentCost = startFval;
for j=1:NN
    if fval(j,1)<=CurrentCost(1) && fval(j,2)<=CurrentCost(2)
        CurrentPara = val(j,:);
        CurrentCost = fval(j,:);
    end
end
if CurrentCost(1)<=25 && CurrentCost(2)<=1
    OutputPara = CurrentPara; OutputFval = CurrentCost;
else
    OutputPara = startPara; OutputFval=startFval;
end
end


