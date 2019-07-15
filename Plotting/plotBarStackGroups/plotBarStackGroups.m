function [h,barpos] = plotBarStackGroups(stackData, groupLabels,cols)
%% Plot a set of stacked bars, but group them according to labels provided.
%%
%% Params: 
%%      stackData is a 3D matrix (i.e., stackData(i, j, k) => (Group, Stack, StackElement)) 
%%      groupLabels is a CELL type (i.e., { 'a', 1 , 20, 'because' };)
%%      cols is a cell array of colours to input to each stack group - each cell should contain
%%      a NElements x Nstacks x 3 array, ordered from bottom up
%%
%% Copyright 2011 Evan Bollig (bollig at scs DOT fsu ANOTHERDOT edu
%%
%% Modified by SWT 11/07/2019
%%
%% 
NumGroupsPerAxis = size(stackData, 1);
NumStacksPerGroup = size(stackData, 2);

if nargin < 3
   cols = []; 
else
   cols = cat(4,cols{:});
end

% Count off the number of bins
groupBins = 1:NumGroupsPerAxis;
MaxGroupWidth = 0.65; % Fraction of 1. If 1, then we have all bars in groups touching
groupOffset = MaxGroupWidth/NumStacksPerGroup;
%figure removed by SWT for plotting into subplots
    hold on; 
for i=1:NumStacksPerGroup

    Y = squeeze(stackData(:,i,:));
    
    % Center the bars:
    
    internalPosCount = i - ((NumStacksPerGroup+1) / 2);
    
    % Offset the group draw positions:
    groupDrawPos = (internalPosCount)* groupOffset + groupBins;
    barpos(i,:) = groupDrawPos;
    
    h(i,:) = bar(Y, 'stacked');
    set(h(i,:),'BarWidth',groupOffset);
    set(h(i,:),'XData',groupDrawPos);
    if ~isempty(cols)
    set(h(i,:),'FaceColor','flat')
        for c = 1:length(h(i,:))
            set(h(i,c),'CData',squeeze(cols(c,i,:,:))');
        end
    end
    set(h(i,:),'LineStyle','none')
end
hold off;
set(gca,'XTickMode','manual');
set(gca,'XTick',1:NumGroupsPerAxis);
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',groupLabels);
end 
