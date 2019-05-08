function [pdifs] = RQA_compareconds(conditions,electrode,measureID,varargin)

load('/Users/Soren/Desktop/SorenMusic/EEGData/EpochedReal/trialFlags.mat')

load('/Users/Soren/Desktop/SorenMusic/EEGData/EpochedReal/crqaVals.mat')

indices = cell(1,length(conditions));

for q = 1:length(conditions)
    for qq = 1:length(trialFlags)
        if ContainsAll(trialFlags{qq},strsplit(conditions{q},'_'))
            indices{q} = [indices{q} qq];
        end
    end
end

xAxis = linspace(-200,300,250);

% for c = 1:2
%    for cc = 1:length(crqaVals)
%         crqaVals{c,cc} = downsample(crqaVals{c,cc},2);
%         crqaVals{c,cc} = crqaVals{c,cc}(1:250,:);
%    end
% end


for q = measureID
    thing1 = squeeze(crqaVals.(electrode)(:,q,indices{1}));
    thing2 = squeeze(crqaVals.(electrode)(:,q,indices{2}));
    thing1means = mean(thing1,1);
    thing2means = mean(thing2,1);
    if ~EasyParse(varargin,'RemoveOutliers','off')
        thing1(:,find(isoutlier(thing1means))) = [];
        thing2(:,find(isoutlier(thing2means))) = [];
    end
    for c = 1:250
        pdifs(c,q) = ranksum(thing1(c,:),thing2(c,:));
    end
end

if ~EasyParse(varargin,'MCorrect','off')
    for q = measureID
        if ~any(isnan(pdifs(:,q)))
            pdifs(:,q) = mafdr(pdifs(:,q));
        end
    end
end

if ~EasyParse(varargin,'Plot','off')
    for q = measureID
        figure
        plot(xAxis,nan_mean(crqaVals.(electrode)(:,q,indices{1}),3));
        hold on;
        plot(xAxis,nan_mean(crqaVals.(electrode)(:,q,indices{2}),3));
        
        patchindex = pdifs(:,q) < 0.05;
        yl = ylim;
        patchstep = patchindex*(yl(2)-yl(1));
        patchstep = patchstep + yl(1);
        ylim(yl)
        
        area(xAxis,patchstep,yl(1),'FaceAlpha',0.3,'LineStyle','none','FaceColor',[0.5 0.5 0.5])
        legend([conditions {'Significant differences'}])
    end
end