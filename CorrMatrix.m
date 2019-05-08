function [corrMat,pMat,nMat] = CorrMatrix(input1,input2,varargin)

% if ~istable(input1) && CheckInput(varargin,'Names')
%    input = array2table(input,'VariableNames',EasyParse(varargin,'Names')); 
% elseif ~istable(input1) && nargin < 3
%     error('You must input either a table or a matrix with a cell array of variable names')
% end

if isempty(input2)
    input2 = input1;
end

if EasyParse(varargin,'Threshold','off')
    threshold = 0;
else
    threshold = 1;
end

if ~CheckInput(varargin,'Cov')
    for c = 1:width(input1)
        for cc = 1:width(input2)
            removeindex = [];
            
            thing1 = input1{:,c};
            thing2 = input2{:,cc};
            
            if ~EasyParse(varargin,'RemoveOutliers','false')
                if ~isequal(unique(input1{:,c}),[0; 1]) && ~isequal(unique(input2{:,cc}),[0; 1])
                    removeindex = [removeindex find(isnan(input1{:,c})) find(isoutlier(input1{:,c}))
                        find(isnan(input2{:,cc})) find(isoutlier(input2{:,cc}))];
                end
                removeindex = unique(removeindex);
                
                thing1(removeindex) = [];
                thing2(removeindex) = [];
            end
            
            
            
            [corrMat(c,cc) pMat(c,cc)] = corr(thing1,thing2,'Type','Spearman');
            if ((pMat(c,cc) > 0.1) || (abs(corrMat(c,cc)) < 0.2)) && threshold
                corrMat(c,cc) = 0;
            end
            nMat(c,cc) = length(thing1);
        end
    end
else
    
    for c = 1:width(input1)
        for cc = 1:width(input2)
            removeindex = [];
            covari = EasyParse(varargin,'Cov');
            
            thing1 = input1{:,c};
            thing2 = input2{:,cc};

            if ~EasyParse(varargin,'RemoveOutliers','false')
                if ~isequal(unique(input1{:,c}),[0; 1]) && ~isequal(unique(input2{:,cc}),[0; 1])
                    removeindex = [find(isnan(input1{:,c}))' find(isoutlier(input1{:,c}))' ...
                        find(isnan(input2{:,cc}))' find(isoutlier(input2{:,cc}))'];
                end
                
                for q = 1:size(covari,2)
                   removeindex = [removeindex find(isnan(covari(:,q)))' find(isoutlier(covari(:,q)))']; 
                end
                
                removeindex = unique(removeindex);
                
                thing1(removeindex) = [];
                thing2(removeindex) = [];
                covari(removeindex,:) = [];
            end
            
            [r,p] = partialcorr(horzcat(thing1,thing2),covari,'Type','Spearman');
            corrMat(c,cc) = r(2,1);
            pMat(c,cc) = p(2,1);
            if ((pMat(c,cc) > 0.05) || (abs(corrMat(c,cc)) < 0.2)) && threshold
                corrMat(c,cc) = 0;
            end
            nMat(c,cc) = length(thing1);
        end
    end
end

if EasyParse(varargin,'MCorrect','on')
    pMat = reshape(mafdr(reshape(pMat,[],1),'BHFDR',true),size(pMat,1),size(pMat,2));
    if ~EasyParse(varargin,'Threshold','off')
        corrMat = corrMat.*(pMat < 0.05);
    end
end

if ~EasyParse(varargin,'Plot','off')
    if ~CheckInput(varargin,'Colormap')
   heatmap(input2.Properties.VariableNames,input1.Properties.VariableNames,corrMat,'Colormap',parula)
    else
           heatmap(input2.Properties.VariableNames,input1.Properties.VariableNames,corrMat,'Colormap',EasyParse(varargin,'Colormap'))
    end
   %h.ColorLimits = [-1 1]; 
end