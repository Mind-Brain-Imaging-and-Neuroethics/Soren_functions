function [slopes] = MixedEffectsPlot(model)

indvar = model.PredictorNames{2};
depvar = model.ResponseName;


data = model.Variables;

if isa(model,'LinearMixedModel')
    Subjects = unique(data.Subject);
    for c = 1:length(Subjects)
        currData = data(find(data.Subject == Subjects(c)),:);
        if ~isempty(currData)
            B = polyfit(currData.(indvar),currData.(depvar),1);
            slopes(c) = B(1);
            a = scatter(currData.(indvar),currData.(depvar),6,'filled');
            hold on
            b = plot(linspace(min(currData.(indvar)),max(currData.(indvar)),200),B(1)*linspace(min(currData.(indvar)),max(currData.(indvar)),200) + B(2),'--','LineWidth',2);
            b.Color = a.CData;
            hold on
        end
    end
    xlabel(indvar,'FontSize',14)
    ylabel(depvar,'FontSize',14)
    
elseif isa(model,'GeneralizedLinearMixedModel') && strcmpi(model.Link.Name,'logit')
    
    for c = 1:length(unique(data.Subject))
        currData = data(find(data.Subject == c),:);
        [coeffs,dev,stats] = glmfit(currData.(indvar),currData.(depvar),'binomial','link','logit');
        logitFit = glmval(coeffs,currData.(indvar),'logit');
        a = scatter(currData.(indvar),currData.(depvar),6,'filled');
        hold on
        b = plot(sort(currData.(indvar)),sort(logitFit),'--','LineWidth',2);
        b.Color = a.CData;
        hold on
    end
    xlabel(indvar,'FontSize',14)
    ylabel(depvar,'FontSize',14)
    
    % for c = 1:length(unique(data.Subject))
    %     currData = data(find(data.Subject == c),:);
    %
    % end
end