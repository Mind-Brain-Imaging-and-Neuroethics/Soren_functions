function TTVplot(datattv,pdifs,pseudottv,varargin)

ttv1 = nan_mean(nan_mean(datattv{1}(:,:,:),1),3);
ttv2 = nan_mean(nan_mean(datattv{2}(:,:,:),1),3);

xax = linspace(1,2*length(ttv1),length(ttv1));

if CheckInput(varargin,'Xaxis')
    xax = EasyParse(varargin,'Xaxis');
end

tpoints = xax;

if ~EasyParse(varargin,'Plot','off')
    if ~EasyParse(varargin,'Plot','all')
        if ~EasyParse(varargin,'PlotPseudo','true')
            
            plot(xax,ttv1,'b');
            
            hold on;
            plot(xax,ttv2,'r');
        else
            plot(xax,ttv1,'b')
            hold on
            plot(xax,ttv2,'r')
            plot(xax,nan_mean(nan_mean(((pseudottv{1}(:,:,:)-mean(pseudottv{1}(:,:,:),2))*100)./mean(pseudottv{1}(:,:,:),2),1),3),'b--');
            hold on
            plot(xax,nan_mean(nan_mean(((pseudottv{2}(:,:,:)-mean(pseudottv{2}(:,:,:),2))*100)./mean(pseudottv{2}(:,:,:),2),1),3),'r--');            
        end
        
        patchindex = pdifs < 0.05;
        yl = ylim;
        if ~EasyParse(varargin,'SigDifs','off')
            patchstep = patchindex*(yl(2)-yl(1));
            patchstep = patchstep + yl(1);
            area(xax,patchstep,yl(1),'FaceAlpha',0.3,'LineStyle','none','FaceColor',[0.5 0.5 0.5])
%             if CheckInput(varargin,'Legend')
%                 legend([EasyParse(varargin,'Legend') {'Significant Differences'}])
%             else
%                 legend({'Condition 1','Condition 2','Significant Differences'})
%             end
            ylim(yl)
        else
%             if CheckInput(varargin,'Legend')
%                 legend([EasyParse(varargin,'Legend') {'Significant Differences'}])
%             else
%                 legend({'Condition 1','Condition 2','Significant Differences'})
%             end
        end
        xlabel('Time (ms)')
        if CheckInput(varargin,'NormRange') || CheckInput(varargin,'Pseudotrial') || CheckInput(varargin,'PseudotrialSplit')
            if CheckInput(varargin,'SlidingWindow')
                ylabel('SD of % change in TTV')
            else
                ylabel('% change in trial-to-trial SD')
            end
        else
            if CheckInput(varargin,'SlidingWindow')
                ylabel('SD of TTV')
            else
                ylabel('Trial-to-trial SD')
            end
        end
    else
        for q = 1:length(fields)
            
            plot(xax,nan_mean(datattv{1}(:,:,q),1));
            
            hold on;
            plot(xax,nan_mean(datattv{2}(:,:,q),1));
            %         patchindex = pdifs < 0.05;
            %         yl = ylim;
            %         patchstep = patchindex*(yl(2)-yl(1));
            %         patchstep = patchstep + yl(1);
            %         %    for c = 1:length(patchstep)
            %         %        patchstep2(((c-1)*2)+1) = patchstep(c);
            %         %        patchstep2(c*2) = patchstep(c);
            %         %    end
            %         area(xax,patchstep,yl(1),'FaceAlpha',0.3,'LineStyle','none','FaceColor',[0.5 0.5 0.5])
%             if CheckInput(varargin,'Legend')
%             legend([EasyParse(varargin,'Legend') {'Significant Differences'}])
%             else
%                 legend({'Condition 1','Condition 2','Significant Differences'})
%             end
            %        ylim(yl)
            xlabel('Time (ms)')
            ylabel('Trial-to-trial SD')
        end
    end
    
    %    for c = 1:length(ttv1)
    %        patch([t(patchindex) fliplr(t(patchindex))], [ones(size(ttv1(patchindex)))*yl(1) ones(size(ttv1(patchindex)))*yl(2)], [0.6 0.4 0.9], 'FaceAlpha',0.3, 'EdgeColor','none')
    %    end
end