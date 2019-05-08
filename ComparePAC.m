function [pDifs,pac1,pac2,pPac1,pPac2] = ComparePAC(EEG,conditions,electrode,varargin)

indices = cell(1,2);

if ~CheckInput(varargin,'Frange')
    Frange = [0.5 3 7 12];
else
    Frange = EasyParse(varargin,'Frange');
end

if ~CheckInput(varargin,'Trange')
    Trange = [-200 200];
else
    Trange = EasyParse(varargin,'Trange');
end

for q = 1:length(conditions)
    for c = 1:length(EEG.event)
        if ContainsAll(EEG.event(c).type,strsplit(conditions{q},'_'))
            indices{q} = [indices{q} EEG.event(c).latency];
        end
    end
end

[pac1,pac2,~,pDifs] = erpac_corr2(EEG.data(electrode,:),EEG.srate,indices{1},indices{2},Trange,Frange(1),Frange(2),Frange(3),Frange(4));

if EasyParse(varargin,'Mcorrect','fdr')
    pDifs = mafdr(pDifs,'BHFDR',true);
elseif EasyParse(varargin,'Mcorrect','bonferroni')
    pDifs = pDifs*length(pDifs);
end

if CheckInput(varargin,'SigMask')
    pac1Mask = pac1;
    pac2Mask = pac2;
    if EasyParse(varargin,'SigMask','Differences')
        for c = 1:length(pDifs)
            if pDifs(c) > 0.05
                pac1Mask(c) = 0;
                pac2Mask(c) = 0;
            end
        end
    else
        [pac1Mask,~,pPac1] = erpac_corr(EEG.data(electrode,:),EEG.srate,indices{1},Trange,Frange(1),Frange(2),Frange(3),Frange(4),200);
        [pac2Mask,~,pPac2] = erpac_corr(EEG.data(electrode,:),EEG.srate,indices{1},Trange,Frange(1),Frange(2),Frange(3),Frange(4),200);
        for c = 1:length(pac1Mask)
           if pPac1(c) > 0.05
               pac1Mask(c) = 0;
           end
           if pPac2(c) > 0.05
              pac2Mask(c) = 0; 
           end
        end
        
        if EasyParse(varargin,'SigMask','All')
            for c = 1:length(pDifs)
                if pDifs(c) > 0.05
                    pac1Mask(c) = 0;
                    pac2Mask(c) = 0;
                end
            end
        end
    end
end


if EasyParse(varargin,'Plot','on')
    figure;
    Ts = 1000/EEG.srate;
    
    if CheckInput(varargin,'SigMask')
        plot([(Trange(1)/Ts):(Trange(2)/Ts)],pac1Mask)
        hold on
        plot([(Trange(1)/Ts):(Trange(2)/Ts)],pac2Mask)
    else
        plot([(Trange(1)/Ts):(Trange(2)/Ts)],pac1)
        hold on
        plot([(Trange(1)/Ts):(Trange(2)/Ts)],pac2)
    end
    
    legend(conditions{1},conditions{2})
    xlabel(['Time point (x' num2str(Ts) ' ms)'])
    ylabel('Circular-linear correlation coefficient')
    title(['Event related phase-amplitude coupling between ' num2str(Frange(1)) '-' num2str(Frange(2)) 'Hz and ' num2str(Frange(3)) '-' num2str(Frange(4)) 'Hz'])
end


