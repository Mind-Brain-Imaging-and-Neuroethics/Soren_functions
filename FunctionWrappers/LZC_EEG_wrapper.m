function [LZCOut] = LZC_EEG_handle(EEG)

LZCOut = zeros(1,EEG.nbchan);

disp(' ')
disp('Computing Lempel-Ziv Complexity...')

for c = 1:EEG.nbchan
    fprintf([num2str(c) ' ']);
    [~,LZCOut(c)] = lzcomplexity_tramas(EEG.data(c,:),'mediana',2,2500,0.9);
end