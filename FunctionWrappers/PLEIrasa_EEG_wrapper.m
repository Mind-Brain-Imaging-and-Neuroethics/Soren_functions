function [IrasaOut] = PLEIrasa_EEG_handle(EEG)

disp(' ')
disp('Computing PLE with IRASA...')

tmp = amri_sig_plawfit(amri_sig_fractal(EEG.data',EEG.srate),[0.5 50]);
IrasaOut = tmp.Beta;
IrasaOut = IrasaOut';
