function [IrasaOut] = PLEIrasa_EEG_handle_2(spec)

disp(' ')
disp('Computing PLE with IRASA...')

tmp = amri_sig_plawfit(spec,[1 40]);
IrasaOut = tmp.Beta;
IrasaOut = IrasaOut';
