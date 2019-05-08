function [PLEOut] = PLE_EEG_handle(EEG)

disp(' ')
disp('Computing power law exponent...')

PLEOut = PedroPLE(EEG.data);