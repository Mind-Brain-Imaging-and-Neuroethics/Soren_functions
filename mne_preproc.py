def autoreject_log(ft_file,out_file):
    import mne
    import autoreject
    import json
    
    # Load the file
    #info = mne.io.read_info(raw_file)
    epochs = mne.read_epochs_fieldtrip(ft_file,info=None)
    
    # Resample the data - we're only going to use the thresholds
    #epochs.resample(400,npad='auto')
    
    # Apply autoreject to find bad channels
    ar = autoreject.AutoReject()
    ar.fit(epochs)
    reject_log = ar.get_reject_log(epochs)

    reject_log = reject_log.bad_epochs
    reject_log = reject_log.tolist()

    # Write to disk
    with open(out_file,'w') as f:
        json.dump(reject_log,f)
    
    return

def autoreject_epochs(ft_file,out_file,log_file):
    import mne
    import autoreject
    import json

    # Load the file
    #info = mne.io.read_info(raw_file)
    epochs = mne.read_epochs_fieldtrip(ft_file,info=None)
    
    # Resample the data
    #epochs.resample(500,npad='auto')

    # Apply autoreject    
    ar = autoreject.AutoReject()
    epochs,reject_log = ar.fit_transform(epochs,return_log=True)

    reject_log = reject_log.bad_epochs
    reject_log = reject_log.tolist()

    # Write to disk
    with open(log_file,'w') as f:
        json.dump(reject_log,f)
    
    #Save data to file
    epochs.save(out_file)
    return

def autoreject_threshold(ft_file,out_file,raw_file=0):
    import sys
    import mne
    import autoreject
    import json

    PY3 = sys.version_info[0] == 3

    if PY3:
        string_types = str,
    else:
        string_types = basestring,

    if isinstance(raw_file,string_types):
        info = mne.io.read_info(raw_file)
        epochs = mne.read_epochs_fieldtrip(ft_file,info)
    else:
        epochs = mne.read_epochs_fieldtrip(ft_file,info=None)
    
    reject = autoreject.get_rejection_threshold(epochs)
    
    with open(out_file,'w') as f:
        json.dump(reject,f)
        
        
    
def movement_correct(raw_file,out_file):
    import mne
    
    
    return
