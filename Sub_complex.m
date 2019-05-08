function [out] = Sub_complex(compin,realn)

len = abs(compin);

newLen = len-realn;

ang = angle(compin);

out = newLen.*exp(1i*ang);





