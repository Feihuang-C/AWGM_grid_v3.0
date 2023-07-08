function [y]=bpfilt(x,dt,lf,hf);
  nyq=0.5/dt;
  wn=[lf/nyq,hf/nyq];
  [b,a]=butter(2,wn);
  y=filtfilt(b,a,x);
