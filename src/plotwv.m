function plotwv(fl,tlim,txt,desamp)
[a,b,c,d]=readsac0(fl);
dt=a(1);
if nargin>=4
    dti=(dt*desamp);
    fn=1/dt/2;
    lf=0.01;
    hf=fn;
    d=bpfilt(d,dt,lf,hf);
    d=decimate(d,desamp);
end
tim=[0:length(d)-1]*dt;
plot(tim,d);
if nargin>=2
    if ~isempty(tlim)
        xlim(tlim);
    end
end

if nargin>=3
    text(tlim(2)/2,max(d)*0.8,txt);
end

xlabel('Time (seconds)');