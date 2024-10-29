function out=pnseq11(N)
m=11;  % No. of stages
%length of the sequence



% taps = [1 (((sign(randn(1,m-1)))+1)./2)];

taps=[1 0 1 1 1 0 0 1 0 1 1];
for J=1:N   
    out(J)=taps(m);
    temp2=xor(taps(11),taps(9));
    temp=taps(1:length(taps)-1);
    taps=[temp2 temp];
end