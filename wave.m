% This code is not complete and is just for test.
% The complete version will be updated after the paper been published.
 
s2=input('number of the plates'); % Input the maxinum number of plates.( more than 2)
d=input('draft');                 % Input the draft of plates.
b=input('distance');              % Input the distance from the plates to the origin.
om=input('frequency');            % Input the range of frequencies.
s1=input('size of frequency');             % Input the sample size of the frequency.
maxm=input('truncation n=');      % Input the truncation.
 
R0=zeros(s2-2,s1);
T0=zeros(s2-2,s1);
A0=zeros(s2-2,s1);
B0=zeros(s2-2,s1);
degA0=zeros(s2-2,s1);
degB0=zeros(s2-2,s1);
degR0=zeros(s2-2,s1);
degT0=zeros(s2-2,s1);
 
 
 
h=0.321;
bf0=1;
g=9.80665;
for jj=3:s2
for i5=1:s1
    jp0=om(i5);
v=jp0^2/g;
F1=@(x)x*tanh(x*h)-v;
 
k=1:maxm+1;
A=zeros(jj-1,maxm+1);         
B=zeros(jj-1,maxm+1);
P=zeros(2*(jj-1)*(maxm+1),2*(jj-1)*(maxm+1));     
Q=zeros(1,2*(jj-1)*(maxm+1));    
K0=fzero(F1,1); 
k(1)=K0;
 

 

