function [x,fval,flag,j]=fzeromanysols(f,range);

if size(range,2)==1
   range=range';
end
j=0;
for i=1:size(range,2)
    ff(i)=f(range(i));
    signf(i)=sign(ff(i));
    if i>1
       if (signf(i-1)==1 & signf(i)==-1) | (signf(i-1)==-1 & signf(i)==1) & ~isnan(signf(i-1).*signf(i)) & isreal(ff(i-1).*ff(i));
         j=j+1;
         x00(1,j)=range(i-1);
         x00(2,j)=range(i);
         
         [x(j),fval(j),flag(j)]=fzero(f,  [x00(1,j) x00(2,j)] );
       end
    end
end

if j==0
    x=[];
    fval=[];
    flag=-999;
end

end
