function descr = R2string(R2,varargin)

R2 = round(R2,2);

if R2<1
       % descr = ['R²_{',myName,'} = ',num2str(R2)];
       descr = ['R²=',num2str(R2)];
    else
      %  descr = ['R²_{',myName,'} > 0.99'];
      descr = 'R²>0.99';
end
end