function descr = R2string(R2,varargin)

R2 = round(R2,2);

if R2<1
       % descr = ['R�_{',myName,'} = ',num2str(R2)];
       descr = ['R�=',num2str(R2)];
    else
      %  descr = ['R�_{',myName,'} > 0.99'];
      descr = 'R�>0.99';
end
end