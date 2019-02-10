function stab = calculateStability(SS, pv )
%0 - unstable
%1 - stable

stab = false(1,size(SS,1));

for i = 1:size(SS,1)
   Jeval = Jac(SS(i,1), SS(i,2));
   %check the eigenvectors
   if all(real(eig(Jeval)) < 0) %the point is stable
      stab(i) = true;
   end
end

    function J = Jac(T,E)
       
        J = zeros(2,2);
        J(1,1) = pv(1)*(1-pv(2)*T)-pv(2)*pv(1)*T-pv(3)*E;
        J(1,2) = -pv(3)*T;
        J(2,1) = pv(4)*E/(pv(5)+T)-pv(4)*E*T/(pv(5)+T)^2-pv(6)*E;
        J(2,2) = pv(4)*T/(pv(5)+T)-pv(6)*T-pv(7);
    end


end

