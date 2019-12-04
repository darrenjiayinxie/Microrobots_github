function Van = adhensive_force(A,i,index)
Van = 0;
if A.shape == 'cuboid_shape'
    % z(17:22) are l_1 to l6. l1,l2:len*hig; l3,l4:wid*hig; l5,l6:Wid*Len
    l1 = A.z(17,i,index);
    l2 = A.z(18,i,index);
    l3 = A.z(19,i,index);
    l4 = A.z(20,i,index);
    l5 = A.z(21,i,index);
    l6 = A.z(22,i,index);
    % normal force
    p_n = A.z(24,i,index);
    
    N = nnz(A.z(17:22,i,index)); % number of unzero l_i
    
    if ( N > 1) % line or point contact 
       area = 0;
       Van =  A.van_constant*area*A.unit*A.unit*A.h;
       
   elseif (N == 1)&&(l1>0||l2>0)&&(p_n>0)% surface contact  heg*len
       heg = A.dim(3);
       len = A.dim(1);
       
       area = heg*len;
       Van =  A.van_constant*area*A.unit*A.unit*A.h;
   elseif (N == 1)&&(l3>0||l4>0)&&(p_n>0)% surface contact  heg*wid
       heg = A.dim(3);
       wid = A.dim(2);

       area = heg*wid;
       Van =  A.van_constant*area*A.unit*A.unit*A.h;
   elseif (N == 1)&&(l5>0||l6>0)&&(p_n>0) % surface contact  len*wid
       
       wid = A.dim(2);
       len = A.dim(1);
       
       area = len*wid;
       Van =  A.van_constant*area*A.unit*A.unit*A.h;
    else
 
       area = 0;
       Van =  A.van_constant*area*A.unit*A.unit*A.h;
        
    end
elseif A.shape == 'spiked_shape'
     % z(17:24) are l_1 to l8. 
     l1 = A.z(17,i,index);
     l2 = A.z(18,i,index);
     l3 = A.z(19,i,index);
     l4 = A.z(20,i,index);
     l5 = A.z(21,i,index);
     l6 = A.z(22,i,index);
     l7 = A.z(23,i,index);
     l8 = A.z(24,i,index);
     p_n = A.z(26,i,index);
     N = nnz(A.z(17:24,i,index)); % number of unzero l_i
     if (N > 1) % line or point contact
       area = 0;
       Van =  A.van_constant*area*A.unit*A.unit_mass*A.h;
       
   elseif (N == 1)&&(l3>0||l4>0)&&(p_n>0)
       area = A.dim(4)*A.dim(5);
       
       Van =  A.van_constant*area*A.unit*A.unit_mass*A.h;
     else
       area = 0;
       Van =  A.van_constant*area*A.unit*A.unit_mass*A.h;  
    end
elseif A.shape == 'spiked_ended'
    % z(17:26) are l_1 to l10. 
     l1 = A.z(17,i,index);
     l2 = A.z(18,i,index);
     l3 = A.z(19,i,index);
     l4 = A.z(20,i,index);
     l5 = A.z(21,i,index);
     l6 = A.z(22,i,index);
     l7 = A.z(23,i,index);
     l8 = A.z(24,i,index);
     l9 = A.z(25,i,index);
     l10 = A.z(26,i,index);
     
     p_n = A.z(28,i,index);
     N = nnz(A.z(17:26,i,index)); % number of unzero l_i
     if (N > 1)
       area = 0;
       Van =  A.van_constant*area*A.unit*A.unit_mass*A.h;
       
   elseif (N == 1)&&(l3>0||l4>0)&&(p_n>0)
       area = A.dim(3)*A.dim(5);
       Van =  A.van_constant*area*A.unit*A.unit_mass*A.h;
   elseif (N == 1)&&(l5 > 0||l6 >0||l7>0||l8>0)&&(p_n>0)
       L = sqrt((A.dim(4))^2+(A.dim(2)/2)^2);
       area = L*A.dim(5);
       Van =  A.van_constant*area*A.unit*A.unit_mass*A.h;  
     else
       area = 0;
       Van =  A.van_constant*area*A.unit*A.unit_mass*A.h;
    end
    
end

end
