function Van = adhensive_force(A,i,index)
Van = 0;
if A.shape == 'cuboid_shape'
    if ((nnz(A.z(17:22,i,index)) > 1) &&((A.z(9,i,index)-A.z(12,i,index)) == 0))
       heg = A.dim(3);
       wid = A.dim(2);
       len = A.dim(1);
       area = 0.05*heg*wid;
       Van =  A.van_constant*area*A.unit*A.unit*A.h;
       
   elseif ((nnz(A.z(17:22,i,index)) == 1)&&((A.z(19,i,index)>0)||(A.z(20,i,index)>0))&&((A.z(9,i,index)-A.z(12,i,index)) == 0))
       heg = A.dim(3);
       wid = A.dim(2);
       len = A.dim(1);
       
       area = heg*wid;
       Van =  A.van_constant*area*A.unit*A.unit*A.h;
   elseif ((nnz(A.z(17:22,i,index)) == 1)&&((A.z(21,i,index)>0)||(A.z(22,i,index)>0))&&((A.z(9,i,index)-A.z(12,i,index)) == 0))
       heg = A.dim(3);
       wid = A.dim(2);
       len = A.dim(1);
       
       area = len*wid;
       Van =  A.van_constant*area*A.unit*A.unit*A.h;
    else
       heg = A.dim(3);
       wid = A.dim(2);
       len = A.dim(1);
       area = 0*heg*wid;
       Van =  A.van_constant*area*A.unit*A.unit*A.h;
        
    end
elseif A.shape == 'spiked_shape'
     if ((nnz(A.z(17:25,i,index)) > 1) &&((A.z(9,i,index)-A.z(12,i,index)) == 0))
       area = 0*A.dim(4)*A.dim(5);
       Van =  A.van_constant*area*A.unit*A.unit_mass*A.h;
       
   elseif ((nnz(A.z(17:25,i,index)) == 1)&&((A.z(19,i,index)>0)||(A.z(20,i,index)>0))&&((A.z(9,i,index)-A.z(12,i,index)) == 0))
       area = A.dim(4)*A.dim(5);
       
       Van =  A.van_constant*area*A.unit*A.unit_mass*A.h;
     else
       area = 0*A.dim(4)*A.dim(5);
       Van =  A.van_constant*area*A.unit*A.unit_mass*A.h;  
    end
elseif A.shape == 'spiked_ended'
     if ((nnz(A.z(17:27,i,index)) > 1)&&((A.z(9,i,index)-A.z(12,i,index)) == 0))
       area = 0*400e-6*100e-6;
       Van =  A.van_constant*area*A.unit*A.unit_mass*A.h;
       
   elseif ((nnz(A.z(17:25,i,index)) == 1)&&((A.z(19,i,index)>0)||(A.z(20,i,index)>0))&&((A.z(9,i,index)-A.z(12,i,index)) == 0))
       area = A.dim(3)*A.dim(5);
       
       Van =  A.van_constant*area*A.unit*A.unit_mass*A.h;
   elseif ((nnz(A.z(17:25,i,index)) == 1)&&((A.z(21,i,index)>0)||(A.z(22,i,index)>0)||(A.z(23,i,index)>0)||(A.z(24,i,index)>0))&&((A.z(9,i,index)-A.z(12,i,index)) == 0))
       L = sqrt((A.dim(4))^2+(A.dim(2)/2)^2);
       area = L*A.dim(5);
       Van =  A.van_constant*area*A.unit*A.unit_mass*A.h;  
     else
       area = 0*A.dim(4)*A.dim(5);
       Van =  A.van_constant*area*A.unit*A.unit_mass*A.h;
    end
    
end

end
