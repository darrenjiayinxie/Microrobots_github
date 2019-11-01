function Z_new = change_initial_guess(A,Z)
unit = A.unit;
unit_mass = A.unit_mass;
e_t = A.ellipsoid(1);
e_o = A.ellipsoid(2);
e_r = A.ellipsoid(3)*unit;

m = A.mass*unit_mass;
g = A.gravity*unit;
h = A.h;
mu = A.cof;

if A.shape == 'spiked_shape'
    R = rand;  
    Z_new = R*Z;

    R = rand;
    if R < 0.5
        Z_new(34) = 0;
        Z_new(70) = 0;
    else
        Z_new(70) = R*Z(70);
    end
    R = rand; 
    if R < 0.5
        Z_new(35) = 0;
        Z_new(71) = 0;
    else
        Z_new(71) = R*Z(71);
    end
    R = rand;  
    if R < 0.5
        Z_new(36) = 0;
        Z_new(72) = 0;
    else
        Z_new(72) = R*Z(72);
    end
    R = rand;  
    if R < 0.5
        Z_new(37) = 0;
        Z_new(73) = 0;
    else
        Z_new(73) = R*Z(73);
    end

    R = rand;  
    if R < 0.5
        Z_new(38) = 0;
        Z_new(74) = 0;
    else
        Z_new(74) = R*Z(74);
    end


    R = rand;  
    if R < 0.5
        Z_new(39) = 0;
        Z_new(75) = 0;
    else
        Z_new(75) = R*Z(75);
    end
elseif A.shape == 'spiked_ended'
    R = rand;    
    

    R = rand;
    if R < 0.5
        Z_new(44) = 0;
        Z_new(98) = 0;
    else
        Z_new(98) = R*Z(98);
    end
    R = rand; 
    if R < 0.5
        Z_new(45) = 0;
        Z_new(99) = 0;
    else
        Z_new(99) = R*Z(99);
    end
    R = rand;  
    if R < 0.5
        Z_new(46) = 0;
        Z_new(100) = 0;
    else
        Z_new(100) = R*Z(100);
    end
    R = rand;  
    if R < 0.5
        Z_new(47) = 0;
        Z_new(101) = 0;
    else
        Z_new(101) = R*Z(101);
    end

    R = rand;  
    if R < 0.5
        Z_new(48) = 0;
        Z_new(102) = 0;
    else
        Z_new(102) = R*Z(102);
    end


    R = rand;  
    if R < 0.5
        Z_new(49) = 0;
        Z_new(103) = 0;
    else
        Z_new(103) = R*Z(103);
    end
    
    R = rand;  
    if R < 0.5
        Z_new(50) = 0;
        Z_new(104) = 0;
    else
        Z_new(104) = R*Z(104);
    end
    
    R = rand;  
    if R < 0.5
        Z_new(51) = 0;
        Z_new(105) = 0;
    else
        Z_new(105) = R*Z(105);
    end
    Z_new(52:97) = rand(46,1);
end
end