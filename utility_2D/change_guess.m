function Z_new = change_guess(A,Z)
unit = A.unit;
unit_mass = A.unit_mass;
e_t = A.ellipsoid(1);
e_o = A.ellipsoid(2);
e_r = A.ellipsoid(3)*unit;

Z_new = Z;

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
Z_new(40:69) = rand(30,1);

end



