function resolve(A,Z,Q,status,time_NCP,i)


if status == 1 
    A.time_NCP(i) = time_NCP;
else
    A.time_NCP(i) = 0;
end
j = 1;


if i == 1
    while status == 0
    j = j+1;
    Z_new = change_initial_guess(A,Z,Q);
    tic
    [A.z(:,i,index),~,~,~,status] = pathmcp(Z_new,l,u,fun);
    time_NCP = toc;
    if status == 1 
        A.time_NCP(i) = time_NCP;
    else
        A.time_NCP(i) = 0;
    end
        if j>=60
            error('Path can not found the solution, change your initial guess');
        end
    end

else
    while status == 0
        j = j+1;
        Z_new = change_guess(A,Z,Q);
        tic
        [A.z(:,i,index),~,~,~,status] = pathmcp(Z_new,l,u,fun);
        time_NCP = toc;
        if status == 1 
            A.time_NCP(i) = time_NCP;
        else
            A.time_NCP(i) = 0;
        end
        if j>=60
            error('Path can not found the solution, change your initial guess');
        end
    end
end


end