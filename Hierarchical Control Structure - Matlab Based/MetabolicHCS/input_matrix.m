function U = input_matrix(PARA_C4b)
light = 1;
temp = 2;
RH = 3;
U = zeros(PARA_C4b.batch_shift*24*3600/PARA_C4b.sample_time*PARA_C4b.n,4);
switch PARA_C4b.usePFC
    case 1
        for i = 1:size(U,1)
            time = rem(i*PARA_C4b.sample_time/3600,24);
            if time <= 5 && time ~= 0
                U(i,light) = 450e-6;
                U(i,temp) = 299;
                U(i,RH) = 0.5;
                U(i,4) = 21;
                U(i,5) = 1;
            elseif time <=16 && time ~= 0
                U(i,light) = 450e-6;
                U(i,temp) = 299;
                U(i,RH) = 0.5;
                U(i,4) = 21;
                U(i,5) = 2;
            else
                U(i,light) = 450e-6;
                U(i,temp) = 299;
                U(i,RH) = 0.5;
                U(i,4) = 21.5;
                U(i,5) = 2;
        %         U(i,light) = 0e-6;
        %         U(i,temp) = 293;
        %         U(i,RH) = 0.7;
             end
        end
    case 2
        for i = 1:size(U,1)
            time = rem(i*PARA_C4b.sample_time/3600,24);
            if time <=16 && time ~= 0
                U(i,light) = 450e-6;
                U(i,temp) = 299;
                U(i,RH) = 0.5;
                U(i,4) = 21;
                U(i,5) = 2;
            else
                U(i,light) = 0e-6;
                U(i,temp) = 295;
                U(i,RH) = 0.7;
                U(i,4) = 21.2;
                U(i,5) = 2;
            end
        end
end