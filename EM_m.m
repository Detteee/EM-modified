function [xbest, fmin,f] = EM_m(obj_f, n_particles, L, U, MAX_iter,flog,verbose)
    %% -- Parameters --
    % Local search
    n_LS_main = 24; n_LS_leaf = 16  ;
    % Random
    n_Random = fix(n_particles * 0.3); %number of Random particles
    n_Move = n_particles - n_Random; %number of particles moved by EM - force
    Random_iter = fix(MAX_iter * 0.1); %max iterations to be done with Random particles
    % alpha_r
    a_r_U = 1; %initial max probability(ratio) of repulsive force 
    a_r_U_min = 0.3; %minimun of max repulsion ratio
    a_r_U_reduction_rate = exp(log(a_r_U_min)/MAX_iter);
    
    if nargin < 6
        flog = 0;verbose = 0;
    elseif nargin < 7
        verbose = 0;
    end
    if flog
        f = zeros(1, MAX_iter);
    end

    %% -- Initialize --
    dim = length(L);
    % Randomly generate particles
    Particles = ones(n_particles, 1) * L + rand(n_particles, dim) .* (U - L);
    fvalue = zeros(n_particles, 1);
    % Calculate objective funtion value
    for i = 1:n_particles
        fvalue(i) = obj_f(Particles(i, :));
    end
    Particles = [Particles(:, 1:dim), fvalue];
    % Sort particles
    Particles = sortrows(Particles, dim + 1);


    %% -- Iteration --
    for iter = 1:MAX_iter
        %% -- Cal q --
        q_denominator = Particles(1, end) - Particles(n_Move + 1, end);
        q = (Particles(1:n_Move + 1, end) - Particles(n_Move + 1, end) * ones(n_Move + 1, 1)) / q_denominator;
        %% -- Move --
        Particles = Move(obj_f, Particles, q, n_Move, a_r_U, dim, L,U);
        %% -- Sort new particles --
        Particles = sortrows(Particles, dim + 1);
        %% -- Random phase(Move Random particles)--
        % %Move Random particles
        Particles = Random_regenerate(obj_f, Particles, dim,L,U, n_Random, n_Move);
        %% -- Sort new particles --
        Particles = sortrows(Particles, dim + 1);
        %% -- Local Search(with the best particle) --
        Particles = l_s(obj_f, Particles, dim, L,U, n_LS_main, n_LS_leaf);

        a_r_U = a_r_U * a_r_U_reduction_rate;

        if iter == Random_iter
            n_Move = n_particles - 1;
            n_Random = 1;
        end


        if verbose
            disp(" ")
            fprintf("iter : %d  ; fbest(x) = %d \n", iter, Particles(1, end));
            Particles(1, 1:dim)
            disp("----------------------------------------------------------")
        end

        if flog
            f(iter) = Particles(1, end);
        end

    end


    %% -- output --
    if flog
        semilogy(f);
    end

    xbest = Particles(1, 1:dim);
    fmin = Particles(1, end);
end

function Particles = Move(obj_f, Particles, q, n_Move, a_r_U, dim,L,U)
    a_r = rand*a_r_U;

    P_temp = Particles(1, 1:dim) + a_r * (Particles(1, 1:dim) - Particles(2, 1:dim)) * q(2);
    P_temp(P_temp > U) = U(P_temp > U);
    P_temp(P_temp < L) = L(P_temp < L);
    tempf = obj_f(P_temp);
    if Particles(1, end) > tempf
        Particles(1, :) = [P_temp,tempf];
    end

    for i = 2:n_Move
        a_r = rand * a_r_U;
        F_a = Particles(i - 1, 1:dim) - Particles(i, 1:dim);
        F_r = Particles(i, 1:dim) - Particles(i + 1, 1:dim);
        P_temp = Particles(i, 1:dim) + F_a * q(i - 1) * q(i) * (1 - a_r) ...
            +F_r * q(i) * q(i + 1) * a_r;
        P_temp(P_temp > U) = U(P_temp > U);
        P_temp(P_temp < L) = L(P_temp < L);
        tempf = obj_f(P_temp);

        % Center
        if Particles(i, end) < tempf
            P_temp = (P_temp + F_a + F_r + Particles(i, 1:dim))/4;
            tempf = obj_f(P_temp);
        end

        Particles(i,:) = [P_temp,tempf];

    end

end

function Particles = Random_regenerate(obj_f, Particles, dim,L,U, n_Random, n_Move)
    if n_Random == 1
        %randomly regenerate the worst particle
        P_random = L + (U - L) * rand;
        Particles(end, :) = [P_random, obj_f(P_random)];
    else
        
        P_random = ones(n_Random, 1) * L + rand(n_Random, dim) .* (U - L);
        f_rand = zeros(n_Random, 1);

        for r = 1:n_Random
            f_rand(r) = obj_f(P_random(r,:));
        end
        Particles(n_Move + 1:end, :) = [P_random,f_rand];
    end
    
end

function Particles = l_s(obj_f, Particles, dim, L,U, n_LS_main, n_LS_leaf)
    LS_radius =  norm(Particles(1, 1:dim) - Particles(randi(2) + 1, 1:dim));

    for i = 1:n_LS_main
        P_LS_main = Particles(1, 1:dim) + LS_radius * (rand(1, dim) - 0.5);
        P_LS_main(P_LS_main > U) = U(P_LS_main > U);
        P_LS_main(P_LS_main < L) = L(P_LS_main < L);
        f_main = obj_f(P_LS_main);

        if Particles(1, end) > f_main
            Particles(2, :) = Particles(1, :);
            Particles(1, :) = [P_LS_main f_main];
        elseif Particles(2, end) > f_main
            Particles(2, :) = [P_LS_main f_main];
        end

        for k = 1:n_LS_leaf
            P_LS_leaf = P_LS_main + LS_radius * (rand(1, dim) - 0.5);
            P_LS_leaf(P_LS_leaf > U) = U(P_LS_leaf > U);
            P_LS_leaf(P_LS_leaf < L) = L(P_LS_leaf < L);
            f_leaf = obj_f(P_LS_leaf);

            if Particles(1, end) >= f_leaf
                Particles(2, :) = Particles(1, :);
                Particles(1, :) = [P_LS_leaf, f_leaf];
                break
            elseif Particles(2, end) > f_leaf
                Particles(2, :) = [P_LS_leaf, f_leaf];
            end
 
        end
    end
end