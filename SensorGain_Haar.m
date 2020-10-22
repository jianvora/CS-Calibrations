%% PARAMETERS

clear; clc; close all;
N = 128; 
r = 0.2;
desirednormx = 100;
tolerance = 0.01; % difference in consecutive estimates of x
st=0;
noisefrac = 0.02;
maxiter = 25;
numdeltas = 1;
r_gain = 0.2;

%numruns = 5;
numruns = 1;
numstarts = 5;

%Mvals = [40,50,60,70,80,90,100,110,120];
Mvals = 90;
svals = 20;
%svals = [10,20,30,40,50,60,70,80];

Mi = N;
H  = haar(N);
i_H = H';


 
avgbigxerrormat = zeros(numel(Mvals), numel(svals));
avgbigdeltaerrormat = zeros(numel(Mvals), numel(svals));
avgbiggainerrormat = zeros(numel(Mvals),numel(svals));

bigxerrormat = zeros(numruns, numel(Mvals), numel(svals));
bigdeltaerrormat = zeros(numruns, numel(Mvals), numel(svals));
biggainerrormat = zeros(numel(Mvals),numel(svals));


tic

for run = 1:numruns
    Mindex = 0;
    for M = Mvals
      %  numdeltas = M/15;
        Mindex=Mindex+1;
        fprintf('\nMindex=%d\n', Mindex);
        sindex = 0;
        for sparse = svals
            sindex = sindex+1;
            fprintf('Run=%d. Mindex=%d. Sindex=%d\n', run, Mindex, sindex);

            %% GENERATE DATA

            %sparse
            x = zeros(N,1);
            x(randsample(N,sparse)) = normrnd(0,5,sparse,1);
            x = x * desirednormx / norm(x,2);

            U = -(Mi-1)/2:(Mi-1)/2;
            inds = sort(randperm(Mi,M));
            U = U(inds);
            
            U = U';
            deltasmall = r*2*(rand(numdeltas,1) - 0.5);
            deltaU = zeros(M,1);
            for i = 1:numdeltas
                segment = floor(M/numdeltas);
                indices = ((i-1)*segment+1):i*segment;
                deltaU(indices) = deltasmall(i);
            end
            deltaU(numdeltas*floor(M/numdeltas):M) = deltasmall(numdeltas);

            Utilde = U + deltaU;
            Delta = diag(deltaU);
            gain = 1-r_gain + 2*r_gain*rand(M,1);
            gainmat=diag(gain);

            ntouse = -(N)/2+1:(N)/2;
            %Assemble F on -(N-1)/2:(N-1)/2;
            k = -(N)/2+1:(N)/2;
            F = zeros(M,N);
            Ftrue = zeros(M,N);
            for row = 1:M
                for col = 1:N
                    F(row,col) = exp(-2j * pi * U(row) * k(col)/N) * (1/sqrt(N));
                    Ftrue(row,col) = exp(-2j * pi * Utilde(row) * ntouse(col)/N) * (1/sqrt(N));
                end
            end

            X = diag(-2j * pi * (-(N-1)/2:(N-1)/2)/N) * (1/sqrt(N)); %derivative of the shifted Fourier
            
            
           
            Ftrue = Ftrue*H;
            F = F*H;
            yExact = gainmat*Ftrue * x;
            
            sdnoise = noisefrac*norm(yExact,1)/N;
            noise = normrnd(0,sdnoise,M,1) + 1j*normrnd(0,sdnoise,M,1);
                      
            
            
            y_measured = yExact + noise;

            
            
            %% RECOVERY
            
           objvector = zeros(numstarts,1);
           simpleobjvector = zeros(numstarts,1);
           actualerrorvec = zeros(numstarts,1);
           
           best_x_estimated = zeros(N,1);
           best_gain_estimated = zeros(M,1);
           best_deltasmall_estimated = zeros(numdeltas,1);
           bestobj = Inf;
           
           
           for starti = 1:numstarts
            fprintf('Run=%d. Mindex=%d. Sindex=%d. Start = %d\n', run, Mindex, sindex, starti);   
               
            % Two step optimization
            x_estimated = normrnd(10,5,N,1);
            gain_estimated = 1-r_gain + 2*r_gain*rand(M,1);
            gain_mat_estimated = diag(gain_estimated);
            %beta_estimated = 0.01*(rand(N,1) - 0.5);
            
            deltasmall_estimated = r*2*(rand(numdeltas,1) - 0.5);
%           deltasmall_estimated = deltasmall;
            initial_delta = deltasmall_estimated;
            initial_delta_error = norm(deltasmall - deltasmall_estimated) / norm(deltasmall);  
            beta_estimated = zeros(M,1);
            for i = 1:numdeltas
                segment = floor(M/numdeltas);
                indices = ((i-1)*segment+1):i*segment;
                beta_estimated(indices) = deltasmall_estimated(i);
            end
            beta_estimated(numdeltas*floor(M/numdeltas):M) = deltasmall_estimated(numdeltas);
            
            x_errors = [];
            beta_errors = [];
            
            %% Unconstrained optimization

            preverror_x = Inf;
            preverror_beta = Inf;
            preverror_gain = Inf;

            differror = Inf;
            iterations = 0;
            differrorsvec = [];
            prev_x_estimated = zeros(N,1);
            
            
            Fprime = F;
            
            fprintf('iter ');
            
            while differror > tolerance && iterations < maxiter
            %while iterations < maxiter
                
                fprintf('%d ', iterations+1);

%1. Delta fixed.Gain fixed. Estimate x.
                
                % Set up Fprime
                for row = 1:M
                    for col = 1:N
                        Fprime(row,col) = exp(-2j * pi * (U(row)+beta_estimated(row)) * ntouse(col)/N) * (1/sqrt(N));
                    end
                end
                Fprime = Fprime*H;
                lambda=20;
                cvx_begin quiet
                    variable x_estimated(N) 
                    minimize norm(x_estimated, 1) + lambda*norm(y_measured - (gain_mat_estimated*Fprime * x_estimated), 2)
                cvx_end
                
%2. x fixed. Delta fixed. Gain to be estimated
                for row = 1:M
                    for col = 1:N
                        Fprime(row,col) = exp(-2j * pi * (U(row)+beta_estimated(row)) * ntouse(col)/N) * (1/sqrt(N));
                    end
                end
                Fprime = Fprime*H;
                for i=1:M
                    gainvec = zeros(M,1);
                    besterritergain = Inf;
                    bestguessedgain=1;
                    for guessedgain = (1-r_gain):(r_gain/50):(1+r_gain)
                        gainvec(i)=guessedgain;
                        erritergain = norm(y_measured(i) - guessedgain*Fprime(i,:)*x_estimated,2);
                        if erritergain < besterritergain
                            %fprintf('%d %f\n',i,guessedbeta)
                            besterritergain = erritergain;
                            bestguessedgain = guessedgain;
                        end
                    end
                    gain_estimated(i) = bestguessedgain; 
                    gain_mat_estimated = diag(gain_estimated);
                end

                

%3. x fixed. Gain fixed. Estimate Delta(actually deltaU=beta)

                for i = 1:numdeltas
                    segment = floor(M/numdeltas);
                    ind = ((i-1)*segment+1):i*segment;

                    yiter = y_measured(ind);
                    giter = gain_estimated(ind);
                    gmatiter = diag(giter);
                    besterriter = Inf;
                    bestguessedbeta = 0;
                    
                    for guessedbeta = -r:(r/50):r
                        
                        Fiter = zeros(segment,N);
                        for row = 1:segment
                            for col = 1:N
                                Fiter(row,col) = exp(-2j * pi * (U(ind(row))+guessedbeta) * ntouse(col)/N) * (1/sqrt(N));
                            end
                        end
                        Fprime = Fprime*H;
                        erriter = norm(yiter - gmatiter*Fiter*x_estimated, 2);
                        
                        if erriter < besterriter
                            %fprintf('%d %f\n',i,guessedbeta)
                            besterriter = erriter;
                            bestguessedbeta = guessedbeta;
                        end
                    end
                    
                    deltasmall_estimated(i) = bestguessedbeta;
                    beta_estimated(ind) = deltasmall_estimated(i);
                end
                
                beta_estimated(numdeltas*floor(M/numdeltas):M) = deltasmall_estimated(numdeltas);
                beta_estimated(beta_estimated>r) = r;
                beta_estimated(beta_estimated<-r) = -r;


                iterations = iterations + 1;
         
                differror = norm((prev_x_estimated/(norm(prev_x_estimated)+1e-6)) - (x_estimated/(norm(x_estimated)+1e-6)))
                %differror = norm(prev_x_estimated - x_estimated)/norm(prev_x_estimated)
                prev_x_estimated = x_estimated;
            end

            ['Total ' num2str(iterations) ' iterations']

            %% RESULTS
            x_error = norm((x/norm(x)) - (x_estimated/norm(x_estimated)),2);
            delta_error = norm(deltasmall - deltasmall_estimated) / norm(deltasmall);
            
            for row = 1:M
                for col = 1:N
                    Fprime(row,col) = exp(-2j * pi * (U(row)+beta_estimated(row)) * ntouse(col)/N) * (1/sqrt(N));
                end
            end
            obj = norm(x_estimated, 1) + lambda*norm(y_measured - (gain_mat_estimated*Fprime *H*x_estimated), 2);
            simpleobj = norm(y_measured - (gain_mat_estimated*Fprime *i_H* x_estimated), 2);
                        
            if obj < bestobj
                best_x_estimated = x_estimated;
                best_deltasmall_estimated = deltasmall_estimated;
                best_gain_estimated = gain_estimated;
                bestobj = obj;
                
            end
            objvector(starti) = obj;
            simpleobjvector(starti) = simpleobj;
            actualerrorvec(starti) = x_error;
           end %multistart loop

            bigxerrormat(run, Mindex, sindex) = norm((x/norm(x)) - (best_x_estimated/norm(best_x_estimated)),2);
            bigdeltaerrormat(run, Mindex, sindex) = norm(deltasmall - best_deltasmall_estimated) / norm(deltasmall);
            biggainerrormat(run, Mindex, sindex) = norm(gain - gain_estimated) / norm(gain);
        end
    end
   avgbigxerrormat = avgbigxerrormat + squeeze(bigxerrormat(run,:,:));
   avgbigdeltaerrormat = avgbigdeltaerrormat + squeeze(bigdeltaerrormat(run,:,:)); 
end

toc
avgbigxerrormat = avgbigxerrormat/numruns;
avgbigdeltaerrormat = avgbigdeltaerrormat/numruns;
save('haarwavelet_onrnd.mat');
% 
figure();
avgbigxerrormat(avgbigxerrormat>1)=1;
image(svals, Mvals, repmat(avgbigxerrormat, [1 1 3]));
ylabel('M (number of measurements)');
xlabel('s (number of non-zero entries)');
% %save(['large_beta_xerrors_' num2str(numdeltas) '_bypasstaylor_multistart100deltas.mat'], 'bigxerrormat');
% 
% max_x_error = max(max(avgbigxerrormat))
% min_x_error = min(min(avgbigxerrormat))
% 
% figure();
% avgbigdeltaerrormat(avgbigdeltaerrormat>1)=1;
% image(svals, Mvals, repmat(avgbigdeltaerrormat ,[1 1 3]));
% ylabel('M (number of measurements)');
% xlabel('s (number of non-zero entries)');
% %save(['large_beta_deltaerrors' num2str(numdeltas) '_bypasstaylor_multistart100deltas.mat'], 'bigdeltaerrormat');
% 
% max_delta_error = max(max(avgbigdeltaerrormat))
% min_delta_error = min(min(avgbigdeltaerrormat))

