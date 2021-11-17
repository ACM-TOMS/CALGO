% MAIN.m
%{
/********************              MAIN.m            **********************
 *                                                                        *
 *                                                                        *
 *                        SEQUENTIAL TALBOT SUITE DE                      *
 *                                                                        *
 *                      DRIVER PROGRAM FOR EXAMPLE 4a                     *
 *                          ACCURACY AND TIME TEST                        *
 *                                                                        *
 *                                                                        *
 *                       AUTHOR: Mariarosaria Rizzardi                    *
 *                                                                        *
 *                  mariarosaria.rizzardi@uniparthenope.it                *
 *                                                                        *
 *                  DiST - Dept. of Science and Technology                *
 *                  "Parthenope" University, Naples (Italy)               *
 *                                                                        *
 **************************************************************************
 *                                                                        *
 * REFERENCES                                                             *
 * ==========                                                             *
 * M. RIZZARDI: "Algorithm xxx: TALBOT SUITE DE: APPLICATION OF MODIFIED  *
 *                              TALBOT'S METHOD TO SOLVE DIFFERENTIAL     *
 *                              PROBLEMS".                                *
 *                                 ACM TRANS. MATH. SOFTWARE, VOL. xx,    *
 *                                 NO. x, month year, (PP ##).            *
 *                                                                        *
 **************************************************************************/
%}

%% MENUS
XTmenu = 3;
% XTmenu = menu('???','x in [0, 1],   t in [  1,  5]', ... % 1
%                     'x in [0,2*pi], t in [  1,  5]', ... % 2
%                     'x in [0, 1],   t in [100,500]', ... % 3
%                     'x in [0,2*pi], t in [100,500]');    % 4
switch XTmenu
    case 1 % x in [0, 1], t in [  1,  5]
        Xmin=0;    Xmax=1;
        Tmin=1;    Tmax=5;
    case 2 % x in [0,2*pi], t in [  1,  5]
        Xmin=0;    Xmax=2*pi;
        Tmin=1;    Tmax=5;
    case 3 % x in [0, 1], t in [100,500]
        Xmin=0;    Xmax=1;
        Tmin=100;  Tmax=500;
    case 4 % x in [0,2*pi], t in [100,500]
        Xmin=0;     Xmax=2*pi;
        Tmin=100;   Tmax=500;
end

TOLmenu = 4;
%TOLmenu = menu('tol = ?','1e-6','1e-8','1e-10','1e-12');
switch TOLmenu
    case 1, tol = 1e-6;
    case 2, tol = 1e-8;
    case 3, tol = 1e-10;
    case 4, tol = 1e-12;
end

%%
problem_string
fprintf('\tTmin = %g,  Tmax = %g,  Xmin = %g,  Xmax = %g,  tol = %e\n', Tmin,Tmax,Xmin,Xmax,tol)
fprintf('====================================================================================\n\n')

%%
if ~exist('choice','var')
    choice = menu('?','accuracy','times');
end

switch choice
        
    %%%%%%%%%%%%%%%%%%%%%%
    case 1 % ACCURACY
    %%%%%%%%%%%%%%%%%%%%%%
        if isempty( dir('rel_err.mex*') )
            error('No mex executable file. Run before:  runme  for accuracy test.')
        end
        NXval=9;  NTval=5;
        fprintf('\nNTval = %d;  Tmin = %g;  Tmax = %g;  NXval = %d;  Xmin = %g;  Xmax = %g;  tol = %e', NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);
        % COMPUTE ERRORS
        RELERR1 = rel_err (1,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);
        RELERR2 = rel_err (2,NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol);

        % DISPLAY ERROR MATRICES
        fprintf('\n====================================================================================\n')
        fprintf('\n\t   Ex. 4a: output from ./1SEQ/LTS3_mex/MAIN.m\n')
        fprintf('\n\tLT samples by solving ODE problems by means of MATLAB ode45.m/ode113.m\n')
        fprintf('\n\t%d t in [%g,%g],    %d x in [%g,%g],    tol = %e\n', NTval,Tmin,Tmax,NXval,Xmin,Xmax,tol)
        fprintf('\n====================================================================================\n\n')
        display_errors

    %%%%%%%%%%%%%%%%%%%%%%
    case 2 % TIMES:  time(h,k) for NTval(h) NVxal(k)
    %%%%%%%%%%%%%%%%%%%%%%
        if isempty( dir('SEQ_PDE_skill_TIMES.mex*') )
            error('No mex executable file. Run before:  runme  for time test.')
        end
        NXval=[5 20 120]'; NTval=[5 20 120]';
        PARtime1 = nan*ones(numel(NTval),numel(NXval)); % SEQ_TalbotSUM1
        LTStime1 = nan*ones(numel(NTval),numel(NXval));
        SUMtime1 = nan*ones(numel(NTval),numel(NXval));
        TOTtime1 = nan*ones(numel(NTval),numel(NXval));

        PARtime2 = nan*ones(numel(NTval),numel(NXval)); % SEQ_TalbotSUM2
        LTStime2 = nan*ones(numel(NTval),numel(NXval));
        SUMtime2 = nan*ones(numel(NTval),numel(NXval));
        TOTtime2 = nan*ones(numel(NTval),numel(NXval));

        for i=1:numel(NTval)
            for j=1:numel(NXval)
                % chiama 2 volte Talbot1 per ammortizzare l'inizializzazione dei cicli for
                [PARtime1(i,j),LTStime1(i,j),SUMtime1(i,j),TOTtime1(i,j)] = SEQ_PDE_skill_TIMES (1,NTval(i),Tmin,Tmax,NXval(j),Xmin,Xmax,tol);
                [PARtime1(i,j),LTStime1(i,j),SUMtime1(i,j),TOTtime1(i,j)] = SEQ_PDE_skill_TIMES (1,NTval(i),Tmin,Tmax,NXval(j),Xmin,Xmax,tol);
                fprintf('\n PARTIAL TIME MATRICES: Talbot1\n'); display_times % DISPLAY PARTIAL TIME MATRICES

                [PARtime2(i,j),LTStime2(i,j),SUMtime2(i,j),TOTtime2(i,j)] = SEQ_PDE_skill_TIMES (2,NTval(i),Tmin,Tmax,NXval(j),Xmin,Xmax,tol);
                fprintf('\n PARTIAL TIME MATRICES: Talbot2\n'); display_times % DISPLAY PARTIAL TIME MATRICES

            end
        end

        %% DISPLAY TIME MATRICES
        fprintf('\n====================================================================================\n')
        fprintf('\n\t   Ex. 4a: output from ./1SEQ/LTS3_mex/MAIN.m\n')
        fprintf('\n\tLT samples by solving ODE problems by means of MATLAB ode45.m/ode113.m\n')
        fprintf('\n\tt in [%g,%g],    x in [%g,%g],    tol = %e\n', Tmin,Tmax,Xmin,Xmax,tol)
        fprintf('\n====================================================================================\n\n')
        fprintf('\n FINAL TIME MATRICES\n'); display_times
        fprintf('\n FINAL TIME PERCENTAGE MATRICES\n'); display_time_perc
end

