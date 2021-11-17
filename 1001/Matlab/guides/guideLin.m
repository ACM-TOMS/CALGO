init;
seti.incNb = 1;
seti.measNb = 2;
seti = setData(seti);
q = seti.qROIexact; % define some contrast
h = 0.1.*(rand(size(q)) + 1i*rand(size(q))); % define some random update h

%% Test for sufficent changes to deal with a nother wave number

if 0
    seti.k = 50;

    % Change specific... for wave number
    seti = setKernel(seti);
    seti = setIncField(seti);
    seti = setMeasKer(seti);
    [FFqMeas,~,~] = forward(seti,q);

    % Comparison...
    seti = setData(seti);
    [FFqMeasComp,~,~] = forward(seti,q);
    norm(FFqMeasComp-FFqMeas)/norm(FFqMeas)
end

%% Compute F(q+h) and the linear approximation F'(q)[h]+F(q):

% kVec = [5, 250];
% kVec = 0:10:150;
kVec = 0:1:150;
rel = zeros(size(kVec));

i = 0;
for k = kVec
    i = i+1;
    
    seti.k = k;
    % Specific changes because changed wave number
    % A expensive alternative is: seti = setData(seti);
    seti = setKernel(seti);
    seti = setIncField(seti);
    seti = setMeasKer(seti);

    % Linearization of $\mathcal{F}(q+h)$ at $q$:
    % $\mathcal{F}(q+h)$ is approximated by its linearization $\mathcal{F}'(q)[h]+\mathcal{F}(q)$

    % $\mathcal{F}(q+h):
    [FFqhMeas,~,~] = forward(seti,q+h);

    % Linearization F'(q)[h] + F(q):
    [FFqMeas,~,~] = forward(seti,q); % $|FFqROI| = \mathcal{F}(q)$:
    [JA,JB] = derivative(seti,q);
    FFqhMeasLin = JA*diag(h)*JB + FFqMeas;

    % Relative error:
    rel(i) = norm(FFqhMeas-FFqhMeasLin)/norm(FFqhMeas);
end

%% Relative error in dependence of wave number k?

figure(1); h = plot(kVec,rel); axis square;
set(gca,'FontSize',20); axis square; set(h,'LineWidth',2); print(1,'-depsc','guideLin.eps');
