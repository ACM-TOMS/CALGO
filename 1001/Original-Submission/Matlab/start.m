%% start
% Main file of framework |IPscatt| in MATLAB, i.e. 
% *Computational Framework in MATLAB for the Inverse Medium Problem in Scattering*.
%
% *Warning*: |start| clears variables and stores figures and files in the
% folder |output|. The routine |eval| or |feval| is used in |setInput| and
% |setContrast|.
%
%
%% Quick start
%
% *Adviced usage of MATLAB*
%
% The toolbox |IPscatt| generates several figures and data, which are saved in the folder |output|. 
% We recommend to start MATLAB with the option nodisplay to avoid annoying pop-ups of figures. 
%
%   matlab -nodisplay
%
% *Simple Example*
%
% To run the computational framework, i.e. solve the direct scattering
% problem, add some noise to the data and start the reconstruction process 
% you can just type in
%
%   start
%
% The framework |IPscatt| will set some default parameters, do 
% computation, and save the results in the folder |output|.
%
% *Change parameters*
%
% To compute with other parameters you can set them in a file, e.g. |example|,
% in the folder |inseti| and refer to this file by Parameter |inseti|, e.g.
%
%   inseti = 'example';
%   start
%
% *Attention in usage of the framework*
%
% * |start| saves figures and files in the folder |output| without demand.
%   (Adapt |setInput| and set |out| to 0 or 1 instead of 2 to change this
% behaviour).
% * |start| calls |setInput|, in which *pre-defined variables* (except |closed|, |test|, and |inseti|) *are deleted*.
% * Therefore you *must define* the struct *|seti| in* such a *file* (and not in terminal).
% * The name of a *new file in |inseti|* must differ from existing functions
% (because the file in |inseti| is called by |eval|).
% * The code was written and tested with MATLAB in version *R2016a*.
% 
%
%% Input "Arguments"
% There are no classical input arguments because start is not a function.
%
% From default deviating parameters in struct |seti| are set as described
% in the section Quick Start.
%
% Some input Parameters are described in <example.html>.
%
%
%% Output "Arguments"
% There are no classical output arguments because start is not a function.
%
% *Output results in folder*
%
% Output results are saved in a folder created after |start| like this
%
%   output/20160907T111711_example_test
%
% Structure of foldername: |seti.dirOutput| + / + |seti.dirDatetime| + |_| + |seti.inseti| + |seti.dirSuffix|
% 
% * seti.dirOutput is set to "output"
% * seti.dirDatetime is date and time separated by the character "T"
% * seti.inseti is parameter |inseti| set before |start| (otherwise
% seti.inseti is set to "noinseti")
% * seti.dirSuffix can set in the corresponding file in |inseti|
%
% *Content of output folder*
%
% * figures with prefix |fig_|, see table in the next section (possible
% output formats are png, eps and fig, see <example.html>)
% * |diary.log|: contains terminal output.
% * |inseti_...|: contains a copy of the parameters set in the file
% committed by |inseti|.
% * save_Fmeas.mat: contains exact data $\mathcal{F}_\mathrm{meas}$ and
% noisy simulated data $\mathcal{F}_\mathrm{meas}^\delta$ as 
% |FmeasExact| and |FmeasDelta|.
% * save_qROIcomp_iOutStop.mat: contains the reconstructed contrast
% |qROIcomp| and stop index |iOutStop|.
% * |save_dis.mat|, |save_err.mat|, and |save_dif.mat| are saved if
% |seti.savedata| equals 1. They contain relative discrepancy, error, and difference of 
% computed contrast to its predecessor.
%
%
%% Output: Figures
% A table of possible plotted and saved figures and their number.
%
% <html>
% <table>
% <tr><td><em>Fig. no.</em></td><td><em>Content</em></td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td></td><td><strong>01-10 Preliminary definitions...</strong></td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td></td></td><td><code>expSetup.m</code></td></tr>
% <tr><td>01</td><td>source and measurement points plot...</td></tr>
% <tr><td></td><td><code>setContrast.m</code></td></tr>
% <tr><td>02</td><td>predefined contrast (real part)</td></tr>
% <tr><td>03</td><td>predefined contrast (imag part)</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td>04-10</td><td>currently unused</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td></td><td><strong>11-20 Reconstruction</strong></td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td>11</td><td>Subplots (<code>subplots.m</code>)<br>
%   (1,1) ground truth contrast (real part)<br>
%   (1,2) reconstructed contrast (real part)<br>
%   (2,1) relative error and relative discrepancy<br>
%   (2,2) value of Tikhonov function (function to minimize)</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td></td><td><code>plotAndSaveFigures.m</code></td></tr>
% <tr><td>12</td><td>discrepancy and error (as in subplot)</td></tr>
% <tr><td>13</td><td>Tikhonov functional (as in subplots): M1, M2, MT</td></tr>
% <tr><td>14</td><td>reconstructed contrast (real part)</td></tr>
% <tr><td>15</td><td>reconstructed contrast (imag part)</td></tr>
% <tr><td>16</td><td>difference of reconstructed and true contrast (abs) (switched off by <code>if 0</code>)</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td>17-20</td><td>currently unused</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td></td><td><strong>21-30 Additional for 3D reconstruction (<code>plotAndSaveFigures.m</code>)</strong><br> (contrast in a sectional plane through the scatterer)</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td>21-23</td><td>sectional plane through the predefined contrast (1st, 2nd, 3rd direction)</td></tr>
% <tr><td>24-26</td><td>sectional plane through the reconstructed contrast (1st, 2nd, 3rd direction)</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td>27-30</td><td>currently unused</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td></td><td><strong>31-40 PDA specific additional reconstruction plots</strong><br> over inner iterations (<code>pdaPlot.m</code>)</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td>31</td><td>Parts of minimization functional: fd, fs, fg, fp</td></tr>
% <tr><td>32</td><td>Error in pda and relLinDisInPda</td></tr>
% <tr><td>33</td><td>pda lin. disc. and nonlin. disc.</td></tr>
% <tr><td>34</td><td>parts of min. func. M1 = F(Kh), M2 = G(h) (switched off by <code>if 0</code>)</td></tr>
% <tr><td>35</td><td>ThetaiOut (inner tolerance principle) (switched off by <code>if 0</code>)</td></tr>
% <tr><td>36</td><td>tau, sigma values (switched off by <code>if 0</code>, because currently sigma and tau are not saved)</td></tr>
% <tr><td colspan=2></td></tr>
% <tr><td>37-40</td><td>currently unused</td></tr>
% </table>
% </html>
%
% *For Output figures see also:*
%
% * <pdaPlot.html>
% * <plotAndSaveFigures.html>
% * <plotExpSetup.html>
% * <plotPredefinedContrast.html>
% * <subplots.html>
%
%% Advanced starts
%
% It is possible to run the computational framework with various input
% parameters. Three cases are prepared (others can add similarly).
%
% * |varalpha.m|: various input parameters for the regularization parameter
% $\alpha$ (sparsity).
% * |varbeta.m|: various input parameters for the regularization parameter
% $\beta$ (total variation).
% * |vardelta.m|: various input parameters for the relative noise level $\delta$.
% * |vartol.m|: various input parameters for the tolerance in GMRES.
%
% *Usage of advanced starts*
%
% We demontrate it in case of various $\alpha$ inputs.
%
% In file varalpha:
% 
%   alpha = [100; 500; 1000];
%
% In terminal:
%
%   inseti = 'example'; % setting of seti.alpha inside is overwrittten by varalpha
%   varalpha
%
% A folder in |output| with suffix |varalpha| is created. Reconstruction
% results for $\alpha$ = 100, 500, and 1000 are stored inside.
%
%

%% Code: Part 1: Initialization, set input and start diary
%
% * |init|: <init.html>: Addpath and set parameter |closed| in |init.m|, which differs between
% closed and public version of code (new developments are not published yet).
% * |setInput|: <setInput.html>: Creates folder for output; clears
% variables except of |closed|, |test|, and |inseti|; deals with various input parameters (|varalpha|, |varbeta|, |vardelta|).
% * start diary log (output of terminal is stored in diary.log)
%

disp(' ')
disp('################ -------- start of IPscatt ------- ################')
disp(' ');
disp('               Computational Framework in MATLAB for               ')
disp('              the Inverse Medium Problem in Scattering             ')
disp(' ')

% Initialization
disp(' ')
disp('## init -- initialization -----------------------------------------')
disp(' ')
init;

% Set input and make directories for output
disp(' ')
disp('## setInput -- set input and make directories for output ----------')
disp(' ')
setInput;

% Diary log: start
disp(' ')
disp('## diary -- start diary log ---------------------------------------')
disp(' ')
diary(sprintf('%s/diary%s.log',seti.dirname,seti.fileSuffix)); % store all output in diary
diary on; % later diary off

% Set output depth out
dispDepth = 4; % 0: nearly no display messages; 4: default; 5: too much information
out = 2; % 0: no figures, 1: figures, 2: save figures or files

%% Code: Part 2: Set data
% 
% |setData.m| (<setData.html>) does the following:
%
% * Load real-world data (if set).
% * Set geometry and simulation (grid, kernel, experimental set-up,
% predefined contrast, general settings for figures)
% * Compute exact and noisy data.

% Set data
disp(' ')
disp('## setData -- geometry and experimental set-up, noisy data -------')
disp(' ')
seti = setData(seti,dispDepth,out); % geometry and simulation and fresnel etc. is inside setData...

%% Code: Part 3: Variational reconstruction
%
% * |setRecon.m| (<setRecon.html>): define functions for reconstruction.
% * |recon.m| (<recon.html>): reconstruction process: minimize the 
% Tikhonov functional of the linearized problem by primal-dual algorithm
% and update the contrast $q$ (in public version).
% (The functionals are given in the section "More About...").
%

% Settings for variational reconstruction
disp(' ')
disp('## setRecon -- variational reconstruction -------------------------')
disp(' ')
seti = setRecon(seti,dispDepth);

% Variational reconstruction (process)
disp(' ')
disp('## recon -- variational reconstruction -------------------------')
disp(' ')
seti = recon(seti,dispDepth,out);

%% Code: Part 4: End
%
% * Define end of diray log.
% * Save workspace (if seti.saveworkspace was set to 1)

% Diary log: end
disp(' ')
disp('## diary -- end diary log -----------------------------------------')
disp(' ')
diary off;

% Save workspace
seti = checkfield(seti,'saveworkspace',0);
if seti.saveworkspace == 1 && out >= 2
    disp(' ')
    disp('## save workspace ---------------------------------------------')
    disp(' ')
    save(sprintf('%s/workspace.mat',seti.dirname));
end

%% More About... Theory: The direct and inverse scattering from inhomogeneous media
%
% We give a brief description of the direct and inverse scattering problem.
% A reference is e.g. Chapter 8 in [1].
%
% *The direct scattering problem*
%
% We consider a time-harmonic *incident wave* $u^\mathrm{i}: \bf{R}^d \to \bf{C}$ in dimension $d = 2, 3$ 
% with time-dependence $\exp(-\mathrm{i} \omega t)$, where $\omega > 0$ is the angular frequency. 
% 
% When we omit the time-dependence, the incident field satisfies the
% *Helmholtz equation* 
% 
% $$ \Delta u^\mathrm{i}(x) + k^2 u^\mathrm{i}(x) = 0, \quad x \in \bf{R}^d $$
% 
% with constant wave number $k>0$.
%
% The scattering object is described by a refractive index function
% $n: \bf{R}^d \to \bf{C}$.
% This function equals $1$ outside a bounded domain.
% 
% When the incident wave interacts with a scattering object 
% a *scattered wave* $u^\mathrm{s}$ is generated. 
% The *total wave*
%
% $$ (1) \quad u^\mathrm{t} := u^\mathrm{i} + u^\mathrm{s} $$
%
% admits the *Helmholtz equation*
%
% $$ (2) \quad \Delta u^\mathrm{t} + k^2 n^2 u^\mathrm{t} = 0 \quad \textrm{in } \bf{R}^d $$.
%
% The scattered wave $u^\mathrm{s}$ satisfies 
% *Sommerfeld's radiation condition*
% 
% $$ (3) \quad \lim\limits_{|x|\to\infty} |x|^{(d-1)/2}\left(
% \frac{\partial}{\partial |x|} - \mathrm{i} k\right) u^\mathrm{s}(x) = 0, $$
%
% uniformly in all directions 
% $\hat x = x/|x| \in \bf{S} := \{ y \in \bf{R}^d: \, |y| =1\}$, 
% where $|y|:= \sqrt{y_1^2 + \ldots + y_d^2}$.
%
% The *direct scattering problem* is to find a solution $u$ to equations (1)-(3).
% 
% We will not solve the above PDE. Instead we will solve the
% "Lippmann-Schwinger equation", which is a equivalent reformulation.
%
% Further, instead of the refractive index we will consider the contrast
%
% $$ q := n^2-1 \quad \textrm{in } \bf{R}^d. $$
%
% *The inverse scattering problem*
%
% Keep in mind that $u^\mathrm{s} = u^\mathrm{t}-u^\mathrm{i}$. 
% If we know the incident wave $u^\mathrm{i}$ and measure the total wave $u^\mathrm{t}$, 
% we can compute the scattered wave $u^\mathrm{s}$.
% 
% Then the scattered wave is used as data for the *inverse scattering problem*: 
% Reconstruct the contrast function _q_ from several measurements of the scattered wave.
%
% For this we will need *multi-static measurements*, i.e. several
% transmitters propagate incident waves one after another and the generated
% scattered fields are measured by several receivers. (In real-world total 
% fields are measured and scattered fields are computed subtracting the 
% corresponding incident fields.)
%
% *More theoretical background*
%
% The exact techniques how to tackle the direct and inverse scattering from 
% inhomogeneous media are described in [2]. In the next section we have a
% brief look at the function we minimize for reconstruction.
%
% *Function to minimize...*
%
% *... We would like to* find a contrast $q$, that minimizes the Tikhonov
% functional of the non-linearized problem, i.e.:
%
% $$ \min_{q \in X} 
%      \underbrace{\frac{1}{2}\|\mathcal{F}(q)-F_\mathrm{meas}^\delta\|_\mathrm{F}^2}_{\mathrm{discrepancy}}
%    + \underbrace{\alpha \|q\|_\mathrm{1}}_{\mathrm{sparsity\ penalty}}
%    + \underbrace{\beta \| \nabla q \|_\mathrm{1}}_{\mathrm{total\ variation\ penalty}}
%    + \underbrace{
%          \delta_{[a,b]}( \mathrm{Re}(q) ) +
%          \delta_{[c,d]}(\mathrm{Im}(q) )}_{\mathrm{penalty\ for\ physical\ bounds}
%          }.
%    $$
%
% * $\|\cdot\|_\mathrm{F}$ is the Frobenius norm.
% * Actually all norms are weighted. (We ommit it to keep notation simple.)
% * $\delta_{[a,b]}(x)$ is a indicator function, i.e. 
%   is 0, if all entries of the vector x are between $a$ and $b$,
%   and is $\infty$, otherwise.
%
% *... What we do* is to minimize the Tikhonov functional of the linearized
% problem, i.e. (this is done in <pda.html>)
%
% $$ \min_{h \in X} 
%      \underbrace{\frac{1}{2}\|\mathcal{F}'(q)[h]+\mathcal{F}(q) - F_\mathrm{meas}^\delta\|_\mathrm{F}^2}_{\mathrm{discrepancy\ (linearized\ problem)}}
%    + \underbrace{\alpha \|q+h\|_\mathrm{1}}_{\mathrm{sparsity\ penalty}}
%    + \underbrace{\beta \| \nabla (q+h) \|_\mathrm{1}}_{\mathrm{total\ variation\ penalty}}
%    + \underbrace{
%          \delta_{[a,b]}( \mathrm{Re}(q+h) ) +
%          \delta_{[c,d]}(\mathrm{Im}(q+h) )}_{\mathrm{penalty\ for\ physical\ bounds}
%          }
%    $$
%
% and then update $q := q + h$ (this is done in <minPda.html>).
%
% Next, the Tikhonov functional of the linearized problem
% is minimized again. (This minimization is done by application of primal-dual algorithm,
% see [3], on this problem, see [2].)
%
% The whole process minimizes the Tikhonov functional of the non-linearized
% inversion problem and leads to the seeked contrast $q$.
%
%
%% Authors
% The authors of |IPscatt| are:
%
% * Florian B&uuml;rgel (Center for Industrial Mathematics, University of Bremen, Germany, fbuergel@uni-bremen.de.)
% * Kamil S. Kazimierski (Institute of Mathematics and Scientific Computing, University of Graz, Austria, kazimier@uni-graz.at. The Institute of Mathematics and Scientific Computing is a member of NAWI Graz (www.nawigraz.at) and BioTechMed Graz (www.biotechmed.at).)
% * Armin Lechleiter (Center for Industrial Mathematics, University of Bremen, Germany.)
%
% Note that the source code of <gmresKelley.html> in the folder |3rdparty| was not written by the authors and contains the corresponding license information.
%
%% References
%
% * [1] David Colton and Rainer Kress. _Inverse Acoustic and Electromagnetic Scattering Theory_. Springer, New York, 2013.
% * [2] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
% * [3] Antonin Chambolle and Thomas Pock. A first-order primal-dual algorithm for convex problems with applications to imaging. _Journal of Mathematical Imaging and Vision_, 40(1):120-145, 2011.
%
%% See Also
%
% * <example.html>
%
% * <varalpha.html>
% * <varbeta.html>
% * <vardelta.html>
% * <vartol.html>
%
% * <init.html>
% * <setInput.html>
% * <setData.html>
% * <setRecon.html>
% * <recon.html>
