%% pdaPlot
% Plots and saves specifically for primal-dual algorithm.
%
% Reserved numbers: 31-40. Used numbers 31-36, see <start.html>.
% 
%% Syntax
%
%   pdaPlot(iOut,seti,pdas,out)
%
%% Description
% |pdaPlot(iOut,seti,pdas)| plots and saves primal-dual algorithm specific
% figures of outer iteration |iOut|.

%% Input Arguments
%
% * iOut    : number of outer iteration
% * seti    : structural array
% * pdas    : structural array with specific results from primal-dual
%             algorithm, see <recon.html> for details.
% * out     : Output depth: generate plots (1), generate plots and save them (2).
%
% Specific fields in |seti| influencing the figures are described in 
% <setGeomSim.html> in the section "Subfunction: setFigureSettings".
%
%% Output Arguments
%
% Figures are plotted and saved, see <start.html>
%
%% See Also
%
% * <start.html>
% * <savePngFig.html>
% * <setGeomSim.html>
% * <recon.html>
%
%% Code

function pdaPlot(iOut,seti,pdas,out)

ps = pdas.pdaStopInd;
x = 1:ps;
clear ps;

%%
% *figure 31: terms fd, fs, fg, fp*
figure(31);
set(gcf,'Visible',seti.plotVisible);
[ax,h1,h2] = plotyy([x',x',x'],[pdas.minf.fs(x),pdas.minf.fg(x),pdas.minf.fp(x)],x',pdas.minf.fd(x));
legend([h1;h2],'fs (left axis)','fg (left axis)','fp (left axis)','fd (right axis)');
if ~seti.plotPublish
    xlabel('pda (inner) iterations');
    title('Parts of min. functional');
end
axis(ax,'square');
if out == 2
    savePngFig(31,iOut,seti);
end

%%
% *figure 32: pda err and relLinDisInPda*
figure(32);
set(gcf,'Visible',seti.plotVisible);
[ax,h1,h2] = plotyy(x,pdas.errInPda(x)',x,pdas.relLinDisInPda(x)');
legend([h1;h2],'rel. error','relLinDis = disLinInPda/dis');
if ~seti.plotPublish
    xlabel('pda iterations');
    if pdas.ThetaiOutV(iOut) ~= 0
        fig11title = sprintf('Pda inner iterations: rel. error, rel. lin. dis (ThetaiOut = %1.2g)',pdas.ThetaiOutV(iOut));
    else
        fig11title = sprintf('Pda inner iterations: rel. error, rel. lin. dis');
    end
    title(fig11title);
end
axis(ax,'square');
if out == 2
    savePngFig(32,iOut,seti);
end

%%
% *figure 33: pda lin. and nonlin. disc.*

% later: maybe only interesting if ibreak == true...
figure(33);
set(gcf,'Visible',seti.plotVisible);
[ax,h1,h2] = plotyy(1:iOut,pdas.relDis(1:iOut)',1:iOut,pdas.pdaNv(1:iOut)');
legend([h1;h2],'rel dis = lin. dis / nonlin. dis','max. inner iterations');
if ~seti.plotPublish
    xlabel('outer iterations');
    title('Pda outer iterations: lin. and non-lin. dis');
end
axis(ax,'square');
if out == 2
    savePngFig(33,iOut,seti);
end

%%
% *figure 34: pda inner iterations: parts of min. func. M1 = F(Kh), M2 = G(h)*

if 0
    figure(34);
    set(gcf,'Visible',seti.plotVisible);
    %plot(x,pdas.MTvN(x),x,pdas.M1vN(x),x,pdas.M2vN(x));
    [ax,h1,h2] = plotyy(x',pdas.M1vN(x),[x',x'],[pdas.M2vN(x),pdas.MTvN(x)]);
    if ~seti.plotPublish
        legend([h1;h2],'F(Kh)','G(h)','F(Kh)+G(h)');
        xlabel('pda iterations');
        title('Pda inner iterations: F(Kh) and G(h)');
    end
    axis(ax,'square');
    if out == 2
        savePngFig(34,iOut,seti);
    end
end

%%
% *figure 35: ThetaiOut (inner tolerance principle)*

if 0
    figure(35);
    set(gcf,'Visible',seti.plotVisible);
    plot(1:iOut,pda.ThetaiOutV(1:iOut))
    if ~seti.plotPublish
        xlabel('outer iterations');
        title('ThetaiOut (inner tolerance principle)');
    end
    axis square;
    if out == 2
        savePngFig(35,iOut,seti);
    end
end

%%
% *figure 36: tau, sigma values plot*

if 0
    % Currently tauVal and sigmaVal are not stored.
    figure(36);
    set(gcf,'Visible',seti.plotVisible);
    plot(x,pda.tauVal(x),x,pda.sigmaVal(x))
    legend('tau','sigma');
    if ~seti.plotPublish
        xlabel('pda (inner) iterations');
        title('Values of sepsizes tau and sigma');
    end
    axis square;
    if out == 2
        savePngFig(36,iOut,seti);
    end
end

end
