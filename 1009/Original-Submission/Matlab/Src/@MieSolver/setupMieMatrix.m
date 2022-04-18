%-----------------------------------
% method of the MieSolve class to setup
% the Mie system matrix
%-----------------------------------

% Copyright 2014-2019 Stuart C. Hawkins
% 	
% This file is part of MIESOLVER.
% 
% MIESOLVER is free software: you can redistribute it and/or modify	
% it under the terms of the GNU General Public License as published by	
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% MIESOLVER is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with MIESOLVER.  If not, see <http://www.gnu.org/licenses/>.


function [matrix,indexes,localIndexes,rows,columns,bcrows] = setupMieMatrix(self)

% check that there are some scatterers
if self.num_scatterers==0
    error('This MieSolver has no scatterers.')
end

% work out dimensions of matrix and setup the matrix
% in sparse format
[nm,nb,aa] = self.matrixDimensions();
matrix = sparse(sum(nm),sum(nm));

% initalise indices into the matrix
columns = [];
bcrows = [];
rows = [];

% loop through the scatterers
for scatk=1:self.num_scatterers
    
    % setup wavefunction indexes
    n=-self.nmax(scatk):self.nmax(scatk);
    n=n(:);
    
    % get scatterer for concise code
    scatterer = self.scatterer{scatk};
    
    % get the key to the different variables
    variables=self.get_variables(scatterer);
    
    % the bcs array acts as a key to the different rows in the
    % matrix. For each boundary we have a row corresponding to
    % Dirichlet BC and a row corresponding to Neumann BC.
    % Depending on the BC we wish to implement we can switch off
    % one of these by setting the entry to '0'
    bcs=repmat(['D';'N'],scatterer.getNumRegion-1,1);
    if ~isfinite(scatterer.getRefractiveIndex(1))
        if scatterer.getbc(1) == boundary_condition.soft
            bcs(2)='0';
        elseif scatterer.getbc(1) == boundary_condition.hard
            bcs(1)='0';
        elseif scatterer.getbc(1) == boundary_condition.robin
            bcs(1)='0';
        end
    end
    
    % ii{j} gives indexes to block j of the matrix
    for j=1:2*scatterer.getNumRegion()
        localIndexes{scatk}{j} = (j-1)*nb(scatk)+1:j*nb(scatk);
        indexes{j}=aa{scatk}(localIndexes{scatk}{j});
    end
    
    % copy kwave to shorten code
    kwave=self.kwave;
    
    % loop through the BCs... each BC gives a block row in the
    % matrix
    for i=1:scatterer.getNumRegion-1
                
        % get the refractive index either side of the BC
        inner_m=scatterer.getRefractiveIndex(i);
        outer_m=scatterer.getRefractiveIndex(i+1);
        rad=scatterer.getRadius(i);
        
        % get the density either side of the BC
        inner_rho=scatterer.getDensity(i);
        outer_rho=scatterer.getDensity(i+1);

        % determine the coefficient in the J' and H' matrices
        if scatterer.getbc(i) ~= boundary_condition.transmission

            % Note: we must not scale the matrices J' and H' in the 
            % Robin BC case because the J and H matrices are not scaled.
            
            % then du/dn are not scaled... but pick up m from the
            % derivative of J(m k r)
            inner_alpha = inner_m;
            outer_alpha = outer_m;
            
        else
            
            % the coefficident depends on
            % whether we have TM or TE type polarisation
            if strcmp(self.transmissionType,'TE')
                % TE polarisation du/dn is continuous across the interface...
                % also pick up m from derivative of J(m k r)
                inner_alpha = inner_m;
                outer_alpha = outer_m;
            elseif strcmp(self.transmissionType,'TM')
                % TM polarisation (1/m^2) du/dn is continuous across the interface...
                % also pick up m from derivative of J(m k r)
                inner_alpha = 1/inner_m;
                outer_alpha = 1/outer_m;
            elseif strcmp(self.transmissionType,'acoustic')
                % check density has been set
                if isnan(inner_rho) || isnan(outer_rho)
                    error('Density is needed for acoustic transmission BC but is not set')
                end
                % acoustic (1/rho) du/dn is continuous across the interface...
                % also pick up m from derivative of J(m k r)
                inner_alpha = inner_m/inner_rho;
                outer_alpha = outer_m/outer_rho;
            elseif strcmp(self.transmissionType,'custom')                
                % check transmission function has been set
                if isempty(self.transmissionFun)
                    error('Use transmissionCustom method to set the transmission BC function')
                end
                
                % self.transmissionFun(m) du/dn is continuous across the interface...
                % also pick up m from derivative of J(m k r)
                inner_alpha = inner_m * self.tranmissionFun(inner_m);
                outer_alpha = outer_m * self.tranmissionFun(outer_m);
                
            else
                error('Transmission type %s not recognised',self.transmissionType);
            end
        
        end
        
        % set the matrix entries
        % J  H  | -J  -H
        % JD HD | -JD -HD
        matrix(indexes{2*i-1},indexes{2*i-1}) = diag(besselj(abs(n),inner_m*kwave*rad));
        matrix(indexes{2*i},indexes{2*i-1}) = diag(inner_alpha*kwave*besseljd(abs(n),inner_m*kwave*rad));
        matrix(indexes{2*i-1},indexes{2*i}) = diag(besselh(abs(n),inner_m*kwave*rad));
        matrix(indexes{2*i},indexes{2*i}) = diag(inner_alpha*kwave*besselhd(abs(n),inner_m*kwave*rad));
        matrix(indexes{2*i-1},indexes{2*(i+1)-1}) = -diag(besselj(abs(n),outer_m*kwave*rad));
        matrix(indexes{2*i},indexes{2*(i+1)-1}) = -diag(outer_alpha*kwave*besseljd(abs(n),outer_m*kwave*rad));
        matrix(indexes{2*i-1},indexes{2*(i+1)}) = -diag(besselh(abs(n),outer_m*kwave*rad));
        matrix(indexes{2*i},indexes{2*(i+1)}) = -diag(outer_alpha*kwave*besselhd(abs(n),outer_m*kwave*rad));
        
    end
    
    % deal with Robin case by modifying the JD row to
    % JD + 1i*lambda* J
    if scatterer.getbc(1) == boundary_condition.robin
        matrix(indexes{2},aa{scatk}) = matrix(indexes{2},aa{scatk}) +  1i*scatterer.getRobinParameter()*matrix(indexes{1},aa{scatk});
    end
    
    % work out which columns we need based on the key to the
    % variables
    for k=1:2*scatterer.getNumRegion()
        if variables(k)=='J' || variables(k)=='H'
            columns=[columns,indexes{k}];
        end
    end
    
    % work out which columns correspond to values given by the
    % incident field coefficients
    for k=1:2*scatterer.getNumRegion
        if variables(k)=='I'
            bcrows=[bcrows,indexes{k}];
        end
    end
    
    % work our which rows we need based on the BC key
    for k=1:2*(scatterer.getNumRegion-1)
        if bcs(k)=='D' || bcs(k)=='N'
            rows=[rows,indexes{k}];
        end
    end
    
end
