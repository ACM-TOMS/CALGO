Revision of SFSDP111 from SFSDP100
                                                        July 2009

As a default, this new versionof SFSDP calls sedumiwrap, which is a 
MATLAB version of SDPA, to solve an SDP relaxation problem. If users would 
ike to use SeDuMi instead of SDPA, they need to modify the MATLAB program 
SFSDP.m by replacing the line 

SDPsolverDefault = 'sdpa'; 

by 

SDPsolverDefault = 'sedumi'; 

Or users can set a new parameter pars.SDPsolver = 'sedumi' when they 
execute SFSDPplus or SFSDP. 
 
