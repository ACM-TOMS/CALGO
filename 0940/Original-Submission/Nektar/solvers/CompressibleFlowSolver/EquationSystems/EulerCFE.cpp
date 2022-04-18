///////////////////////////////////////////////////////////////////////////////
//
// File EulerCFE.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Euler equations in conservative variables without artificial
// diffusion
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <CompressibleFlowSolver/EquationSystems/EulerCFE.h>

namespace Nektar
{
    string EulerCFE::className = 
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "EulerCFE", EulerCFE::create, 
            "Euler equations in conservative variables without "
            "artificial diffusion.");
    
    EulerCFE::EulerCFE(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : CompressibleFlowSystem(pSession)
    {
    }

    void EulerCFE::v_InitObject()
    {
        CompressibleFlowSystem::v_InitObject();

        if(m_session->DefinesSolverInfo("PROBLEMTYPE"))
        {
            std::string ProblemTypeStr = m_session->GetSolverInfo("PROBLEMTYPE");
            int i;
            for(i = 0; i < (int) SIZE_ProblemType; ++i)
            {
                if(NoCaseStringCompare(ProblemTypeMap[i],ProblemTypeStr) == 0)
                {
                    m_problemType = (ProblemType)i;
                    break;
                }
            }
        }
        else
        {
            m_problemType = (ProblemType)0;
        }

        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs     (&EulerCFE::DoOdeRhs,        this);
            m_ode.DefineProjection (&EulerCFE::DoOdeProjection, this);
        }
        else
        {
            ASSERTL0(false, "Implicit CFE not set up.");
        }
    }
    
    EulerCFE::~EulerCFE()
    {

    }

    void EulerCFE::v_PrintSummary(std::ostream &out)
    {
        CompressibleFlowSystem::v_PrintSummary(out);
        out << "\tProblem Type    : " << ProblemTypeMap[m_problemType] << endl;
    }

    void EulerCFE::v_SetInitialConditions(NekDouble initialtime, bool dumpInitialConditions)
    {
        switch(m_problemType)
        {
        case eIsentropicVortex:
        {
            SetIsenTropicVortex(initialtime);
        }
        break;
        default:
        {
            EquationSystem::v_SetInitialConditions(initialtime,false);
        }
        break;
        }

        if(dumpInitialConditions)
        {
            // dump initial conditions to file
            std::string outname = m_sessionName + "_initial.chk";
            WriteFld(outname);
        }
    }

    void EulerCFE::DoOdeRhs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int i;
        int nvariables = inarray.num_elements();
        int npoints    = GetNpoints();

        Array<OneD, Array<OneD, NekDouble> > advVel;
        
        m_advection->Advect(nvariables, m_fields, advVel, inarray, outarray);

        for (i = 0; i < nvariables; ++i)
        {
            Vmath::Neg(npoints,outarray[i],1);
        }
    }
    
    void EulerCFE::DoOdeProjection(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
            Array<OneD,       Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        int i;
        int nvariables = inarray.num_elements();

        switch(m_projectionType)
        {
        case MultiRegions::eDiscontinuous:
        {
            // Just copy over array
            int npoints = GetNpoints();

            for(i = 0; i < nvariables; ++i)
            {
                Vmath::Vcopy(npoints,inarray[i],1,outarray[i],1);
            }
            SetBoundaryConditions(outarray,time);
        }
        break;
        case MultiRegions::eGalerkin:
        case MultiRegions::eMixed_CG_Discontinuous:
        {
            ASSERTL0(false,"No Continuous Galerkin for Euler equations");
            break;
        }
        default:
            ASSERTL0(false,"Unknown projection scheme");
            break;
        }
    }

    //----------------------------------------------------
    void EulerCFE::SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble> > &inarray, NekDouble time)
    {

        int nvariables = m_fields.num_elements();
        int nq         = inarray[0].num_elements();
        int cnt = 0;

        // loop over Boundary Regions
        for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
        {
            // Wall Boundary Condition
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eWall)
            {
                WallBoundary(n,cnt,inarray);
            }

            // Symmetric Boundary Condition
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eSymmetry)
            {
                SymmetryBoundary(n,cnt,inarray);
            }

            // Time Dependent Boundary Condition (specified in meshfile)
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eTimeDependent)
            {
                for (int i = 0; i < nvariables; ++i)
                {
                    m_fields[i]->EvaluateBoundaryConditions(time);
                }
            }

            // Isentropic Vortex Boundary Condition
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eIsentropicVortex)
            {
                SetBoundaryIsentropicVortex(n,time,cnt,inarray);
            }

            // Ringleb Flow Inflow and outflow Boundary Condition
            if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eRinglebFlow)
            {
                SetBoundaryRinglebFlow(n,time,cnt,inarray);
            }

            cnt +=m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
        }
    }

    void EulerCFE::v_EvaluateExactSolution(unsigned int field,
            Array<OneD, NekDouble> &outfield,
            const NekDouble time)
    {
        switch(m_problemType)
        {
        case eIsentropicVortex:
        {
            GetExactIsenTropicVortex(field, outfield, time);
        }
        break;
        case eRinglebFlow:
        {
            GetExactRinglebFlow(field, outfield);
        }
        break;
        default:
        {
    
        }
        break;
        }
    }
    
    void EulerCFE::GetExactIsenTropicVortex(int field, Array<OneD, NekDouble> &outarray, NekDouble time)
    {
        int nTotQuadPoints  = GetTotPoints();
    
        Array<OneD, NekDouble> rho(nTotQuadPoints,100.0);
        Array<OneD, NekDouble> rhou(nTotQuadPoints);
        Array<OneD, NekDouble> rhov(nTotQuadPoints);
        Array<OneD, NekDouble> E(nTotQuadPoints);
        Array<OneD, NekDouble> x(nTotQuadPoints);
        Array<OneD, NekDouble> y(nTotQuadPoints);
        Array<OneD, NekDouble> z(nTotQuadPoints);
    
        m_fields[0]->GetCoords(x,y,z);
    
        //---------------------------------
        // flow parameters
        NekDouble x0   = 5.0;
        NekDouble y0   = 0.0;
        NekDouble beta  = 5.0;
        NekDouble u0    = 1.0;
        NekDouble v0    = 0.0;
        NekDouble gamma = m_gamma;
        NekDouble r;
    
        for (int i = 0; i < nTotQuadPoints; ++i)
        {
            r       = sqrt( pow(x[i]-u0*time-x0, 2.0) + pow(y[i]-v0*time-y0, 2.0));
            rho[i]  = pow( (1.0-((gamma-1.0)/(16.0*gamma*M_PI*M_PI))*beta*beta*exp(2.0*(1.0-r*r))), (1.0/(gamma-1.0)) );
            rhou[i] = rho[i] * (1.0 - beta*exp(1.0-r*r)*((y[i]-y0)/(2.0*M_PI)));
            rhov[i] = rho[i] * (beta*exp(1.0-r*r)*((x[i]-x0)/(2.0*M_PI)));
            E[i]    = (pow(rho[i],gamma)/(gamma-1.0)) + 0.5*rho[i]*(pow(rhou[i]/rho[i],2.0)+pow(rhov[i]/rho[i],2.0));
        }
    
        switch (field){
        case 0:
            outarray = rho;
            break;
        case 1:
            outarray = rhou;
            break;
        case 2:
            outarray = rhov;
            break;
        case 3:
            outarray = E;
            break;
        }
    
    }
    
    void EulerCFE::SetIsenTropicVortex(NekDouble initialtime)
    {
        int nTotQuadPoints  = GetTotPoints();

        Array<OneD, NekDouble> rho(nTotQuadPoints,100.0);
        Array<OneD, NekDouble> rhou(nTotQuadPoints);
        Array<OneD, NekDouble> rhov(nTotQuadPoints);
        Array<OneD, NekDouble> E(nTotQuadPoints);
        Array<OneD, NekDouble> x(nTotQuadPoints);
        Array<OneD, NekDouble> y(nTotQuadPoints);
        Array<OneD, NekDouble> z(nTotQuadPoints);

        m_fields[0]->GetCoords(x,y,z);

        //---------------------------------
        // flow parameters
        NekDouble x0   = 5.0;
        NekDouble y0   = 0.0;
        NekDouble beta  = 5.0;
        NekDouble u0    = 1.0;
        NekDouble v0    = 0.0;
        NekDouble gamma = m_gamma;
        NekDouble time  = initialtime;
        NekDouble r;

        for (int i = 0; i < nTotQuadPoints; ++i)
        {
            r       = sqrt( pow(x[i]-u0*time-x0, 2.0) + pow(y[i]-v0*time-y0, 2.0));
            rho[i]  = pow( (1.0-((gamma-1.0)/(16.0*gamma*M_PI*M_PI))*beta*beta*exp(2.0*(1.0-r*r))), (1.0/(gamma-1.0)) );
            rhou[i] = rho[i] * (1.0 - beta*exp(1.0-r*r)*((y[i]-y0)/(2.0*M_PI)));
            rhov[i] = rho[i] * (beta*exp(1.0-r*r)*((x[i]-x0)/(2.0*M_PI)));
            E[i]    = (pow(rho[i],gamma)/(gamma-1.0)) + 0.5*rho[i]*(pow(rhou[i]/rho[i],2.0)+pow(rhov[i]/rho[i],2.0));
        }

        m_fields[0]->SetPhys(rho);
        m_fields[1]->SetPhys(rhou);
        m_fields[2]->SetPhys(rhov);
        m_fields[3]->SetPhys(E);
    
        // forward transform to fill the modal coeffs
        for(int i = 0; i < m_fields.num_elements(); ++i)
        {
            m_fields[i]->SetPhysState(true);
            m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
        }
    }
    
    
    void EulerCFE::SetBoundaryIsentropicVortex(int bcRegion, NekDouble time,int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int nvariables      = physarray.num_elements();
        int nTraceNumPoints = GetTraceTotPoints();
        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);

        // get physical values of the forward trace (from exp to phys)
        for (int i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
            m_fields[i]->ExtractTracePhys(physarray[i],Fwd[i]);
        }

        for(int e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
        {

            int npoints = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
            int id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e);
            int id2 = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondTraceToGlobalTraceMap(cnt++));

            Array<OneD,NekDouble> x(npoints,0.0);
            Array<OneD,NekDouble> y(npoints,0.0);
            Array<OneD,NekDouble> z(npoints,0.0);

            m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetCoords(x,y,z);

            //---------------------------------
            // flow parameters
            NekDouble x0   = 5.0;
            NekDouble y0   = 0.0;
            NekDouble beta  = 5.0;
            NekDouble u0    = 1.0;
            NekDouble v0    = 0.0;
            NekDouble gamma = m_gamma;
            //NekDouble time  = m_time;
            NekDouble r;

            for(int i = 0; i < npoints; i++)
            {
                int kk = id2+i;
                r           = sqrt( pow(x[i]-u0*time-x0, 2.0) + pow(y[i]-v0*time-y0, 2.0));
                Fwd[0][kk]  = pow( (1.0-((gamma-1.0)/(16.0*gamma*M_PI*M_PI))*beta*beta*exp(2.0*(1.0-r*r))), (1.0/(gamma-1.0)) );
                Fwd[1][kk] = Fwd[0][kk] * (1.0 - beta*exp(1.0-r*r)*((y[i]-y0)/(2.0*M_PI)));
                Fwd[2][kk] = Fwd[0][kk] * (beta*exp(1.0-r*r)*((x[i]-x0)/(2.0*M_PI)));
                Fwd[3][kk] = (pow(Fwd[0][kk],gamma)/(gamma-1.0)) + 0.5*Fwd[0][kk]*(pow(Fwd[1][kk]/Fwd[0][kk],2.0)+pow(Fwd[2][kk]/Fwd[0][kk],2.0));
            }
            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Vcopy(npoints,&Fwd[i][id2], 1,&(m_fields[i]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
            }
        }
    }

    void EulerCFE::GetExactRinglebFlow(int field, Array<OneD, NekDouble> &outarray)
    {
        int nTotQuadPoints  = GetTotPoints();

        Array<OneD, NekDouble> rho(nTotQuadPoints,100.0);
        Array<OneD, NekDouble> rhou(nTotQuadPoints);
        Array<OneD, NekDouble> rhov(nTotQuadPoints);
        Array<OneD, NekDouble> E(nTotQuadPoints);
        Array<OneD, NekDouble> x(nTotQuadPoints);
        Array<OneD, NekDouble> y(nTotQuadPoints);
        Array<OneD, NekDouble> z(nTotQuadPoints);

        m_fields[0]->GetCoords(x,y,z);


        //---------------------------------
        // flow parameters
        NekDouble theta    = M_PI/4.0;
        NekDouble kExt     = 0.7;
        NekDouble V        = kExt*sin(theta);
        NekDouble dV,dtheta;
        NekDouble toll     = 1.0e-8;
        NekDouble errV     = 1.0;
        NekDouble errTheta = 1.0;
        NekDouble gamma = m_gamma;
        NekDouble gamma_1_2 = (gamma-1.0)/2.0;
        NekDouble c,k,phi,r,J,VV,pp,sint,P,ss;
        NekDouble J11,J12,J21,J22,det;
        NekDouble Fx, Fy;
        NekDouble xi,yi;
        NekDouble par1;

        for (int i = 0; i < nTotQuadPoints; ++i)
        {
            while((abs(errV) > toll) || (abs(errTheta)>toll))
            {
                VV = V*V;
                sint = sin(theta);
                c = sqrt(1.0-gamma_1_2*VV);
                k = V/sint;
                phi = 1.0/k;
                pp = phi*phi;
                J = 1.0/c + 1.0/(3.0*c*c*c) + 1.0/(5.0*c*c*c*c*c) - 0.5*log((1.0+c)/(1.0-c));
                r = pow(c,1.0/gamma_1_2);
                xi = 1.0/(2.0*r)*(1.0/VV-2.0*pp)+J/2.0;
                yi = phi/(r*V)*sqrt(1.0-VV*pp);
                par1 = 25.0-5.0*VV;
                ss = sint*sint;

                Fx = xi - x[i];
                Fy = yi - y[i];

                J11 = 39062.5/pow(par1,3.5)*(1.0/VV-2.0/VV*ss)*V+1562.5/pow(par1,2.5)*(-2.0/(VV*V)+4.0/(VV*V)*ss)+12.5/pow(par1,1.5)*V+312.5/pow(par1,2.5)*V+7812.5/pow(par1,3.5)*V-0.25*(-1.0/pow(par1,0.5)*V/(1.0-0.2*pow(par1,0.5))-(1.0+0.2*pow(par1,0.5))/pow((1.0-0.2*pow(par1,0.5)),2.0)/pow(par1,0.5)*V)/(1.0+0.2*pow(par1,0.5))*(1.0-0.2*pow(par1,0.5));
                J12 = -6250.0/pow(par1,2.5)/VV*sint*cos(theta);
                J21 = -6250.0/(VV*V)*sint/pow(par1,2.5)*pow((1.0-ss),0.5)+78125.0/V*sint/pow(par1,3.5)*pow((1.0-ss),0.5);
                // the matrix is singular when theta = pi/2
                if(abs(y[i])<toll && abs(cos(theta))<toll)
                {
                    J22 = -39062.5/pow(par1,3.5)/V+3125/pow(par1,2.5)/(VV*V)+12.5/pow(par1,1.5)*V+312.5/pow(par1,2.5)*V+7812.5/pow(par1,3.5)*V-0.25*(-1.0/pow(par1,0.5)*V/(1.0-0.2*pow(par1,0.5))-(1.0+0.2*pow(par1,0.5))/pow((1.0-0.2*pow(par1,0.5)),2.0)/pow(par1,0.5)*V)/(1.0+0.2*pow(par1,0.5))*(1.0-0.2*pow(par1,0.5));

                    // dV = -dV/dx * Fx
                    dV = -1.0/J22*Fx;
                    dtheta = 0.0;
                    theta = M_PI/2.0;
                }
                else
                {
                    J22 = 3125.0/VV*cos(theta)/pow(par1,2.5)*pow((1.0-ss),0.5)-3125.0/VV*ss/pow(par1,2.5)/pow((1.0-ss),0.5)*cos(theta);
                    det = -1.0/(J11*J22-J12*J21);

                    // [dV dtheta]' = -[invJ]*[Fx Fy]'
                    dV     = det*( J22*Fx-J12*Fy);
                    dtheta = det*(-J21*Fx+J11*Fy);
                }

                V = V + dV;
                theta = theta + dtheta;
    
                errV     = abs(dV);
                errTheta = abs(dtheta);
    
            }
    
            c = sqrt(1.0-gamma_1_2*VV);
            r = pow(c,1.0/gamma_1_2);
    
            rho[i] = r;
            rhou[i] = rho[i] * V * cos(theta);
            rhov[i] = rho[i] * V * sin(theta);
            P = (c*c)*rho[i]/gamma;
            E[i]    = P/(gamma-1.0) + 0.5*(rhou[i]*rhou[i]/rho[i]+rhov[i]*rhov[i]/rho[i]);
    
            // Resetting the guess value
            errV = 1.0;
            errTheta = 1.0;
            theta    = M_PI/4.0;
            V        = kExt*sin(theta);
        }
    
        switch (field){
        case 0:
            outarray = rho;
            break;
        case 1:
            outarray = rhou;
            break;
        case 2:
            outarray = rhov;
            break;
        case 3:
            outarray = E;
            break;
        default:
            ASSERTL0(false,"Error in variable number!");
            break;
        }
    }
    
    void EulerCFE::SetBoundaryRinglebFlow(int bcRegion, NekDouble time, int cnt, Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int nvariables      = physarray.num_elements();
        int nTraceNumPoints = GetTraceTotPoints();
        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);

        // get physical values of the forward trace (from exp to phys)
        for (int i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTraceNumPoints);
            m_fields[i]->ExtractTracePhys(physarray[i],Fwd[i]);
        }

        for(int e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
        {

            int npoints = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
            int id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e) ;
            //int id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondCoeffsToGlobalCoeffsMap(cnt+e));
            int id2 = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondTraceToGlobalTraceMap(cnt++));

            Array<OneD,NekDouble> x0(npoints,0.0);
            Array<OneD,NekDouble> x1(npoints,0.0);
            Array<OneD,NekDouble> x2(npoints,0.0);

            m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetCoords(x0,x1,x2);

            //---------------------------------
            // flow parameters
            NekDouble theta    = M_PI/4.0;
            NekDouble kExt     = 0.7;
            NekDouble V        = kExt*sin(theta);
            NekDouble dV,dtheta;
            NekDouble toll     = 1.0e-8;
            NekDouble errV     = 1.0;
            NekDouble errTheta = 1.0;
            NekDouble gamma = m_gamma;
            NekDouble gamma_1_2 = (gamma-1.0)/2.0;
            NekDouble c,k,phi,J,VV,pp,sint,P,ss;
            NekDouble J11,J12,J21,J22,det;
            NekDouble Fx, Fy;
            NekDouble xi,yi;
            NekDouble par1;
            NekDouble r;

            // Loop on all the points of that edge
            for(int j = 0; j < npoints; j++)
            {

                while((abs(errV) > toll) || (abs(errTheta)>toll))
                {
                    VV = V*V;
                    sint = sin(theta);
                    c = sqrt(1.0-gamma_1_2*VV);
                    k = V/sint;
                    phi = 1.0/k;
                    pp = phi*phi;
                    J = 1.0/c + 1.0/(3.0*c*c*c) + 1.0/(5.0*c*c*c*c*c) - 0.5*log((1.0+c)/(1.0-c));
                    r = pow(c,1.0/gamma_1_2);
                    xi = 1.0/(2.0*r)*(1.0/VV-2.0*pp)+J/2.0;
                    yi = phi/(r*V)*sqrt(1.0-VV*pp);
                    par1 = 25.0-5.0*VV;
                    ss = sint*sint;

                    Fx = xi - x0[j];
                    Fy = yi - x1[j];

                    J11 = 39062.5/pow(par1,3.5)*(1.0/VV-2.0/VV*ss)*V+1562.5/pow(par1,2.5)*(-2.0/(VV*V)+4.0/(VV*V)*ss)+12.5/pow(par1,1.5)*V+312.5/pow(par1,2.5)*V+7812.5/pow(par1,3.5)*V-0.25*(-1.0/pow(par1,0.5)*V/(1.0-0.2*pow(par1,0.5))-(1.0+0.2*pow(par1,0.5))/pow((1.0-0.2*pow(par1,0.5)),2.0)/pow(par1,0.5)*V)/(1.0+0.2*pow(par1,0.5))*(1.0-0.2*pow(par1,0.5));
                    J12 = -6250.0/pow(par1,2.5)/VV*sint*cos(theta);
                    J21 = -6250.0/(VV*V)*sint/pow(par1,2.5)*pow((1.0-ss),0.5)+78125.0/V*sint/pow(par1,3.5)*pow((1.0-ss),0.5);
                    // the matrix is singular when theta = pi/2
                    if(abs(x1[j])<toll && abs(cos(theta))<toll)
                    {
                        J22 = -39062.5/pow(par1,3.5)/V+3125/pow(par1,2.5)/(VV*V)+12.5/pow(par1,1.5)*V+312.5/pow(par1,2.5)*V+7812.5/pow(par1,3.5)*V-0.25*(-1.0/pow(par1,0.5)*V/(1.0-0.2*pow(par1,0.5))-(1.0+0.2*pow(par1,0.5))/pow((1.0-0.2*pow(par1,0.5)),2.0)/pow(par1,0.5)*V)/(1.0+0.2*pow(par1,0.5))*(1.0-0.2*pow(par1,0.5));

                        // dV = -dV/dx * Fx
                        dV = -1.0/J22*Fx;
                        dtheta = 0.0;
                        theta = M_PI/2.0;
                    }
                    else
                    {
                        J22 = 3125.0/VV*cos(theta)/pow(par1,2.5)*pow((1.0-ss),0.5)-3125.0/VV*ss/pow(par1,2.5)/pow((1.0-ss),0.5)*cos(theta);
                        det = -1.0/(J11*J22-J12*J21);

                        // [dV dtheta]' = -[invJ]*[Fx Fy]'
                        dV     = det*( J22*Fx-J12*Fy);
                        dtheta = det*(-J21*Fx+J11*Fy);
                    }

                    V = V + dV;
                    theta = theta + dtheta;

                    errV     = abs(dV);
                    errTheta = abs(dtheta);

                }

                c  = sqrt(1.0-gamma_1_2*VV);
                int kk = id2+j;
                NekDouble timeramp = 200.0;
                std::string restartstr = "RESTART";
                if(time<timeramp &&
                        !(m_session->DefinesFunction("InitialConditions") &&
                          m_session->GetFunctionType("InitialConditions", 0) == LibUtilities::eFunctionTypeFile))
                {
                    Fwd[0][kk]  = pow(c,1.0/gamma_1_2)*exp(-1.0+time/timeramp);
                    Fwd[1][kk] = Fwd[0][kk] * V * cos(theta)*exp(-1+time/timeramp);
                    Fwd[2][kk] = Fwd[0][kk] * V * sin(theta)*exp(-1+time/timeramp);
                }
                else
                {
                    Fwd[0][kk]  = pow(c,1.0/gamma_1_2);
                    Fwd[1][kk] = Fwd[0][kk] * V * cos(theta);
                    Fwd[2][kk] = Fwd[0][kk] * V * sin(theta);
                }
                P  = (c*c)*Fwd[0][kk]/gamma;
                Fwd[3][kk]  = P/(gamma-1.0) + 0.5*(Fwd[1][kk]*Fwd[1][kk]/Fwd[0][kk]+Fwd[2][kk]*Fwd[2][kk]/Fwd[0][kk]);

                errV = 1.0;
                errTheta = 1.0;
                theta    = M_PI/4.0;
                V        = kExt*sin(theta);

            }
            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Vcopy(npoints,&Fwd[i][id2], 1,&(m_fields[i]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
            }
        }
    }

    void EulerCFE::SetInitialRinglebFlow(void)
    {

        int nbnd    = m_fields[0]->GetBndConditions().num_elements();// Get number of different boundaries in the input file

        // Loop on all the edges of the input file
        for(int bcRegion=0; bcRegion < nbnd; ++bcRegion)
        {

            int npoints = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetNpoints();
            Array<OneD,NekDouble> x0(npoints,0.0);
            Array<OneD,NekDouble> x1(npoints,0.0);
            Array<OneD,NekDouble> x2(npoints,0.0);

            Array<OneD, NekDouble> rho(npoints,0.0);
            Array<OneD, NekDouble> rhou(npoints,0.0);
            Array<OneD, NekDouble> rhov(npoints,0.0);
            Array<OneD, NekDouble> E(npoints,0.0);

            m_fields[0]->GetBndCondExpansions()[bcRegion]->GetCoords(x0,x1,x2);
            //---------------------------------
            // flow parameters
            NekDouble theta    = M_PI/4.0;
            NekDouble kExt     = 0.7;
            NekDouble V        = kExt*sin(theta);
            NekDouble dV,dtheta;
            NekDouble toll     = 1.0e-8;
            NekDouble errV     = 1.0;
            NekDouble errTheta = 1.0;
            NekDouble gamma = m_gamma;
            NekDouble gamma_1_2 = (gamma-1.0)/2.0;
            NekDouble c,k,phi,J,VV,pp,sint,P,ss;
            NekDouble J11,J12,J21,J22,det;
            NekDouble Fx, Fy;
            NekDouble xi,yi;
            NekDouble par1;
            NekDouble r;

            // Loop on all the points of that edge
            for(int j = 0; j < npoints; j++)
            {

                while((abs(errV) > toll) || (abs(errTheta)>toll))
                {

                    VV = V*V;
                    sint = sin(theta);
                    c = sqrt(1.0-gamma_1_2*VV);
                    k = V/sint;
                    phi = 1.0/k;
                    pp = phi*phi;
                    J = 1.0/c + 1.0/(3.0*c*c*c) + 1.0/(5.0*c*c*c*c*c) - 0.5*log((1.0+c)/(1.0-c));
                    r = pow(c,1.0/gamma_1_2);
                    xi = 1.0/(2.0*r)*(1.0/VV-2.0*pp)+J/2.0;
                    yi = phi/(r*V)*sqrt(1.0-VV*pp);
                    par1 = 25.0-5.0*VV;
                    ss = sint*sint;

                    Fx = xi - x0[j];
                    Fy = yi - x1[j];

                    J11 = 39062.5/pow(par1,3.5)*(1.0/VV-2.0/VV*ss)*V+1562.5/pow(par1,2.5)*(-2.0/(VV*V)+4.0/(VV*V)*ss)+12.5/pow(par1,1.5)*V+312.5/pow(par1,2.5)*V+7812.5/pow(par1,3.5)*V-0.25*(-1.0/pow(par1,0.5)*V/(1.0-0.2*pow(par1,0.5))-(1.0+0.2*pow(par1,0.5))/pow((1.0-0.2*pow(par1,0.5)),2.0)/pow(par1,0.5)*V)/(1.0+0.2*pow(par1,0.5))*(1.0-0.2*pow(par1,0.5));
                    J12 = -6250.0/pow(par1,2.5)/VV*sint*cos(theta);
                    J21 = -6250.0/(VV*V)*sint/pow(par1,2.5)*pow((1.0-ss),0.5)+78125.0/V*sint/pow(par1,3.5)*pow((1.0-ss),0.5);

                    // the matrix is singular when theta = pi/2
                    if(abs(x1[j])<toll && abs(cos(theta))<toll)
                    {

                        J22 = -39062.5/pow(par1,3.5)/V+3125/pow(par1,2.5)/(VV*V)+12.5/pow(par1,1.5)*V+312.5/pow(par1,2.5)*V+7812.5/pow(par1,3.5)*V-0.25*(-1.0/pow(par1,0.5)*V/(1.0-0.2*pow(par1,0.5))-(1.0+0.2*pow(par1,0.5))/pow((1.0-0.2*pow(par1,0.5)),2.0)/pow(par1,0.5)*V)/(1.0+0.2*pow(par1,0.5))*(1.0-0.2*pow(par1,0.5));

                        // dV = -dV/dx * Fx
                        dV = -1.0/J22*Fx;
                        dtheta = 0.0;
                        theta = M_PI/2.0;
                    }
                    else
                    {

                        J22 = 3125.0/VV*cos(theta)/pow(par1,2.5)*pow((1.0-ss),0.5)-3125.0/VV*ss/pow(par1,2.5)/pow((1.0-ss),0.5)*cos(theta);
                        det = -1.0/(J11*J22-J12*J21);

                        // [dV dtheta]' = -[invJ]*[Fx Fy]'
                        dV     = det*( J22*Fx-J12*Fy);
                        dtheta = det*(-J21*Fx+J11*Fy);
                    }

                    V = V + dV;
                    theta = theta + dtheta;

                    errV     = abs(dV);
                    errTheta = abs(dtheta);

                }

                c  = sqrt(1.0-gamma_1_2*VV);
                rho[j]  = pow(c,1.0/gamma_1_2) * exp(-1.0);
                rhou[j] = rho[j] * V * cos(theta)* exp(-1.0);
                rhov[j] = rho[j] * V * sin(theta)* exp(-1.0);
                P  = (c*c)*rho[j]/gamma;
                E[j]  = P/(gamma-1.0) + 0.5*(rhou[j]*rhou[j]/rho[j]+rhov[j]*rhov[j]/rho[j]);

                errV = 1.0;
                errTheta = 1.0;
                theta    = M_PI/4.0;
                V        = kExt*sin(theta);

            }

            m_fields[0]->GetBndCondExpansions()[bcRegion]->SetPhys(rho);
            m_fields[1]->GetBndCondExpansions()[bcRegion]->SetPhys(rhou);
            m_fields[2]->GetBndCondExpansions()[bcRegion]->SetPhys(rhov);
            m_fields[3]->GetBndCondExpansions()[bcRegion]->SetPhys(E);

            // forward transform to fill the modal coeffs
            for(int i = 0; i < m_fields.num_elements(); ++i)
            {
                m_fields[i]->GetBndCondExpansions()[bcRegion]->FwdTrans_BndConstrained(m_fields[i]->GetBndCondExpansions()[bcRegion]->GetPhys(),m_fields[i]->GetBndCondExpansions()[bcRegion]->UpdateCoeffs());
            }

        }

        // Calculation of the initial internal values as a weighted average over the distance from the boundaries
        int nq = m_fields[0]->GetNpoints();

        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);

        // get the coordinates (assuming all fields have the same discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);

        for(int j = 0; j < nq; j++)
        {
            NekDouble Dist = 0.0;
            NekDouble rho  = 0.0;
            NekDouble rhou = 0.0;
            NekDouble rhov = 0.0;
            NekDouble E    = 0.0;
            NekDouble SumDist = 0.0;

            // Calculation of all the distances
            // Loop on all the edges of the input file
            for(int bcRegion=0; bcRegion < nbnd; ++bcRegion)
            {

                int npoints = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetNpoints(); // get quadrature points on the edge

                Array<OneD,NekDouble> xb0(npoints,0.0);
                Array<OneD,NekDouble> xb1(npoints,0.0);
                Array<OneD,NekDouble> xb2(npoints,0.0);

                m_fields[0]->GetBndCondExpansions()[bcRegion]->GetCoords(xb0,xb1,xb2);

                for(int k=0; k<npoints; k++)
                {
                    Dist = sqrt((xb0[k]-x0[j])*(xb0[k]-x0[j])+(xb1[k]-x1[j])*(xb1[k]-x1[j]));
                    SumDist += Dist;
                    rho  += Dist*(m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys())[k];
                    rhou += Dist*(m_fields[1]->GetBndCondExpansions()[bcRegion]->GetPhys())[k];
                    rhov += Dist*(m_fields[2]->GetBndCondExpansions()[bcRegion]->GetPhys())[k];
                    E    += Dist*(m_fields[3]->GetBndCondExpansions()[bcRegion]->GetPhys())[k];
                }
            }

            rho  = rho/SumDist;
            rhou = rhou/SumDist;
            rhov = rhov/SumDist;
            E    = E/SumDist;

            (m_fields[0]->UpdatePhys())[j] = rho;
            (m_fields[1]->UpdatePhys())[j] = rhou;
            (m_fields[2]->UpdatePhys())[j] = rhov;
            (m_fields[3]->UpdatePhys())[j] = E;

        }

        for(int i = 0 ; i < m_fields.num_elements(); i++)
        {
            m_fields[i]->SetPhysState(true);
            m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
        }
    
        // dump initial conditions to file
        std::string outname = m_sessionName + "_initialRingleb.chk";
        //ofstream outfile(outname.c_str());
        WriteFld(outname);
    }

}
