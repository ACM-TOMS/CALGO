///////////////////////////////////////////////////////////////////////////////
//
// File IncNavierStokes.cpp
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
// Description: Incompressible Navier Stokes class definition built on
// ADRBase class
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>
#include <cstdio>
#include <cstdlib>
#include <iomanip>

namespace Nektar
{

    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */
    IncNavierStokes::IncNavierStokes(const LibUtilities::SessionReaderSharedPtr& pSession):
        //EquationSystem(pSession),
        UnsteadySystem(pSession),
        m_steadyStateSteps(0),
        m_subSteppingScheme(false)
    {
    }

    void IncNavierStokes::v_InitObject()
    {
        
        int i,j;
        int numfields = m_fields.num_elements();
        std::string velids[] = {"u","v","w"};
        
        // Set up Velocity field to point to the first m_expdim of m_fields; 
        m_velocity = Array<OneD,int>(m_spacedim);
        
        for(i = 0; i < m_spacedim; ++i)
        {
            for(j = 0; j < numfields; ++j)
            {
                std::string var = m_boundaryConditions->GetVariable(j);
                if(NoCaseStringCompare(velids[i],var) == 0)
                {
                    m_velocity[i] = j;
                    break;
                }
                
                if(j == numfields)
                {
                    std::string error = "Failed to find field: " + var; 
                    ASSERTL0(false,error.c_str());
                }
            }
        }
        
        // Set up equation type enum using kEquationTypeStr
        for(i = 0; i < (int) eEquationTypeSize; ++i)
        {
            bool match;
            m_session->MatchSolverInfo("EQTYPE",kEquationTypeStr[i],match,false);
            if(match)
            {
                m_equationType = (EquationType)i; 
                break;
            }
        }
        ASSERTL0(i != eEquationTypeSize,"EQTYPE not found in SOLVERINFO section");
        
        // This probably should to into specific implementations 
        // Equation specific Setups 
        switch(m_equationType)
        {
        case eSteadyStokes: 
        case eSteadyOseen: 
        case eSteadyNavierStokes:
        case eSteadyLinearisedNS: 
            break;
        case eUnsteadyNavierStokes:
        case eUnsteadyStokes:
            m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
            m_session->LoadParameter("IO_EnergySteps", m_energysteps, 0);
            m_session->LoadParameter("SteadyStateSteps", m_steadyStateSteps, 0);
            m_session->LoadParameter("SteadyStateTol", m_steadyStateTol, 1e-6);
            
            // check to see if any user defined boundary condition is
            // indeed implemented
            
            for(int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
            {	
                // Time Dependent Boundary Condition (if no user
                // defined then this is empty)
                if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() != SpatialDomains::eNoUserDefined)
                {
                    if (m_fields[0]->GetBndConditions()[n]->GetUserDefined() != SpatialDomains::eTimeDependent)
                    {                     	     
                        if(m_fields[0]->GetBndConditions()[n]->GetUserDefined() != SpatialDomains::eI)
                        {  	 	 
                            ASSERTL0(false,"Unknown USERDEFINEDTYPE boundary condition");
                        }
                    }
                }
            }
            break;
        case eNoEquationType:
        default:
            ASSERTL0(false,"Unknown or undefined equation type");
        }
        
        m_session->LoadParameter("Kinvis", m_kinvis);
        
        if (m_equationType == eUnsteadyNavierStokes || m_equationType == eSteadyNavierStokes)
        {
            std::string vConvectiveType = "Convective";
            if (m_session->DefinesTag("AdvectiveType"))
            {
                vConvectiveType = m_session->GetTag("AdvectiveType");
            }
            m_advObject = GetAdvectionTermFactory().CreateInstance(vConvectiveType, m_session, m_graph);
        }
	
        if (m_equationType == eUnsteadyLinearisedNS)// || m_equationType == eSteadyNavierStokes)
        {
            std::string vConvectiveType = "Linearised";
            if (m_session->DefinesTag("AdvectiveType"))
            {
                //vConvectiveType = m_session->GetTag("Linearised");
                vConvectiveType = m_session->GetTag("AdvectiveType");
            }
            m_advObject = GetAdvectionTermFactory().CreateInstance(vConvectiveType, m_session, m_graph);
        }
		
		if(m_equationType == eUnsteadyStokes)
		{
			std::string vConvectiveType = "NoAdvection";
			m_advObject = GetAdvectionTermFactory().CreateInstance(vConvectiveType, m_session, m_graph);
		}
	
#if 0 // Not required if building on an UnsteadySystem rather than an EquationSystem
        // Set up filters
        LibUtilities::FilterMap::const_iterator x;
        LibUtilities::FilterMap f = m_session->GetFilters();
        for (x = f.begin(); x != f.end(); ++x)
        {
            m_filters.push_back(SolverUtils::GetFilterFactory().CreateInstance(x->first, m_session, x->second));
        }
#endif
    }

    IncNavierStokes::~IncNavierStokes(void)
    {

    }
    
    void IncNavierStokes::AdvanceInTime(int nsteps)
    {
        int i,n;
        int phystot  = m_fields[0]->GetTotPoints();
        int n_fields = m_fields.num_elements();
        static int nchk = 0;
		
        if(m_HomogeneousType == eHomogeneous1D)
        {
            for(i = 0; i < n_fields; ++i)
            {
                m_fields[i]->HomogeneousFwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdatePhys());
                m_fields[i]->SetWaveSpace(true);
                m_fields[i]->SetPhysState(false);
            }
        }
	
        // Set up wrapper to fields data storage. 
        Array<OneD, Array<OneD, NekDouble> >  fields(m_nConvectiveFields);
        for(i = 0; i < m_nConvectiveFields; ++i)
        {
            fields[i]  = m_fields[i]->UpdatePhys();
        }
		
        // Initialise NS solver which is set up to use a GLM method
        // with calls to EvaluateAdvection_SetPressureBCs and
        // SolveUnsteadyStokesSystem
        m_integrationSoln = m_integrationScheme[m_intSteps-1]->InitializeScheme(m_timestep, fields, m_time, m_integrationOps);
		
        std::string   mdlname = m_session->GetSessionName() + ".mdl";
        std::ofstream mdlFile;

        if (m_energysteps && m_comm->GetRank() == 0)
        {
            mdlFile.open(mdlname.c_str());
        }
        
        std::vector<SolverUtils::FilterSharedPtr>::iterator x;
        for (x = m_filters.begin(); x != m_filters.end(); ++x)
        {
            (*x)->Initialise(m_fields, m_time);
        }

        //Time advance
        for(n = 0; n < nsteps; ++n)
        {

            if(m_subSteppingScheme)
            {
                SubStepSaveFields(n);
                SubStepAdvance(n);
            }

            // Advance velocity fields
            fields = m_integrationScheme[min(n,m_intSteps-1)]->TimeIntegrate(m_timestep, m_integrationSoln, m_integrationOps);
            
            m_time += m_timestep;
            
            // Write out current time step
            if(m_infosteps && !((n+1)%m_infosteps) && m_comm->GetRank() == 0)
            {
                cout << "Step: " << n+1 << "  Time: " << m_time << endl;
            }

            // Write out energy data to file
            if(m_energysteps && !((n+1)%m_energysteps))
            {
                if(m_HomogeneousType != eNotHomogeneous)
                {
                    if(m_HomogeneousType == eHomogeneous1D)
                    {
                        int colrank = m_comm->GetColumnComm()->GetRank();
                        int nproc   = m_comm->GetColumnComm()->GetSize();
                        int locsize = m_npointsZ/nproc/2;
                        int ensize  = m_npointsZ/nproc/2;
                        
                        Array<OneD, NekDouble> energy    (locsize,0.0);
                        Array<OneD, NekDouble> energy_tmp(locsize,0.0);
                        Array<OneD, NekDouble> tmp;
			
                        // Calculate modal energies.
                        for(i = 0; i < m_nConvectiveFields; ++i)
                        {
                            energy_tmp = m_fields[i]->HomogeneousEnergy();
                            Vmath::Vadd(locsize,energy_tmp,1,energy,1,energy,1);
                        }
                        
                        // Send to root process.
                        if (colrank == 0)
                        {
                            int j, m = 0;
                            
                            for (j = 0; j < energy.num_elements(); ++j, ++m)
                            {
                                mdlFile << setw(10) << m_time 
                                        << setw(5)  << m
                                        << setw(18) << energy[j] << endl;
                            }
                            
                            for (i = 1; i < nproc; ++i)
                            {
                                m_comm->GetColumnComm()->Recv(i, energy);
                                for (j = 0; j < energy.num_elements(); ++j, ++m)
                                {
                                    mdlFile << setw(10) << m_time 
                                            << setw(5)  << m
                                            << setw(18) << energy[j] << endl;
                                }
                            }
                        }
                        else
                        {
                            m_comm->GetColumnComm()->Send(0, energy);
                        }
                    }
                    else
                    {
                        ASSERTL0(false,"3D Homogeneous 2D energy dumping not implemented yet");
                    }
                }
                else
                {
                    NekDouble energy = 0.0;
                    for(i = 0; i < m_nConvectiveFields; ++i)
                    {
                        m_fields[i]->SetPhys(fields[i]);
                        m_fields[i]->SetPhysState(true);
                        NekDouble norm = L2Error(i, true);
                        energy += norm*norm;
                    }
                    mdlFile << m_time << "   " << 0.5*energy << endl;
                }
                
            }
            
            // dump data in m_fields->m_coeffs to file. 
            if(m_checksteps && n&&(!((n+1)%m_checksteps)))
            {
                if(m_HomogeneousType == eHomogeneous1D)
                {
                    for(i = 0; i< n_fields; i++)
                    {
                        m_fields[i]->SetWaveSpace(false);
                        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),m_fields[i]->UpdatePhys());
                        m_fields[i]->SetPhysState(true);
                    }
                    nchk++;
                    Checkpoint_Output(nchk);
                    for(i = 0; i< n_fields; i++)
                    {
                        m_fields[i]->SetWaveSpace(true);
                        m_fields[i]->HomogeneousFwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdatePhys());
                        m_fields[i]->SetPhysState(false);
                    }
                }
                else
                {
                    for(i = 0; i < m_nConvectiveFields; ++i)
                    {
                        m_fields[i]->SetPhys(fields[i]);
                        m_fields[i]->SetPhysState(true);
                    }
                    nchk++;
                    Checkpoint_Output(nchk);
                }
            }
            
            
            if(m_steadyStateSteps && n && (!((n+1)%m_steadyStateSteps)))
            {
                if(CalcSteadyState() == true)
                {
                    cout << "Reached Steady State to tolerance " << m_steadyStateTol << endl; 
                    break;
                }
            }

            // Transform data into coefficient space
            if (m_filters.size() > 0)
            {
                for (i = 0; i < m_nConvectiveFields; ++i)
                {
                    m_fields[i]->FwdTrans_IterPerExp(fields[i],
                                                     m_fields[i]->UpdateCoeffs());
                    m_fields[i]->SetPhysState(false);
                }
            }

            std::vector<SolverUtils::FilterSharedPtr>::iterator x;
            for (x = m_filters.begin(); x != m_filters.end(); ++x)
            {
                (*x)->Update(m_fields, m_time);
            }
        }
	
        if(m_HomogeneousType == eHomogeneous1D)
        {
            for(i = 0 ; i< n_fields ; i++)
            {
                m_fields[i]->SetWaveSpace(false);
				m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),m_fields[i]->UpdatePhys());
                m_fields[i]->SetPhysState(true);
            }
        }
        else 
        {
            for(i = 0; i < m_nConvectiveFields; ++i)
            {
                m_fields[i]->SetPhys(fields[i]);
                m_fields[i]->SetPhysState(true);
            }
        }
	
        
        if (m_energysteps)
        {
            mdlFile.close();
        }
        
        for (x = m_filters.begin(); x != m_filters.end(); ++x)
        {
            (*x)->Finalise(m_fields, m_time);
        }
    }
    

    void IncNavierStokes::SubStepSaveFields(const int nstep)
    {
        int i,n;
        int nvel = m_velocity.num_elements();
        int npts = m_fields[0]->GetTotPoints();

        // rotate fields 
        int nblocks = m_previousVelFields.num_elements()/nvel; 
        Array<OneD, NekDouble> save; 
        
        for(n = 0; n < nvel; ++n)
        {
            save = m_previousVelFields[(nblocks-1)*nvel+n];
            
            for(i = nblocks-1; i > 0; --i)
            {
                m_previousVelFields[i*nvel+n] = m_previousVelFields[(i-1)*nvel+n];
            }
            
            m_previousVelFields[n] = save; 
        }
        
        // Put previous field 
        for(i = 0; i < nvel; ++i)
        {
            m_fields[m_velocity[i]]->BwdTrans(m_fields[m_velocity[i]]->GetCoeffs(),
                                              m_fields[m_velocity[i]]->UpdatePhys());
            Vmath::Vcopy(npts,m_fields[m_velocity[i]]->GetPhys(),1,
                         m_previousVelFields[i],1);   
        }

        if(nstep == 0)// initialise all levels with first field 
        {
            for(n = 0; n < nvel; ++n)
            {
                for(i = 1; i < nblocks; ++i)
                {
                    Vmath::Vcopy(npts,m_fields[m_velocity[n]]->GetPhys(),1,
                                 m_previousVelFields[i*nvel+n],1);   
                    
                }
            }
        }
    }


    void IncNavierStokes::SubStepAdvance(const int nstep)
    {
        int n;
        NekDouble dt; 
        int nsubsteps,minsubsteps;
        NekDouble time = m_time;
        Array<OneD, Array<OneD, NekDouble> > fields,velfields;
        static int ncalls = 1;
        int  nint    = min(ncalls++,m_intSteps);
        Array<OneD, NekDouble> CFL(m_fields[0]->GetExpSize(),m_cfl);
        dt = GetTimeStep(m_fields[0]->EvalBasisNumModesMaxPerExp(),CFL,0.0);

        nsubsteps = (m_timestep > dt)? ((int)(m_timestep/dt)+1):1; 
        m_session->LoadParameter("MinSubSteps",minsubsteps,0);
        nsubsteps = max(minsubsteps,nsubsteps);

        dt = m_timestep/nsubsteps;
        
        if(m_infosteps && !((nstep+1)%m_infosteps) && m_comm->GetRank() == 0)
        {
            cout << "Sub-integrating using "<< nsubsteps << " steps over Dt = " 
                 << m_timestep << " (SubStep CFL=" << m_cfl << ")"<< endl;
        }

        for(int m = 0; m < nint; ++m)
        {
            // We need to update the fields held by the m_integrationSoln
            fields = m_integrationSoln->UpdateSolutionVector()[m];
            
            // Initialise NS solver which is set up to use a GLM method
            // with calls to EvaluateAdvection_SetPressureBCs and
            // SolveUnsteadyStokesSystem
            LibUtilities::TimeIntegrationSolutionSharedPtr 
                SubIntegrationSoln = m_subStepIntegrationScheme->InitializeScheme(dt, fields, time, m_subStepIntegrationOps);
            
            for(n = 0; n < nsubsteps; ++n)
            {
                fields = m_subStepIntegrationScheme->TimeIntegrate(dt, SubIntegrationSoln, m_subStepIntegrationOps);
            }
            
            // Reset time integrated solution in m_integrationSoln 
            m_integrationSoln->SetSolVector(m,fields);
        }
    }


    /** 
     * Explicit Advection terms used by SubStepAdvance time integration
     */
    void IncNavierStokes::SubStepAdvection(const Array<OneD, const Array<OneD, NekDouble> > &inarray,  Array<OneD, Array<OneD, NekDouble> > &outarray, const NekDouble time)
    {
        int i;
        int nVariables     = inarray.num_elements();
        int nQuadraturePts = inarray[0].num_elements();

        /// Get the number of coefficients
        int ncoeffs = m_fields[0]->GetNcoeffs(); 
        
        /// Define an auxiliary variable to compute the RHS 
        Array<OneD, Array<OneD, NekDouble> > WeakAdv(nVariables);
        WeakAdv[0] = Array<OneD, NekDouble> (ncoeffs*nVariables);
        for(i = 1; i < nVariables; ++i)
        {
            WeakAdv[i] = WeakAdv[i-1] + ncoeffs;
        }

        //cout << "Time: " << time << endl;
        Array<OneD, Array<OneD, NekDouble> > Velfields(m_velocity.num_elements());
        Array<OneD, int> VelIds(m_velocity.num_elements());

        Velfields[0] = Array<OneD, NekDouble> (nQuadraturePts*m_velocity.num_elements());
        
        for(i = 1; i < m_velocity.num_elements(); ++i)
        {
            Velfields[i] = Velfields[i-1] + nQuadraturePts; 
        }
        SubStepExtrapoloteField(fmod(time,m_timestep), Velfields);

        m_advObject->DoAdvection(m_fields, Velfields, inarray, outarray, time);
        
        for(i = 0; i < nVariables; ++i)
        {
            m_fields[i]->IProductWRTBase(outarray[i],WeakAdv[i]); 
            // negation requried due to sign of DoAdvection term to be consistent
           Vmath::Neg(ncoeffs, WeakAdv[i], 1);
        }
        
        AddAdvectionPenaltyFlux(Velfields, inarray, WeakAdv);
        
        /// Operations to compute the RHS
        for(i = 0; i < nVariables; ++i)
        {
            // Negate the RHS
            Vmath::Neg(ncoeffs, WeakAdv[i], 1);

            /// Multiply the flux by the inverse of the mass matrix
            m_fields[i]->MultiplyByElmtInvMass(WeakAdv[i], WeakAdv[i]);
            
            /// Store in outarray the physical values of the RHS
            m_fields[i]->BwdTrans(WeakAdv[i], outarray[i]);
        }
        
        //add the force
        if(m_session->DefinesFunction("BodyForce"))
        {
            if(m_SingleMode || m_HalfMode)
            {
                for(int i = 0; i < m_nConvectiveFields; ++i)
                {
                    m_forces[i]->SetWaveSpace(true);	
                    m_forces[i]->BwdTrans(m_forces[i]->GetCoeffs(),
                                          m_forces[i]->UpdatePhys());
                }
            }
 
            int nqtot = m_fields[0]->GetTotPoints();
            for(int i = 0; i < m_nConvectiveFields; ++i)
            {
                Vmath::Vadd(nqtot,outarray[i],1,(m_forces[i]->GetPhys()),1,outarray[i],1);
            }
        }
    }



    void IncNavierStokes::v_GetFluxVector(const int i, 
                                          Array<OneD, Array<OneD, NekDouble> > &physfield,
                                            Array<OneD, Array<OneD, NekDouble> > &flux)
    {
        ASSERTL1(flux.num_elements() == m_velocity.num_elements(),"Dimension of flux array and velocity array do not match");

        for(int j = 0; j < flux.num_elements(); ++j)
        {
            Vmath::Vmul(GetNpoints(), physfield[i], 1, m_fields[m_velocity[j]]->GetPhys(), 1, flux[j], 1);
        }
    }

    
    
    void IncNavierStokes::v_NumericalFlux(Array<OneD, Array<OneD, NekDouble> > &physfield, 
                                          Array<OneD, Array<OneD, NekDouble> > &numflux)
    {
        /// Counter variable
        int i;

        /// Number of trace points
        int nTracePts   = GetTraceNpoints();
        
        /// Number of spatial dimensions
        int nDimensions = m_spacedim;

        /// Forward state array
        Array<OneD, NekDouble> Fwd(2*nTracePts);
        
        /// Backward state array
        Array<OneD, NekDouble> Bwd = Fwd + nTracePts;
        
        /// Normal velocity array
        Array<OneD, NekDouble> Vn (nTracePts, 0.0);
        
        // Extract velocity field along the trace space and multiply by trace normals
        for(i = 0; i < nDimensions; ++i)
        {
            m_fields[0]->ExtractTracePhys(m_fields[m_velocity[i]]->GetPhys(), Fwd);
            Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, Fwd, 1, Vn, 1, Vn, 1);
        }

        /// Compute the numerical fluxes at the trace points
        for(i = 0; i < numflux.num_elements(); ++i)
        {
            /// Extract forwards/backwards trace spaces
            m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd, Bwd);

            /// Upwind between elements
            m_fields[i]->GetTrace()->Upwind(Vn, Fwd, Bwd, numflux[i]);

            /// Calculate the numerical fluxes multipling Fwd or Bwd
            /// by the normal advection velocity
            Vmath::Vmul(nTracePts, numflux[i], 1, Vn, 1, numflux[i], 1);
        }
    }


    void IncNavierStokes::AddAdvectionPenaltyFlux(const Array<OneD, const Array<OneD, NekDouble> > &velfield, 
                                                  const Array<OneD, const Array<OneD, NekDouble> > &physfield, 
                                                  Array<OneD, Array<OneD, NekDouble> > &Outarray)
    {
        ASSERTL1(physfield.num_elements() == Outarray.num_elements(),"Physfield and outarray are of different dimensions");

        int i;

        /// Number of trace points
        int nTracePts   = GetTraceNpoints();
        
        /// Number of spatial dimensions
        int nDimensions = m_spacedim;

        /// Forward state array
        Array<OneD, NekDouble> Fwd(3*nTracePts);
        
        /// Backward state array
        Array<OneD, NekDouble> Bwd = Fwd + nTracePts;

        /// upwind numerical flux state array
        Array<OneD, NekDouble> numflux = Bwd + nTracePts;
        
        /// Normal velocity array
        Array<OneD, NekDouble> Vn (nTracePts, 0.0);
        
        // Extract velocity field along the trace space and multiply by trace normals
        for(i = 0; i < nDimensions; ++i)
        {
            m_fields[0]->ExtractTracePhys(m_fields[m_velocity[i]]->GetPhys(), Fwd);
            Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, Fwd, 1, Vn, 1, Vn, 1);
        }
        
        for(i = 0; i < physfield.num_elements(); ++i)
        {
            /// Extract forwards/backwards trace spaces
            /// Note: Needs to have correct i value to get boundary conditions
            m_fields[i]->GetFwdBwdTracePhys(physfield[i], Fwd, Bwd);
            
            /// Upwind between elements
            m_fields[0]->GetTrace()->Upwind(Vn, Fwd, Bwd, numflux);

            /// Construct difference between numflux and Fwd,Bwd
            Vmath::Vsub(nTracePts, numflux, 1, Fwd, 1, Fwd, 1);
            Vmath::Vsub(nTracePts, numflux, 1, Bwd, 1, Bwd, 1);

            /// Calculate the numerical fluxes multipling Fwd, Bwd and
            /// numflux by the normal advection velocity
            Vmath::Vmul(nTracePts, Fwd, 1, Vn, 1, Fwd, 1);
            Vmath::Vmul(nTracePts, Bwd, 1, Vn, 1, Bwd, 1);

            m_fields[0]->AddFwdBwdTraceIntegral(Fwd,Bwd,Outarray[i]);
        }
    }

    
    /** 
     * Projection used by SubStepAdvance time integration
     */
    void IncNavierStokes::SubStepProjection(const Array<OneD, const Array<OneD, NekDouble> > &inarray,  Array<OneD, Array<OneD, NekDouble> > &outarray, const NekDouble time)
    {
        ASSERTL1(inarray.num_elements() == outarray.num_elements(),"Inarray and outarray of different sizes ");

        for(int i = 0; i < inarray.num_elements(); ++i)
        {
            Vmath::Vcopy(inarray[i].num_elements(),inarray[i],1,outarray[i],1);
        }
    }
    
    /* extrapolate field using equally time spaced field un,un-1,un-2, (at
       dt intervals) to time n+t at order Ord */
    void IncNavierStokes::SubStepExtrapoloteField(NekDouble toff, Array< OneD, Array<OneD, NekDouble> > &ExtVel)
    {
        int npts = m_fields[0]->GetTotPoints();
        int nvel = m_velocity.num_elements();
        int i,j;
        Array<OneD, NekDouble> l(4);

        int ord = m_intSteps;
            
        // calculate Lagrange interpolants
        Vmath::Fill(4,1.0,l,1);
        
        for(i = 0; i <= ord; ++i)
        {
            for(j = 0; j <= ord; ++j)
            {
                if(i != j)
                {
                    l[i] *= (j*m_timestep+toff);
                    l[i] /= (j*m_timestep-i*m_timestep);
                }
            }
        }

        for(i = 0; i < nvel; ++i)
        {
            Vmath::Smul(npts,l[0],m_previousVelFields[i],1,ExtVel[i],1);
            
            for(j = 1; j <= ord; ++j)
            {
                Blas::Daxpy(npts,l[j],m_previousVelFields[j*nvel+i],1,
                            ExtVel[i],1);
            }
        }
    }
    

    // Evaluation -N(V) for all fields except pressure using m_velocity
    void IncNavierStokes::EvaluateAdvectionTerms(const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
                                                 Array<OneD, Array<OneD, NekDouble> > &outarray, 
                                                 Array<OneD, NekDouble> &wk)
    {
        int i;
        int nvariables = inarray.num_elements();
        int nqtot      = m_fields[0]->GetTotPoints();
        int VelDim     = m_velocity.num_elements();
        Array<OneD, Array<OneD, NekDouble> > velocity(VelDim);
        Array<OneD, NekDouble > Deriv;
        for(i = 0; i < VelDim; ++i)
        {
            velocity[i] = inarray[m_velocity[i]]; 
        }
        // Set up Derivative work space; 
        if(wk.num_elements())
        {
            ASSERTL0(wk.num_elements() >= nqtot*VelDim,"Workspace is not sufficient");            
            Deriv = wk;
        }
        else
        {
            Deriv = Array<OneD, NekDouble> (nqtot*VelDim);
        }

        m_advObject->DoAdvection(m_fields,m_nConvectiveFields, 
                                 m_velocity,inarray,outarray,m_time,Deriv);
    }
    
    //time dependent boundary conditions updating
    
    void IncNavierStokes::SetBoundaryConditions(NekDouble time)
    {
        int  nvariables = m_fields.num_elements();
	
        for (int i = 0; i < nvariables; ++i)
        {
            for(int n = 0; n < m_fields[i]->GetBndConditions().num_elements(); ++n)
            {	
                if(m_fields[i]->GetBndConditions()[n]->GetUserDefined() == SpatialDomains::eTimeDependent)
                {
                    m_fields[i]->EvaluateBoundaryConditions(time);
                }
            }
        }
    }
    
    // Decide if at a steady state if the discrerte L2 sum of the
    // coefficients is the same as the previous step to within the
    // tolerance m_steadyStateTol;
    bool IncNavierStokes::CalcSteadyState(void)
    {
        static NekDouble previousL2 = 0.0;
        bool returnval = false;

        NekDouble L2 = 0.0;
        
        // calculate L2 discrete summation 
        int ncoeffs = m_fields[0]->GetNcoeffs(); 
        
        for(int i = 0; i < m_fields.num_elements(); ++i)
        {
            L2 += Vmath::Dot(ncoeffs,m_fields[i]->GetCoeffs(),1,m_fields[i]->GetCoeffs(),1);
        }
        
        if(fabs(L2-previousL2) < ncoeffs*m_steadyStateTol)
        {
            returnval = true;
        }

        previousL2 = L2;

        return returnval;
    }

} //end of namespace

