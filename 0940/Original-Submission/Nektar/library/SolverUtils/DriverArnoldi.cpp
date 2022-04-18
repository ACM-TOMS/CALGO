///////////////////////////////////////////////////////////////////////////////
//
// File DriverArnoldi.cpp
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
// Description: Base Driver class for the stability solver
//
///////////////////////////////////////////////////////////////////////////////
#include <iomanip>

#include <SolverUtils/DriverArnoldi.h>

namespace Nektar
{
    namespace SolverUtils
    {
        /**
         *
         */
        DriverArnoldi::DriverArnoldi(const LibUtilities::SessionReaderSharedPtr pSession)
            : Driver(pSession)
        {
            m_session->LoadParameter("IO_InfoSteps", m_infosteps, 1);
        };


        /**
         *
         */
        DriverArnoldi::~DriverArnoldi()
        {
        };


        /**
         *
         */
        void DriverArnoldi::v_InitObject(ostream &out)
        {
            Driver::v_InitObject(out);
    	
            m_session->MatchSolverInfo("SolverType","VelocityCorrectionScheme",m_timeSteppingAlgorithm, false);

            if(m_timeSteppingAlgorithm)
            {
                m_period  = m_session->GetParameter("TimeStep")* m_session->GetParameter("NumSteps");
                m_nfields = m_equ[0]->UpdateFields().num_elements() - 1;
			
				if(m_session->DefinesSolverInfo("ModeType") && 
				   (m_session->GetSolverInfo("ModeType")=="SingleMode"|| m_session->GetSolverInfo("ModeType")=="HalfMode") )
                {
                    for(int i = 0; i < m_nfields; ++i)
                    {
                        m_equ[0]->UpdateFields()[i]->SetWaveSpace(true);
                    }
                }
			
            }
            else
            {
                m_period  = 1.0;
                ASSERTL0(m_session->DefinesFunction("BodyForce"),"A BodyForce section needs to be defined for this solver type");
                m_nfields = m_equ[0]->UpdateFields().num_elements();

                for(int i = 0; i < m_nfields; ++i)
                {
                    m_equ[0]->UpdateForces()[i]->SetWaveSpace(true);
                }

            }

		
            /*
              if(m_session->DefinesSolverInfo("SingleMode")==false)
              {
              for(int i = 0; i < m_nfields; ++i)
              {
						
              m_equ[0]->UpdateFields()[i]->SetWaveSpace(true);
              }
              }*/

            m_session->LoadParameter("kdim",  m_kdim,  16);
            m_session->LoadParameter("nvec",  m_nvec,  2);
            m_session->LoadParameter("nits",  m_nits,  500);
            m_session->LoadParameter("evtol", m_evtol, 1e-06);

            m_session->LoadParameter("realShift", m_realShift, 0.0);
            m_equ[0]->SetLambda(m_realShift);

            m_session->LoadParameter("imagShift", m_imagShift, 0.0);
        }

        void DriverArnoldi::ArnoldiSummary(std::ostream &out)
        {

            if(m_session->DefinesSolverInfo("SingleMode"))
            {
                out << "\tSingle Fourier mode    : true " << endl;
                ASSERTL0(m_session->DefinesSolverInfo("Homogeneous"),"Expected a homogeneous expansion to be defined with single mode");
            }
            else
            {
                out << "\tSingle Fourier mode    : false " << endl;
            }
            if(m_session->DefinesSolverInfo("BetaZero"))           
            {
                out << "\tBeta set to Zero       : true (overrides LHom)" << endl;
            }
            else
            {
                out << "\tBeta set to Zero       : false " << endl;
            }

            if(m_timeSteppingAlgorithm)
            {
                out << "\tEvolution operator     : " << m_session->GetSolverInfo("EvolutionOperator") << endl;
            }
            else
            {
                out << "\tShift (Real,Imag)      : " << m_realShift <<"," << m_imagShift <<  endl;
            }
            out << "\tKrylov-space dimension : " << m_kdim << endl;
            out << "\tNumber of vectors      : " << m_nvec << endl;
            out << "\tMax iterations         : " << m_nits << endl;
            out << "\tEigenvalue tolerance   : " << m_evtol << endl;
            out << "=======================================================================" << endl;
        }

        /**
         * Copy Arnoldi array to field variables which depend from
         * either the m_fields or m_forces
         */
        void DriverArnoldi::CopyArnoldiArrayToField(Array<OneD, NekDouble> &array)
        {
		
            Array<OneD, MultiRegions::ExpListSharedPtr> fields;
            if(m_timeSteppingAlgorithm)
            {
                fields = m_equ[0]->UpdateFields();
            }
            else
            {
                fields = m_equ[0]->UpdateForces();
            }

            int nq = fields[0]->GetNcoeffs();

            for (int k = 0; k < m_nfields; ++k)
            {
                Vmath::Vcopy(nq, &array[k*nq], 1, &fields[k]->UpdateCoeffs()[0], 1);
                fields[k]->SetPhysState(false);
            }
        };
	
        /**
         * Copy field variables which depend from either the m_fields
         * or m_forces array the Arnoldi array
         */
        void DriverArnoldi::CopyFieldToArnoldiArray(Array<OneD, NekDouble> &array)
        {

            Array<OneD, MultiRegions::ExpListSharedPtr> fields;
            fields = m_equ[m_nequ-1]->UpdateFields();
            for (int k = 0; k < m_nfields; ++k)
            {
                int nq = fields[0]->GetNcoeffs();
                Vmath::Vcopy(nq,  &fields[k]->GetCoeffs()[0], 1, &array[k*nq], 1);
                fields[k]->SetPhysState(false);
            
            }
        };
    
	
        /**
         * Initialisation for the transient growth
         */
        void DriverArnoldi::CopyFwdToAdj()
        {
            Array<OneD, MultiRegions::ExpListSharedPtr> fields;
        
            if(m_timeSteppingAlgorithm)
            {
                fields = m_equ[0]->UpdateFields();
                int nq = fields[0]->GetNcoeffs();
            
            
                for (int k=0 ; k < m_nfields; ++k)
                {
                    Vmath::Vcopy(nq,  &fields[k]->GetCoeffs()[0], 1,&m_equ[1]->UpdateFields()[k]->UpdateCoeffs()[0], 1);
                
                }
            }
            else
            {
                ASSERTL0(false,"Transient Growth non available for Coupled Solver");
            
            }
        };

        void DriverArnoldi::WriteFld(std::string file, Array<OneD, Array<OneD, NekDouble> > coeffs)
        {
        
            Array<OneD, std::string>  variables(m_nfields);
        
            ASSERTL1(coeffs.num_elements() >= m_nfields, "coeffs is not of the correct length");
            for(int i = 0; i < m_nfields; ++i)
            {
                variables[i] = m_equ[0]->GetVariable(i);
            }
        
            m_equ[0]->WriteFld(file,m_equ[0]->UpdateFields()[0], coeffs, variables);
        }


        void DriverArnoldi::WriteFld(std::string file, Array<OneD, NekDouble> coeffs)
        {
        
            Array<OneD, std::string>  variables(m_nfields);
            Array<OneD, Array<OneD, NekDouble> > fieldcoeffs(m_nfields);
        
            int ncoeffs = m_equ[0]->UpdateFields()[0]->GetNcoeffs();
            ASSERTL1(coeffs.num_elements() >= ncoeffs*m_nfields,"coeffs is not of sufficient size");

            for(int i = 0; i < m_nfields; ++i)
            {
                variables[i] = m_equ[0]->GetVariable(i);
                fieldcoeffs[i] = coeffs + i*ncoeffs;
            }
        
            m_equ[0]->WriteFld(file,m_equ[0]->UpdateFields()[0], fieldcoeffs, variables);
        }

        void DriverArnoldi::WriteEvs(ostream &evlout, const int i,  const NekDouble re_ev, const NekDouble im_ev, NekDouble resid)
        {
            if(m_timeSteppingAlgorithm)
            {
                NekDouble abs_ev = hypot (re_ev, im_ev);
                NekDouble ang_ev = atan2 (im_ev, re_ev);

                evlout << setw(2)  << i
                       << setw(12) << abs_ev
                       << setw(12) << ang_ev
                       << setw(12) << log (abs_ev) / m_period
                       << setw(12) << ang_ev       / m_period;
            
                if(resid != NekConstants::kNekUnsetDouble)
                {
                    evlout << setw(12) << resid;
                }
                evlout << endl;
            }
            else
            {
                NekDouble invmag = 1.0/(re_ev*re_ev + im_ev*im_ev);
            
                evlout << setw(2)  <<  i
                       << setw(14) <<  re_ev
                       << setw(14) <<  im_ev
                       << setw(14) <<  -re_ev*invmag + m_realShift
                       << setw(14) <<   im_ev*invmag;
            
                if(resid != NekConstants::kNekUnsetDouble)
                {
                    evlout << setw(12) << resid; 
                }
                evlout << endl;
            }
        }
    }
}
