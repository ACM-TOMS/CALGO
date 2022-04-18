///////////////////////////////////////////////////////////////////////////////
//
// File ExpListHomogeneous1D.cpp
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
// Description: An ExpList which is homogeneous in 1-direction
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpListHomogeneous1D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declaration for typedefs
        ExpListHomogeneous1D::ExpListHomogeneous1D():
            ExpList(),
            m_homogeneousBasis(LibUtilities::NullBasisSharedPtr),
            m_lhom(1),
            m_homogeneous1DBlockMat(MemoryManager<Homo1DBlockMatrixMap>::AllocateSharedPtr())
        {
        }

        ExpListHomogeneous1D::ExpListHomogeneous1D(const LibUtilities::SessionReaderSharedPtr
                &pSession,const LibUtilities::BasisKey &HomoBasis, const NekDouble lhom, const bool useFFT, const bool dealiasing):
            ExpList(pSession),
            m_lhom(lhom),
            m_useFFT(useFFT),
		    m_dealiasing(dealiasing),
            m_homogeneous1DBlockMat(MemoryManager<Homo1DBlockMatrixMap>::AllocateSharedPtr())
        {
            ASSERTL2(HomoBasis != LibUtilities::NullBasisKey,"Homogeneous Basis is a null basis");
            
			m_homogeneousBasis = LibUtilities::BasisManager()[HomoBasis];
			
			m_transposition = MemoryManager<LibUtilities::Transposition>::AllocateSharedPtr(HomoBasis,m_comm->GetColumnComm());
			
			m_planes = Array<OneD,ExpListSharedPtr>(m_homogeneousBasis->GetNumPoints()/m_comm->GetColumnComm()->GetSize());
            
            if(m_useFFT)
            {
                m_FFT = LibUtilities::GetNektarFFTFactory().CreateInstance("NekFFTW", m_homogeneousBasis->GetNumPoints());
            }
			
			if(m_dealiasing)
			{
				ASSERTL0(m_comm->GetColumnComm()->GetSize() == 1,"Remove dealiasing if you want to run in parallel");
				SetPaddingBase();
			}
		}


        /**
         * @param   In          ExpListHomogeneous1D object to copy.
         */
        ExpListHomogeneous1D::ExpListHomogeneous1D(const ExpListHomogeneous1D &In):
            ExpList(In,false),
            m_homogeneousBasis(In.m_homogeneousBasis),
            m_homogeneous1DBlockMat(In.m_homogeneous1DBlockMat),
            m_lhom(In.m_lhom),
            m_useFFT(In.m_useFFT),
		    m_FFT(In.m_FFT),
		    m_dealiasing(In.m_dealiasing),
		    m_padsize(In.m_padsize),
            MatBwdPAD(In.MatBwdPAD),
		    MatFwdPAD(In.MatFwdPAD),
            m_tmpIN(In.m_tmpIN),
            m_tmpOUT(In.m_tmpOUT),
		    m_transposition(In.m_transposition)
        {
            m_planes = Array<OneD, ExpListSharedPtr>(In.m_planes.num_elements());
        }

        /**
         * Destructor
         */
        ExpListHomogeneous1D::~ExpListHomogeneous1D()
        {
        }
	
        void ExpListHomogeneous1D::v_HomogeneousFwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                                         Array<OneD, NekDouble> &outarray, 
                                                         CoeffState coeffstate,
                                                         bool Shuff,
                                                         bool UnShuff)
        {
            // Forwards trans
            Homogeneous1DTrans(inarray,outarray,true,coeffstate,Shuff,UnShuff);
        }
	
        void ExpListHomogeneous1D::v_HomogeneousBwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                                         Array<OneD, NekDouble> &outarray, 
                                                         CoeffState coeffstate,
                                                         bool Shuff,
                                                         bool UnShuff)
        {
            // Backwards trans
            Homogeneous1DTrans(inarray,outarray,false,coeffstate,Shuff,UnShuff);
        }
		
        /**
         * Dealiasing routine
         */
        void ExpListHomogeneous1D::v_DealiasedProd(const Array<OneD, NekDouble> &inarray1,
                                                   const Array<OneD, NekDouble> &inarray2,
                                                   Array<OneD, NekDouble> &outarray, 
                                                   CoeffState coeffstate)
        {
            // inarray1 = first term of the product
            // inarray2 = second term of the product
            // dealiased product stored in outarray
            
            int npoints  = outarray.num_elements(); // number of total physical points
            int nplanes  = m_planes.num_elements(); // number of planes == number of Fourier modes = number of Fourier coeff
            int npencils = npoints/nplanes;         // number of pencils = numebr of physical points per plane
            
            Array<OneD, NekDouble> V1(npoints);
            Array<OneD, NekDouble> V2(npoints);
            Array<OneD, NekDouble> V1V2(npoints);
            Array<OneD, NekDouble> ShufV1(npoints);
            Array<OneD, NekDouble> ShufV2(npoints);
            Array<OneD, NekDouble> ShufV1V2(npoints);
            
            if(m_WaveSpace)
            {
                V1 = inarray1;
                V2 = inarray2;
            }
            else 
            {
                HomogeneousFwdTrans(inarray1,V1,coeffstate);
                HomogeneousFwdTrans(inarray2,V2,coeffstate);
            }
            
            m_transposition->Transpose(V1,ShufV1,false,LibUtilities::eXYtoZ);
            m_transposition->Transpose(V2,ShufV2,false,LibUtilities::eXYtoZ);
            
            /////////////////////////////////////////////////////////////////////////////
            // Creating padded vectors for each pencil
            Array<OneD, NekDouble> PadV1_pencil_coeff(m_padsize,0.0);
            Array<OneD, NekDouble> PadV2_pencil_coeff(m_padsize,0.0);
            Array<OneD, NekDouble> PadRe_pencil_coeff(m_padsize,0.0);
            
            Array<OneD, NekDouble> PadV1_pencil_phys(m_padsize,0.0);
            Array<OneD, NekDouble> PadV2_pencil_phys(m_padsize,0.0);
            Array<OneD, NekDouble> PadRe_pencil_phys(m_padsize,0.0);
            
            NekVector<NekDouble> PadIN_V1(m_padsize,PadV1_pencil_coeff,eWrapper);
            NekVector<NekDouble> PadOUT_V1(m_padsize,PadV1_pencil_phys,eWrapper);
            
            NekVector<NekDouble> PadIN_V2(m_padsize,PadV2_pencil_coeff,eWrapper);
            NekVector<NekDouble> PadOUT_V2(m_padsize,PadV2_pencil_phys,eWrapper);
            
            NekVector<NekDouble> PadIN_Re(m_padsize,PadRe_pencil_phys,eWrapper);
            NekVector<NekDouble> PadOUT_Re(m_padsize,PadRe_pencil_coeff,eWrapper);
            
            //Looping on the pencils
            for(int i = 0 ; i< npencils ; i++)
            {
                //Copying the i-th pencil pf lenght N into a bigger
                //pencil of lenght 2N We are in Fourier space
                Vmath::Vcopy(nplanes,&(ShufV1[i*nplanes]),1,&(PadV1_pencil_coeff[0]),1);
                Vmath::Vcopy(nplanes,&(ShufV2[i*nplanes]),1,&(PadV2_pencil_coeff[0]),1);
                //Moving to physical space using the padded system
                PadOUT_V1 = (*MatBwdPAD)*PadIN_V1;
                PadOUT_V2 = (*MatBwdPAD)*PadIN_V2;
                
                //Perfroming the vectors multiplication in physical space on the padded system
                Vmath::Vmul(m_padsize,PadV1_pencil_phys,1,PadV2_pencil_phys,1,PadRe_pencil_phys,1);
                
                //Moving back the result (V1*V2)_phys in Fourier space, padded system
                PadOUT_Re = (*MatFwdPAD)*PadIN_Re;
                
                //Copying the first half of the padded pencil in the full vector (Fourier space)
                Vmath::Vcopy(nplanes,&(PadRe_pencil_coeff[0]),1,&(ShufV1V2[i*nplanes]),1);
            }
            
            if(m_WaveSpace)
            {
                m_transposition->Transpose(ShufV1V2,outarray,false,LibUtilities::eZtoXY);				
            }
            else 
            {
                m_transposition->Transpose(ShufV1V2,V1V2,false,LibUtilities::eZtoXY);
                //Moving the results in physical space for the output
                HomogeneousBwdTrans(V1V2,outarray,coeffstate);
            }
        }
	
        
        /**
         * Forward transform
         */
        void ExpListHomogeneous1D::v_FwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, CoeffState coeffstate )
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            
            for(int n = 0; n < m_planes.num_elements(); ++n)
            {
                m_planes[n]->FwdTrans(inarray+cnt, tmparray = outarray + cnt1,
                                      coeffstate);
                cnt   += m_planes[n]->GetTotPoints();
                
                cnt1  += m_planes[n]->GetNcoeffs(); // need to skip ncoeffs
            }
            if(!m_WaveSpace)
            {
                HomogeneousFwdTrans(outarray,outarray,coeffstate);
            }
        }

        /**
         * Forward transform element by element
         */
        void ExpListHomogeneous1D::v_FwdTrans_IterPerExp(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            
			//spectral element FwdTrans plane by plane
            for(int n = 0; n < m_planes.num_elements(); ++n)
            {
                m_planes[n]->FwdTrans_IterPerExp(inarray+cnt, tmparray = outarray + cnt1);

                cnt   += m_planes[n]->GetTotPoints();
                cnt1  += m_planes[n]->GetNcoeffs();
            }
            if(!m_WaveSpace)
            {
                HomogeneousFwdTrans(outarray,outarray);
            }
        }
        
        /**
         * Backward transform
         */
        void ExpListHomogeneous1D::v_BwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, CoeffState coeffstate)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            
            for(int n = 0; n < m_planes.num_elements(); ++n)
            {
                m_planes[n]->BwdTrans(inarray+cnt, tmparray = outarray + cnt1,
                                      coeffstate);
                cnt  += m_planes[n]->GetNcoeffs();
                cnt1 += m_planes[n]->GetTotPoints();
            }
            if(!m_WaveSpace)
            {
                HomogeneousBwdTrans(outarray,outarray);
            }
        }
	
        /**
         * Backward transform element by element
         */
        void ExpListHomogeneous1D::v_BwdTrans_IterPerExp(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            
            for(int n = 0; n < m_planes.num_elements(); ++n)
            {
                m_planes[n]->BwdTrans_IterPerExp(inarray+cnt, tmparray = outarray + cnt1);
                
                cnt    += m_planes[n]->GetNcoeffs();
                cnt1   += m_planes[n]->GetTotPoints();
            }
            if(!m_WaveSpace)
            {
                HomogeneousBwdTrans(outarray,outarray);
            }
        }
        
        /**
         * Inner product
         */
        void ExpListHomogeneous1D::v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, CoeffState coeffstate)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            
            for(int n = 0; n < m_planes.num_elements(); ++n)
            {
                m_planes[n]->IProductWRTBase(inarray+cnt, tmparray = outarray + cnt1,coeffstate);

                cnt1    += m_planes[n]->GetNcoeffs();
                cnt   += m_planes[n]->GetTotPoints();
            }
        }
	
        /**
         * Inner product element by element
         */
        void ExpListHomogeneous1D::v_IProductWRTBase_IterPerExp(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        { 
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            
            for(int n = 0; n < m_planes.num_elements(); ++n)
            {
                m_planes[n]->IProductWRTBase_IterPerExp(inarray+cnt, tmparray = outarray + cnt1);
		
                cnt1  += m_planes[n]->GetNcoeffs();
                cnt   += m_planes[n]->GetTotPoints();
            } 
        }
	
        /**
         * Homogeneous transform Bwd/Fwd (MVM and FFT)
         */
        void ExpListHomogeneous1D::Homogeneous1DTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, 
                                                      bool IsForwards, 
                                                      CoeffState coeffstate,
                                                      bool Shuff,
                                                      bool UnShuff)
        {
            int num_dofs;
            
            if(IsForwards)
            {
                num_dofs = inarray.num_elements();
            }
            else
            {
                num_dofs = outarray.num_elements();
            }
            
            if(m_useFFT)
            {		
                
                int num_points_per_plane = num_dofs/m_planes.num_elements();
                int num_dfts_per_proc    = num_points_per_plane/m_comm->GetColumnComm()->GetSize() + (num_points_per_plane%m_comm->GetColumnComm()->GetSize() > 0);
                
                Array<OneD, NekDouble> fft_in(num_dfts_per_proc*m_homogeneousBasis->GetNumPoints(),0.0);
                Array<OneD, NekDouble> fft_out(num_dfts_per_proc*m_homogeneousBasis->GetNumPoints(),0.0);
		
                if(Shuff)
                {
                    m_transposition->Transpose(inarray,fft_in,false,LibUtilities::eXYtoZ);
                }
                else 
                {
                    Vmath::Vcopy(num_dfts_per_proc*m_homogeneousBasis->GetNumPoints(),inarray,1,fft_in,1);
                    //fft_in = inarray;
                }
                
                if(IsForwards)
                {
                    for(int i = 0 ; i < num_dfts_per_proc ; i++)
                    {
                        m_FFT->FFTFwdTrans(m_tmpIN = fft_in + i*m_homogeneousBasis->GetNumPoints(), m_tmpOUT = fft_out + i*m_homogeneousBasis->GetNumPoints());
                    }
                }
                else 
                {
                    for(int i = 0 ; i < num_dfts_per_proc ; i++)
                    {
                        m_FFT->FFTBwdTrans(m_tmpIN = fft_in + i*m_homogeneousBasis->GetNumPoints(), m_tmpOUT = fft_out + i*m_homogeneousBasis->GetNumPoints());
                    }
                }
		
                if(UnShuff)
                {
                    m_transposition->Transpose(fft_out,outarray,false,LibUtilities::eZtoXY);
                }
                else 
                {
                    Vmath::Vcopy(num_dfts_per_proc*m_homogeneousBasis->GetNumPoints(),fft_out,1,outarray,1);
                    //outarray = fft_out;
                }
            }
            else 
            {
                DNekBlkMatSharedPtr blkmat;
		
                if(num_dofs == m_npoints) //transform phys space
                {
                    if(IsForwards)
                    {
                        blkmat = GetHomogeneous1DBlockMatrix(eForwardsPhysSpace1D);
                    }
                    else
                    {
                        blkmat = GetHomogeneous1DBlockMatrix(eBackwardsPhysSpace1D);
                    }
                }
                else
                {
                    if(IsForwards)
                    {
                        blkmat = GetHomogeneous1DBlockMatrix(eForwardsCoeffSpace1D,coeffstate);
                    }
                    else
                    {
                        blkmat = GetHomogeneous1DBlockMatrix(eBackwardsCoeffSpace1D,coeffstate);
                    }
                }
		
                int nrows = blkmat->GetRows();
                int ncols = blkmat->GetColumns();
		
                Array<OneD, NekDouble> sortedinarray(ncols,0.0);
                Array<OneD, NekDouble> sortedoutarray(nrows,0.0);
		
                if(Shuff)
                {
                    m_transposition->Transpose(inarray,sortedinarray,!IsForwards,LibUtilities::eXYtoZ);
                }
                else 
                {
                    Vmath::Vcopy(ncols,inarray,1,sortedinarray,1);
                    //sortedinarray = inarray;
                }
                
                // Create NekVectors from the given data arrays
                NekVector<NekDouble> in (ncols,sortedinarray,eWrapper);
                NekVector<NekDouble> out(nrows,sortedoutarray,eWrapper);
		
                // Perform matrix-vector multiply.
                out = (*blkmat)*in;
		
                if(UnShuff)
                {
                    m_transposition->Transpose(sortedoutarray,outarray,IsForwards,LibUtilities::eZtoXY);
                }
                else 
                {
                    Vmath::Vcopy(nrows,sortedinarray,1,outarray,1);
                    //outarray = sortedinarray;
                }
                
            }
        }

        DNekBlkMatSharedPtr ExpListHomogeneous1D::GetHomogeneous1DBlockMatrix(Homogeneous1DMatType mattype, CoeffState coeffstate) const
        {
            Homo1DBlockMatrixMap::iterator matrixIter = m_homogeneous1DBlockMat->find(mattype);
            
            if(matrixIter == m_homogeneous1DBlockMat->end())
            {
                return ((*m_homogeneous1DBlockMat)[mattype] =
                        GenHomogeneous1DBlockMatrix(mattype,coeffstate));
            }
            else
            {
                return matrixIter->second;
            }
        }


        DNekBlkMatSharedPtr ExpListHomogeneous1D::GenHomogeneous1DBlockMatrix(Homogeneous1DMatType mattype, CoeffState coeffstate) const
        {
            DNekMatSharedPtr    loc_mat;
            DNekBlkMatSharedPtr BlkMatrix;
            int n_exp = 0;
            int num_trans_per_proc = 0;
            
            
            if((mattype == eForwardsCoeffSpace1D)
               ||(mattype == eBackwardsCoeffSpace1D)) // will operate on m_coeffs
            {
                n_exp = m_planes[0]->GetNcoeffs();
            }
            else
            {
                n_exp = m_planes[0]->GetTotPoints(); // will operatore on m_phys
            }
			
            num_trans_per_proc = n_exp/m_comm->GetColumnComm()->GetSize() + (n_exp%m_comm->GetColumnComm()->GetSize() > 0);

            Array<OneD,unsigned int> nrows(num_trans_per_proc);
            Array<OneD,unsigned int> ncols(num_trans_per_proc);

            if((mattype == eForwardsCoeffSpace1D)||(mattype == eForwardsPhysSpace1D))
            {
                nrows = Array<OneD, unsigned int>(num_trans_per_proc,m_homogeneousBasis->GetNumModes());
                ncols = Array<OneD, unsigned int>(num_trans_per_proc,m_homogeneousBasis->GetNumPoints());
            }
            else
            {
                nrows = Array<OneD, unsigned int>(num_trans_per_proc,m_homogeneousBasis->GetNumPoints());
                ncols = Array<OneD, unsigned int>(num_trans_per_proc,m_homogeneousBasis->GetNumModes());
            }

            MatrixStorage blkmatStorage = eDIAGONAL;
            BlkMatrix = MemoryManager<DNekBlkMat>
                ::AllocateSharedPtr(nrows,ncols,blkmatStorage);

			//Half Mode
			if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeRe || m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeIm)
			{
				StdRegions::StdPointExp StdPoint(m_homogeneousBasis->GetBasisKey());
				
				if((mattype == eForwardsCoeffSpace1D)||(mattype == eForwardsPhysSpace1D))
				{
					StdRegions::StdMatrixKey matkey(StdRegions::eFwdTrans,
													StdPoint.DetExpansionType(),
													StdPoint);
					
					loc_mat = StdPoint.GetStdMatrix(matkey);
				}
				else
				{
					StdRegions::StdMatrixKey matkey(StdRegions::eBwdTrans,
													StdPoint.DetExpansionType(),
													StdPoint);
					
					loc_mat = StdPoint.GetStdMatrix(matkey);
				}
			}
			//other cases
			else 
			{
				StdRegions::StdSegExp StdSeg(m_homogeneousBasis->GetBasisKey());
				
				if((mattype == eForwardsCoeffSpace1D)||(mattype == eForwardsPhysSpace1D))
				{
					StdRegions::StdMatrixKey matkey(StdRegions::eFwdTrans,
													StdSeg.DetExpansionType(),
													StdSeg);
					
					loc_mat = StdSeg.GetStdMatrix(matkey);
				}
				else
				{
					StdRegions::StdMatrixKey matkey(StdRegions::eBwdTrans,
													StdSeg.DetExpansionType(),
													StdSeg);
					
					loc_mat = StdSeg.GetStdMatrix(matkey);
				}				
				
			}

            // set up array of block matrices.
            for(int i = 0; i < num_trans_per_proc; ++i)
            {
                BlkMatrix->SetBlock(i,i,loc_mat);
            }

            return BlkMatrix;
        }

        std::vector<SpatialDomains::FieldDefinitionsSharedPtr> ExpListHomogeneous1D::v_GetFieldDefinitions()
        {
            std::vector<SpatialDomains::FieldDefinitionsSharedPtr> returnval;
            
			// Set up Homogeneous length details.
            Array<OneD,LibUtilities::BasisSharedPtr> HomoBasis(1,m_homogeneousBasis);
            
			std::vector<NekDouble> HomoLen;
            HomoLen.push_back(m_lhom);
			
			std::vector<unsigned int> PlanesIDs;
			
			for(int i = 0; i < m_planes.num_elements(); i++)
			{
				PlanesIDs.push_back(m_transposition->GetPlaneID(i));
			}

            m_planes[0]->GeneralGetFieldDefinitions(returnval, 1, HomoBasis, HomoLen, PlanesIDs);
            
			return returnval;
        }

        void  ExpListHomogeneous1D::v_GetFieldDefinitions(std::vector<SpatialDomains::FieldDefinitionsSharedPtr> &fielddef)
        {
            // Set up Homogeneous length details.
            Array<OneD,LibUtilities::BasisSharedPtr> HomoBasis(1,m_homogeneousBasis);
            
			std::vector<NekDouble> HomoLen;
            HomoLen.push_back(m_lhom);
			
			std::vector<unsigned int> PlanesIDs;
			
			for(int i = 0; i < m_planes.num_elements(); i++)
			{
				PlanesIDs.push_back(m_transposition->GetPlaneID(i));
			}

             // enforce NumHomoDir == 1 by direct call
            m_planes[0]->GeneralGetFieldDefinitions(fielddef,1, HomoBasis,HomoLen,PlanesIDs);
        }


        void ExpListHomogeneous1D::v_AppendFieldData(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, Array<OneD, NekDouble> &coeffs)
        {
            int i,n;
            int ncoeffs_per_plane = m_planes[0]->GetNcoeffs();
            
            // Determine mapping from element ids to location in
            // expansion list
            map<int, int> ElmtID_to_ExpID;
            for(i = 0; i < m_planes[0]->GetExpSize(); ++i)
            {
                ElmtID_to_ExpID[(*m_exp)[i]->GetGeom()->GetGlobalID()] = i;
            }

            for(i = 0; i < fielddef->m_elementIDs.size(); ++i)
            {
                int eid     = ElmtID_to_ExpID[fielddef->m_elementIDs[i]];
                int datalen = (*m_exp)[eid]->GetNcoeffs();

                for(n = 0; n < m_planes.num_elements(); ++n)
                {
                    fielddata.insert(fielddata.end(),&coeffs[m_coeff_offset[eid]+n*ncoeffs_per_plane],&coeffs[m_coeff_offset[eid]+n*ncoeffs_per_plane]+datalen);
                }
            }
        }
		
        void ExpListHomogeneous1D::v_AppendFieldData(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata)
        {
           v_AppendFieldData(fielddef,fielddata,m_coeffs);
        }

        //Extract the data in fielddata into the m_coeff list
        void ExpListHomogeneous1D::v_ExtractDataToCoeffs(
            SpatialDomains::FieldDefinitionsSharedPtr &fielddef,
            std::vector<NekDouble>                    &fielddata,
            std::string                               &field,
            Array<OneD, NekDouble>                    &coeffs)
        {
            int i,n;
            int offset = 0;
            int nzmodes;
            int datalen = fielddata.size()/fielddef->m_fields.size();
            int ncoeffs_per_plane = m_planes[0]->GetNcoeffs();
            
            // Build map of plane IDs lying on this process.
            std::map<int,int> homoZids;
            for (i = 0; i < m_planes.num_elements(); ++i)
            {
                homoZids[m_transposition->GetPlaneID(i)] = i;
            }
            
            for(i = 0; i < fielddef->m_basis.size(); ++i)
            {
                if(fielddef->m_basis[i] == m_homogeneousBasis->GetBasisType())
                {
                    nzmodes = fielddef->m_homogeneousZIDs.size();
                    break;
                }
            }
            ASSERTL1(i != fielddef->m_basis.size(),"Failed to determine number of Homogeneous modes");
            
            // Find data location according to field definition
            for(i = 0; i < fielddef->m_fields.size(); ++i)
            {
                if(fielddef->m_fields[i] == field)
                {
                    break;
                }
                offset += datalen;
            }
            ASSERTL0(i != fielddef->m_fields.size(),
                     "Field " + field + " not found in data file");
            
            // Determine mapping from element ids to location in expansion list.
            map<int, int> ElmtID_to_ExpID;
            for(i = 0; i < m_planes[0]->GetExpSize(); ++i)
            {
                ElmtID_to_ExpID[(*m_exp)[i]->GetGeom()->GetGlobalID()] = i;
            }

            int modes_offset = 0;
            int planes_offset = 0;
            Array<OneD, NekDouble> coeff_tmp;
            
            for(i = 0; i < fielddef->m_elementIDs.size(); ++i)
            {
                int eid = ElmtID_to_ExpID[fielddef->m_elementIDs[i]];
                int datalen = (*m_exp)[eid]->CalcNumberOfCoefficients(
                    fielddef->m_numModes,modes_offset);
                
                if(fielddef->m_uniOrder == true) // reset modes_offset to zero
                {
                    modes_offset = 0;
                }
                
                for(n = 0; n < nzmodes; ++n, offset += datalen)
                {
                    std::map<int,int>::iterator it = homoZids.find(
                        fielddef->m_homogeneousZIDs[n]);
                    
                    // Check to make sure this mode number lies in this field.
                    if (it == homoZids.end())
                    {
                        continue;
                    }
                    
                    planes_offset = it->second;
                    if(datalen == (*m_exp)[eid]->GetNcoeffs())
                    {
                        Vmath::Vcopy(datalen,&fielddata[offset],1,&coeffs[m_coeff_offset[eid]+planes_offset*ncoeffs_per_plane],1);
                    }
                    else // unpack data to new order
                    {
                        (*m_exp)[eid]->ExtractDataToCoeffs(fielddata, offset, fielddef->m_numModes,modes_offset,coeff_tmp = coeffs + m_coeff_offset[eid] + planes_offset*ncoeffs_per_plane);
                    }
                }
            }
        }
		
        //Extract the data in fielddata into the m_coeff list (for 2D files into 3D cases)
        void ExpListHomogeneous1D::v_ExtractDataToCoeffs(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, std::string &field, bool BaseFlow3D)
        {
            int i,n;
            int offset = 0;
            int nzmodes = m_homogeneousBasis->GetNumModes();
            int datalen = fielddata.size()/fielddef->m_fields.size();
            int ncoeffs_per_plane = m_planes[0]->GetNcoeffs();
			
            // Find data location according to field definition
            for(i = 0; i < fielddef->m_fields.size(); ++i)
            {
                if(fielddef->m_fields[i] == field)
                {
                    break;
                }
                offset += datalen;
            }
			
            ASSERTL0(i!= fielddef->m_fields.size(),"Field not found in data file");
			
            // Determine mapping from element ids to location in
            // expansion list
            map<int, int> ElmtID_to_ExpID;
            for(i = 0; i < m_planes[0]->GetExpSize(); ++i)
            {
                ElmtID_to_ExpID[(*m_exp)[i]->GetGeom()->GetGlobalID()] = i;
            }
			
            Vmath::Vcopy(datalen,&fielddata[offset],1,&m_coeffs[0],1);
        }
		

        /**
         * Write Tecplot Files Header
         * @param   outfile Output file name.
         * @param   var                 variables names
         */
        void ExpListHomogeneous1D::v_WriteTecplotHeader(std::ofstream &outfile, std::string var)
        {
            if(GetExp(0)->GetCoordim() == 1)
            {
                outfile << "Variables = x, y";
            }
            else
            {
                outfile << "Variables = x, y, z";
            }
            outfile << ", "<< var << std::endl << std::endl;
        }

        /**
         * Write Tecplot Files Field
         * @param   outfile    Output file name.
         * @param   expansion  Expansion that is considered
         */
        void ExpListHomogeneous1D::v_WriteTecplotField(std::ofstream &outfile, int expansion)
        {
            int npoints_per_plane = m_planes[0]->GetTotPoints();

            for(int n = 0; n < m_planes.num_elements(); ++n)
            {
                (*m_exp)[expansion]->SetPhys(m_phys+m_phys_offset[expansion]+
                                             n*npoints_per_plane);
                (*m_exp)[expansion]->WriteTecplotField(outfile);
            }
        }

        void ExpListHomogeneous1D::v_WriteVtkPieceData(std::ofstream &outfile, int expansion,
                                        std::string var)
        {
            int i;
            int nq = (*m_exp)[expansion]->GetTotPoints();
            int npoints_per_plane = m_planes[0]->GetTotPoints();

            // printing the fields of that zone
            outfile << "        <DataArray type=\"Float32\" Name=\""
                    << var << "\">" << endl;
            outfile << "          ";
            for (int n = 0; n < m_planes.num_elements(); ++n)
            {
                const Array<OneD, NekDouble> phys = m_phys + m_phys_offset[expansion] + n*npoints_per_plane;
                for(i = 0; i < nq; ++i)
                {
                    outfile << (fabs(phys[i]) < NekConstants::kNekZeroTol ? 0 : phys[i]) << " ";
                }
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
        }
		
        void ExpListHomogeneous1D::v_PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                                               Array<OneD, NekDouble> &out_d0,
                                               Array<OneD, NekDouble> &out_d1, 
                                               Array<OneD, NekDouble> &out_d2)
        {
            int nT_pts = inarray.num_elements();          //number of total points = n. of Fourier points * n. of points per plane (nT_pts)
            int nP_pts = nT_pts/m_planes.num_elements();    //number of points per plane = n of Fourier transform required (nP_pts)
            
            Array<OneD, NekDouble> temparray(nT_pts);
            Array<OneD, NekDouble> outarray(nT_pts);
            Array<OneD, NekDouble> tmp1;
            Array<OneD, NekDouble> tmp2;
            Array<OneD, NekDouble> tmp3;            
            
            for(int i = 0; i < m_planes.num_elements(); i++)
            {
                m_planes[i]->PhysDeriv(inarray + i*nP_pts ,tmp2 = out_d0 + i*nP_pts , tmp3 = out_d1 + i*nP_pts );
            }
            
            if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourier || m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierSingleMode || 
               m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeRe || m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeIm)			
            {
                if(m_WaveSpace)
                {
                    temparray = inarray;
                }
                else 
                { 
                    HomogeneousFwdTrans(inarray,temparray);
                }
                
                NekDouble sign = -1.0;
                NekDouble beta;
		
                //Half Mode
				if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeRe)
				{
					beta = sign*2*M_PI*(m_transposition->GetK(0))/m_lhom;
					
					Vmath::Smul(nP_pts,beta,temparray,1,outarray,1);
				}
				else if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeIm)
				{
					beta = -sign*2*M_PI*(m_transposition->GetK(0))/m_lhom;
					
					Vmath::Smul(nP_pts,beta,temparray,1,outarray,1);
				}
				
				//Fully complex
				else
				{
					for(int i = 0; i < m_planes.num_elements(); i++)
					{
						beta = -sign*2*M_PI*(m_transposition->GetK(i))/m_lhom;
						
						Vmath::Smul(nP_pts,beta,tmp1 = temparray + i*nP_pts,1,tmp2 = outarray + (i-int(sign))*nP_pts,1);
						
						sign = -1.0*sign;
					}
				}
		
                if(m_WaveSpace)
                {
                    out_d2 = outarray;
                }
                else 
                {
                    HomogeneousBwdTrans(outarray,out_d2);
                }
            }
            else 
            {
                ASSERTL0(m_comm->GetColumnComm()->GetSize() == 1,"Parallelisation in the homogeneous direction implemented just for Fourier basis");
		
                if(m_WaveSpace)
                {
                    
                    ASSERTL0(false,"Semi-phyisical time-stepping not implemented yet for non-Fourier basis");
                }
                else 
                {
                    StdRegions::StdSegExp StdSeg(m_homogeneousBasis->GetBasisKey());
                    
                    m_transposition->Transpose(inarray,temparray,false,LibUtilities::eXYtoZ);
                    
                    for(int i = 0; i < nP_pts; i++)
                    {
                        StdSeg.PhysDeriv(temparray + i*m_planes.num_elements(), tmp2 = outarray + i*m_planes.num_elements());
                    }
                    
                    m_transposition->Transpose(outarray,out_d2,false,LibUtilities::eZtoXY);
                    
                    Vmath::Smul(nT_pts,2.0/m_lhom,out_d2,1,out_d2,1);					
                }
            }
        }
	
        void ExpListHomogeneous1D::v_PhysDeriv(Direction edir,
                                               const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &out_d)
            
        {
            int nT_pts = inarray.num_elements();        //number of total points = n. of Fourier points * n. of points per plane (nT_pts)
            int nP_pts = nT_pts/m_planes.num_elements();  //number of points per plane = n of Fourier transform required (nP_pts)
            
            int dir= (int)edir;
            
            Array<OneD, NekDouble> temparray(nT_pts);
            Array<OneD, NekDouble> outarray(nT_pts);
            Array<OneD, NekDouble> tmp1;
            Array<OneD, NekDouble> tmp2;
            
            if (dir < 2)
            {
                for(int i=0; i < m_planes.num_elements(); i++)
                {
                    m_planes[i]->PhysDeriv(edir, inarray + i*nP_pts ,tmp2 = out_d + i*nP_pts);
                }
            }
            else
            {
                if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourier || m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierSingleMode || 
                   m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeRe || m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeIm)	
                {
                    if(m_WaveSpace)
                    {
                        temparray = inarray;
                    }
                    else 
                    { 
                        HomogeneousFwdTrans(inarray,temparray);
                    }
                    
                    NekDouble sign = -1.0;
                    NekDouble beta;
                    
                    //HalfMode
					if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeRe)
					{
						beta = 2*M_PI*(m_transposition->GetK(0))/m_lhom;
						
						Vmath::Smul(nP_pts,beta,temparray,1,outarray,1);
					}
					else if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeIm)
					{
						beta = -2*M_PI*(m_transposition->GetK(0))/m_lhom;
						
						Vmath::Smul(nP_pts,beta,temparray,1,outarray,1);
					}
					//Fully complex
					else
					{
						for(int i = 0; i < m_planes.num_elements(); i++)
						{
							beta = -sign*2*M_PI*(m_transposition->GetK(i))/m_lhom;
							
							Vmath::Smul(nP_pts,beta,tmp1 = temparray + i*nP_pts,1,tmp2 = outarray + (i-int(sign))*nP_pts,1);
							
							sign = -1.0*sign;
						}
					}
                    if(m_WaveSpace)
                    {
                        out_d = outarray;
                    }
                    else 
                    {
                        HomogeneousBwdTrans(outarray,out_d);
                    }
                }
                else 
                {
                    ASSERTL0(m_comm->GetColumnComm()->GetSize() == 1,"Parallelisation in the homogeneous direction implemented just for Fourier basis");
                    
                    if(m_WaveSpace)
                    {
                        ASSERTL0(false,"Semi-phyisical time-stepping not implemented yet for non-Fourier basis");
                    }
                    else 
                    {
                        StdRegions::StdSegExp StdSeg(m_homogeneousBasis->GetBasisKey());
                        
                        m_transposition->Transpose(inarray,temparray,false,LibUtilities::eXYtoZ);
                        
                        for(int i = 0; i < nP_pts; i++)
                        {
                            StdSeg.PhysDeriv(temparray + i*m_planes.num_elements(), tmp2 = outarray + i*m_planes.num_elements());
                        }
			
                        m_transposition->Transpose(outarray,out_d,false,LibUtilities::eZtoXY);
                        
                        Vmath::Smul(nT_pts,2.0/m_lhom,out_d,1,out_d,1);
                    }
                }
            }
        }
        
        void ExpListHomogeneous1D::PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD, NekDouble> &out_d0,
                                             Array<OneD, NekDouble> &out_d1, 
                                             Array<OneD, NekDouble> &out_d2)
            
        {
            v_PhysDeriv(inarray,out_d0,out_d1,out_d2);
        }
	
        void ExpListHomogeneous1D::PhysDeriv(Direction edir,
                                             const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD, NekDouble> &out_d)
        {
            v_PhysDeriv(edir,inarray,out_d);
        }
		
        /*
         * Setting the Padding base for dealisaing
         */
        void ExpListHomogeneous1D::SetPaddingBase(void)
        {
            NekDouble size = 1.5*m_homogeneousBasis->GetNumPoints();
            m_padsize = int(size);
            
            const LibUtilities::PointsKey Ppad(m_padsize,LibUtilities::eFourierEvenlySpaced);
            const LibUtilities::BasisKey  Bpad(LibUtilities::eFourier,m_padsize,Ppad);
            
            m_paddingBasis = LibUtilities::BasisManager()[Bpad];
            
            StdRegions::StdSegExp StdSeg(m_paddingBasis->GetBasisKey());
            
            StdRegions::StdMatrixKey matkey1(StdRegions::eFwdTrans,StdSeg.DetExpansionType(),StdSeg);
            StdRegions::StdMatrixKey matkey2(StdRegions::eBwdTrans,StdSeg.DetExpansionType(),StdSeg);
            
            MatFwdPAD = StdSeg.GetStdMatrix(matkey1);
            MatBwdPAD = StdSeg.GetStdMatrix(matkey2);
        }
	
        Array<OneD, unsigned int> ExpListHomogeneous1D::v_GetZIDs(void)
        {
            return m_transposition->GetPlanesIDs();
        }
    } //end of namespace
} //end of namespace


/**
* $Log: v $
*
**/

