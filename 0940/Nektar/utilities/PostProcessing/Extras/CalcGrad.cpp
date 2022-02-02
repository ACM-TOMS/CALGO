/**
 * This function calculate the vorticity vector starting from an .fld file.
 * It is meant to be used with solutions produced by the incompressible Navier-Stokes solver.
 * To use it with solutions coming form another solver further generalisations are required.
 */
#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList0D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList1DHomogeneous2D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    int i,j;

    if(argc != 4)
    {
        fprintf(stderr,"Usage: ./CalcGrad ifield file.xml file.fld\n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc-1, argv+1);


    int ifield = atoi(argv[argc-3]); 
    
    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-2]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(vSession);//meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // Import field file.
    string fieldfile(argv[argc-1]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    graphShPt->Import(fieldfile,fielddef,fielddata);
    bool useFFT = false;
    bool dealiasing = false;
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    int expdim  = graphShPt->GetMeshDimension();
    int nfields = fielddef[0]->m_fields.size();
    int graddim;
    if(expdim == 1)
    {
        if(fielddef[0]->m_numHomogeneousDir == 2)//3D Homogeneous 2D
        {
            graddim = 3;
        }
        else // 1D
        {
            graddim = 1;
        }
    }
    else if(expdim ==2)
    {
        if(fielddef[0]->m_numHomogeneousDir == 1)// 3D Homogeneous 1D
        {
            graddim = 3;
        }
        else //2D
        {
            graddim = 2;
        }
	
    }
    else // Full 3D
    {
        graddim = 3;
    }

    
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nfields + graddim);
    
    switch(expdim)
    {
    case 1:
        {
            ASSERTL0(fielddef[0]->m_numHomogeneousDir <= 2,"Quasi-3D approach is only set up for 1 or 2 homogeneous directions");
            
            if(fielddef[0]->m_numHomogeneousDir == 1)
            {
                MultiRegions::ExpList2DHomogeneous1DSharedPtr Exp2DH1;
                
                // Define Homogeneous expansion
                int nplanes = fielddef[0]->m_numModes[1];
                
                // choose points to be at evenly spaced points at
                const LibUtilities::PointsKey Pkey(nplanes+1,LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[1],nplanes,Pkey);
                NekDouble ly = fielddef[0]->m_homogeneousLengths[0];
                
                Exp2DH1 = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey,ly,useFFT,dealiasing,graphShPt);
                Exp[0] = Exp2DH1;

                for(i = 1; i < nfields + graddim; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>::AllocateSharedPtr(*Exp2DH1);
                }
            }
            else if(fielddef[0]->m_numHomogeneousDir == 2)
            {
                MultiRegions::ExpList3DHomogeneous2DSharedPtr Exp3DH2;
		
                // Define Homogeneous expansion
                int nylines = fielddef[0]->m_numModes[1];
                int nzlines = fielddef[0]->m_numModes[2];
		
                // choose points to be at evenly spaced points at
                const LibUtilities::PointsKey PkeyY(nylines+1,LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  BkeyY(fielddef[0]->m_basis[1],nylines,PkeyY);
                
                const LibUtilities::PointsKey PkeyZ(nzlines+1,LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  BkeyZ(fielddef[0]->m_basis[2],nzlines,PkeyZ);
                
                NekDouble ly = fielddef[0]->m_homogeneousLengths[0];
                NekDouble lz = fielddef[0]->m_homogeneousLengths[1];
		
                Exp3DH2 = MemoryManager<MultiRegions::ExpList3DHomogeneous2D>::AllocateSharedPtr(vSession,BkeyY,BkeyZ,ly,lz,useFFT,dealiasing,graphShPt);
                Exp[0] = Exp3DH2;
		
                for(i = 1; i < nfields + graddim; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous2D>::AllocateSharedPtr(*Exp3DH2);
                }
            }
            else
            {
                MultiRegions::ExpList1DSharedPtr Exp1D;
                Exp1D = MemoryManager<MultiRegions::ExpList1D>
                    ::AllocateSharedPtr(vSession,graphShPt);
                Exp[0] = Exp1D;
                for(i = 1; i < nfields + graddim; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList1D>
                        ::AllocateSharedPtr(*Exp1D);
                }
            }
        }
        break;
    case 2:
        {
            ASSERTL0(fielddef[0]->m_numHomogeneousDir <= 1,"NumHomogeneousDir is only set up for 1");

            if(fielddef[0]->m_numHomogeneousDir == 1)
            {
                MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;

                // Define Homogeneous expansion
                int nplanes = fielddef[0]->m_numModes[2];
                
                // choose points to be at evenly spaced points at
                // nplanes + 1 points
                const LibUtilities::PointsKey Pkey(nplanes+1,LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[2],nplanes,Pkey);
                NekDouble lz = fielddef[0]->m_homogeneousLengths[0];

                Exp3DH1 = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey,lz,useFFT,dealiasing,graphShPt);
                Exp[0] = Exp3DH1;

                for(i = 1; i < nfields + graddim; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(*Exp3DH1);
                }
            }
            else
            {
                MultiRegions::ExpList2DSharedPtr Exp2D;
                Exp2D = MemoryManager<MultiRegions::ExpList2D>
                    ::AllocateSharedPtr(vSession,graphShPt);
                Exp[0] =  Exp2D;
                
                for(i = 1; i < nfields + graddim; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList2D>
                        ::AllocateSharedPtr(*Exp2D);
                }
            }
        }
        break;
    case 3:
        {
            MultiRegions::ExpList3DSharedPtr Exp3D;
            Exp3D = MemoryManager<MultiRegions::ExpList3D>
                ::AllocateSharedPtr(vSession,graphShPt);
            Exp[0] =  Exp3D;

            for(i = 1; i < nfields + graddim; ++i)
            {
                Exp[i] = MemoryManager<MultiRegions::ExpList3D>
                    ::AllocateSharedPtr(*Exp3D);
            }
        }
        break;
    default:
        ASSERTL0(false,"Expansion dimension not recognised");
        break;
    }

    //----------------------------------------------
    // Copy data from field file
    for(j = 0; j < nfields; ++j)
    {
        for(unsigned int i = 0; i < fielddata.size(); ++i)
        {
            Exp[j]->ExtractDataToCoeffs(fielddef [i],
                                        fielddata[i],
                                        fielddef [i]->m_fields[j]);
        }
        Exp[j]->BwdTrans(Exp[j]->GetCoeffs(),Exp[j]->UpdatePhys());
    }
    //----------------------------------------------
    
    int nq = Exp[0]->GetNpoints();
    
    Array<OneD, NekDouble> Dy(nq);
    Array<OneD, NekDouble> Dz(nq);
    Array<OneD, NekDouble> Dx(nq);

    switch(expdim)
    {
    case 1:
        {
            if(fielddef[0]->m_numHomogeneousDir == 1)
            {
                ASSERTL0(false,"Not implemented yet");
            }
            else if(fielddef[0]->m_numHomogeneousDir == 2)
            {
                Exp[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],
                                  Exp[ifield]->GetPhys(),Dx);
                Exp[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],
                                  Exp[ifield]->GetPhys(),Dy);
                Exp[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],
                                  Exp[ifield]->GetPhys(),Dz);
                Exp[nfields  ]->FwdTrans(Dx,Exp[nfields  ]->UpdateCoeffs());
                Exp[nfields+1]->FwdTrans(Dy,Exp[nfields+1]->UpdateCoeffs());
                Exp[nfields+2]->FwdTrans(Dz,Exp[nfields+2]->UpdateCoeffs());

            }
            else 
            {
                ASSERTL0(false,"Not implemented yet");
            }
        }
        break;
    case 2:
        {
            if(fielddef[0]->m_numHomogeneousDir == 1)
            {
                Exp[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],Exp[ifield]->GetPhys(),Dx);
                Exp[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],Exp[ifield]->GetPhys(),Dy);
                Exp[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],Exp[ifield]->GetPhys(),Dz);
                Exp[nfields  ]->FwdTrans(Dx,Exp[nfields  ]->UpdateCoeffs());
                Exp[nfields+1]->FwdTrans(Dy,Exp[nfields+1]->UpdateCoeffs());
                Exp[nfields+2]->FwdTrans(Dz,Exp[nfields+2]->UpdateCoeffs());
            }
            else 
            {
                Exp[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],Exp[ifield]->GetPhys(),Dx);
                Exp[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],Exp[ifield]->GetPhys(),Dy);
                Exp[nfields  ]->FwdTrans(Dx,Exp[nfields  ]->UpdateCoeffs());
                Exp[nfields+1]->FwdTrans(Dy,Exp[nfields+1]->UpdateCoeffs());
            }
        }
        break;
    case 3:
        {
            Exp[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],Exp[ifield]->GetPhys(),Dx);
            Exp[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],Exp[ifield]->GetPhys(),Dy);
            Exp[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],Exp[ifield]->GetPhys(),Dz);
            Exp[nfields  ]->FwdTrans(Dx,Exp[nfields  ]->UpdateCoeffs());
            Exp[nfields+1]->FwdTrans(Dy,Exp[nfields+1]->UpdateCoeffs());
            Exp[nfields+2]->FwdTrans(Dz,Exp[nfields+2]->UpdateCoeffs());
        }
        break;
    default:
        {
            ASSERTL0(false,"Expansion dimension not recognised");
        }
        break;
    }
	
    //-----------------------------------------------
    // Write solution to file with additional computed fields
    string   fldfilename(argv[2]);
    string   out = fldfilename.substr(0, fldfilename.find_last_of("."));
    string   endfile("_with_grad.fld");
    out += endfile;
    std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
                                                = Exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for(j = 0; j < nfields + graddim; ++j)
    {
        for(i = 0; i < FieldDef.size(); ++i)
        {
            if (j >= nfields)
            {
                if(j == nfields)
                {
                    std::string val = FieldDef[i]->m_fields[ifield] + "_x";
                    FieldDef[i]->m_fields.push_back(val.c_str());
                }
                else if(j == nfields+1)
                {
                    std::string val = FieldDef[i]->m_fields[ifield] + "_y";
                    FieldDef[i]->m_fields.push_back(val.c_str());
                }
                else
                {
                    std::string val = FieldDef[i]->m_fields[ifield] + "_z";
                    FieldDef[i]->m_fields.push_back(val.c_str());
                }
            }
            else
            {
                FieldDef[i]->m_fields.push_back(fielddef[i]->m_fields[j]);
            }
            Exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
        }
    }
    graphShPt->Write(out, FieldDef, FieldData);
    //-----------------------------------------------

    return 0;
}

