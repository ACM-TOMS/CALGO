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
    NekDouble scal1,scal2;

    if(argc != 6)
    {
        fprintf(stderr,"Usage: AddModeTo2DFld scal1 scal2 2Dfieldfile1 fieldfile2 outfield\n"
                "\t produces scal1*2Dfieldfiel1 + scal2*fieldfile2 in outfield\n" );
        exit(1);
    }

    scal1  = boost::lexical_cast<double>(argv[argc-5]);
    scal2  = boost::lexical_cast<double>(argv[argc-4]);

    //default meshgraph
    SpatialDomains::MeshGraph graph; 

    //----------------------------------------------
    // Import fieldfile1.
    string fieldfile1(argv[argc-3]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef1;
    vector<vector<NekDouble> > fielddata1;
    graph.Import(fieldfile1,fielddef1,fielddata1);
    //----------------------------------------------

    //----------------------------------------------
    // Import fieldfile2.
    string fieldfile2(argv[argc-2]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef2;
    vector<vector<NekDouble> > fielddata2;
    graph.Import(fieldfile2,fielddef2,fielddata2);
    //----------------------------------------------

    vector<vector<NekDouble> > combineddata;

    ASSERTL0(fielddata1.size() == fielddata2.size(),"Inner has different size");
    //----------------------------------------------
    // Add fielddata2 to fielddata1 using m_fields definition to align data. 

    for(int i = 0; i < fielddata2.size(); ++i)
    {
        ASSERTL0(fielddef2[i]->m_numHomogeneousDir == 1,"Expected second fld to have one homogeneous direction");
        ASSERTL0(fielddef2[i]->m_numModes[2] == 2,"Expected Fourier field to have 2 modes");

        int j;
        int datalen1 = fielddata1[i].size()/fielddef1[i]->m_fields.size();
        int datalen2 = fielddata2[i].size()/fielddef2[i]->m_fields.size();

        ASSERTL0(datalen1*2 == datalen2,"Data per fields is note compatible");

        
        // Determine the number of coefficients per element
        int ncoeffs;
        switch(fielddef2[i]->m_shapeType)
        {
        case SpatialDomains::eTriangle:
            ncoeffs = StdRegions::StdTriData::getNumberOfCoefficients(fielddef2[i]->m_numModes[0], fielddef2[i]->m_numModes[1]);
            break;
        case SpatialDomains::eQuadrilateral:
            ncoeffs = fielddef2[i]->m_numModes[0]*fielddef2[i]->m_numModes[1];
            break;
        default:
                ASSERTL0(false,"Shape not recognised");
            break;
        }
        
        // array for zero packing 
        Array<OneD,NekDouble> Zero(ncoeffs,0.0);
        
        
        // scale first field
        for(j = 0; j < fielddata1[i].size(); ++j)
        {
            fielddata1[i][j] *= scal1;
            
        }
        
        // scale second field
        for(j = 0; j < fielddata2[i].size(); ++j)
        {
            fielddata2[i][j] *= scal2;
            
        }

        std::vector<NekDouble>::iterator vec_iter; 

        vector<NekDouble> newdata;
        vec_iter = fielddata2[i].begin();

        for(int k = 0; k < fielddef2[i]->m_fields.size(); ++k)
        {
            
            // get location of 2D field information in order of field2 ordering
            int offset = 0;
            for(j = 0; j < fielddef1[i]->m_fields.size(); ++j)
            {
                if(fielddef1[i]->m_fields[j] == fielddef2[i]->m_fields[k])
                {
                    break;
                }
                offset  += datalen1;
            }
            
            if(j != fielddef1[i]->m_fields.size())
            {
                for(int n = 0; n < fielddef2[i]->m_elementIDs.size(); ++n)
                {
                    // Real zero component 
                    newdata.insert(newdata.end(),&(fielddata1[i][offset+n*ncoeffs]),
&(fielddata1[i][offset+n*ncoeffs]) + ncoeffs);
                    // Imaginary zero component; 
                    newdata.insert(newdata.end(),&Zero[0],&Zero[0] + ncoeffs);
                    
                    // Put orginal mode in here. 
                    newdata.insert(newdata.end(),vec_iter, vec_iter+2*ncoeffs);
                    vec_iter += 2*ncoeffs;
                }
            }
            else
            {
                
                for(int n = 0; n < fielddef2[i]->m_elementIDs.size(); ++n)
                {
                    // Real & Imag zero component 
                    newdata.insert(newdata.end(),&Zero[0],&Zero[0] + ncoeffs);
                    newdata.insert(newdata.end(),&Zero[0],&Zero[0] + ncoeffs);
                    
                    // Put orginal mode in here. 
                    newdata.insert(newdata.end(),vec_iter, vec_iter+2*ncoeffs);
                    vec_iter += 2*ncoeffs;
                }
            }
        }
        combineddata.push_back(newdata);
        fielddef2[i]->m_numModes[2] += 2;
        fielddef2[i]->m_homogeneousZIDs.push_back(2);
        fielddef2[i]->m_homogeneousZIDs.push_back(3);
        
        // check to see if any field in fielddef1[i]->m_fields is
        // not defined in fielddef2[i]->m_fields
        for(int k = 0; k < fielddef1[i]->m_fields.size(); ++k)
        {
            int offset = 0;
            for(j = 0; j < fielddef2[i]->m_fields.size(); ++j)
            {
                if(fielddef1[i]->m_fields[k] == fielddef2[i]->m_fields[j])
                {
                    break;
                }
            }
            
            if(j == fielddef2[i]->m_fields.size())
            {
                cout << "Warning: Field \'" << fielddef1[i]->m_fields[k] << "\' was not included in output " << endl;
            }

        }
        
    }
    //----------------------------------------------

    //-----------------------------------------------
    // Write out datafile. 
    graph.Write(argv[argc-1], fielddef2, combineddata);
    //-----------------------------------------------

    return 0;
}

