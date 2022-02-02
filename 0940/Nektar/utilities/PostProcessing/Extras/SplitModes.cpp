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

    if(argc != 2)
    {
        fprintf(stderr,"Usage: Splitmodes fieldfile \n");
        exit(1);
    }

    //default meshgraph
    SpatialDomains::MeshGraph graph; 


    //----------------------------------------------
    // Import fieldfile2.
    string fieldfile(argv[argc-1]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    graph.Import(fieldfile,fielddef,fielddata);
    //----------------------------------------------


    ASSERTL0(fielddef[0]->m_numModes[2] > 1,"Expected Fourier field to have at least 2 modes");
    
    ASSERTL0(fielddef[0]->m_numHomogeneousDir == 1,"Expected second fld to have one homogeneous direction");
    
    int nmodes = fielddef[0]->m_numModes[2];

    // take off modes from fielddef 
    vector<unsigned int> newNumModes;
    newNumModes.push_back(fielddef[0]->m_numModes[0]);
    newNumModes.push_back(fielddef[0]->m_numModes[1]);
    vector<LibUtilities::BasisType> newBasis;
    newBasis.push_back(fielddef[0]->m_basis[0]);
    newBasis.push_back(fielddef[0]->m_basis[1]);
    for(int i = 0; i < fielddata.size(); ++i)
    {
        fielddef[i]->m_numModes = newNumModes;
        fielddef[i]->m_basis    = newBasis;
        fielddef[i]->m_numHomogeneousDir = 0;
    }


    for(int m = 0; m < nmodes; ++m)
    {
        vector<vector<NekDouble> > writedata;
        
        // outfile 
        string outfile(argv[argc-1]);
        string out = outfile.substr(0, outfile.find_last_of("."));
        char num[16] ="";
        sprintf(num,"%d",(m/2));
        
        if(m%2 == 0)
        {
            out = out + "_mode_" + num + "_real.fld";
        }
        else
        {
            out = out +  "_mode_" + num  + "_imag.fld"; 
        }
        
        for(int i = 0; i < fielddata.size(); ++i)
        {
            int j;
            int datalen = fielddata[i].size()/fielddef[i]->m_fields.size()/nmodes;

            vector<NekDouble> newdata;
        
            // Determine the number of coefficients per element
            int ncoeffs;
            switch(fielddef[i]->m_shapeType)
            {
            case SpatialDomains::eTriangle:
                ncoeffs = StdRegions::StdTriData::getNumberOfCoefficients(fielddef[i]->m_numModes[0], fielddef[i]->m_numModes[1]);
                break;
            case SpatialDomains::eQuadrilateral:
                ncoeffs = fielddef[i]->m_numModes[0]*fielddef[i]->m_numModes[1];
                break;
            default:
                ASSERTL0(false,"Shape not recognised");
                break;
            }
        
            std::vector<NekDouble>::iterator vec_iter; 
            vec_iter = fielddata[i].begin();
        
            for(int k = 0; k < fielddef[i]->m_fields.size(); ++k)
            {
                for(int n = 0; n < fielddef[i]->m_elementIDs.size(); ++n)
                {
                    // Put orginal mode in here. 
                    vec_iter += m*ncoeffs;
                    newdata.insert(newdata.end(),vec_iter,vec_iter+ncoeffs);
                    vec_iter += (nmodes-m)*ncoeffs;
                }
            }
            writedata.push_back(newdata);
        }

        //-----------------------------------------------
        // Write out datafile. 
        graph.Write(out.c_str(), fielddef, writedata);
        //-----------------------------------------------

    }
    //----------------------------------------------


    return 0;
}

