#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ContField3D.h>

using namespace Nektar;

int NoCaseStringCompare(const string & s1, const string& s2);

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    LibUtilities::CommSharedPtr vComm = vSession->GetComm();
    MultiRegions::ContField3DSharedPtr Exp, Fce;
    MultiRegions::ExpListSharedPtr DerExp1, DerExp2, DerExp3;
    int     i, nq,  coordim;
    Array<OneD,NekDouble>  fce;
    Array<OneD,NekDouble>  xc0,xc1,xc2;
    StdRegions::ConstFactorMap factors;
    FlagList flags;

    if(argc < 2)
    {
        fprintf(stderr,"Usage: Helmholtz3D  meshfile [solntype]\n");
        exit(1);
    }

    try
    {
        //----------------------------------------------
        // Read in mesh from input file
        SpatialDomains::MeshGraphSharedPtr graph3D = MemoryManager<SpatialDomains::MeshGraph3D>::AllocateSharedPtr(vSession);
        //----------------------------------------------

        //----------------------------------------------
        // Print summary of solution details
        flags.set(eUseGlobal, true);
        factors[StdRegions::eFactorLambda] = vSession->GetParameter("Lambda");
        const SpatialDomains::ExpansionMap &expansions = graph3D->GetExpansions();
        LibUtilities::BasisKey bkey0 = expansions.begin()->second->m_basisKeyVector[0];
        cout << "Solving 3D Helmholtz:"  << endl;
        cout << "         Lambda     : " << factors[StdRegions::eFactorLambda] << endl;
        cout << "         No. modes  : " << bkey0.GetNumModes() << endl;
        cout << endl;
        //----------------------------------------------

        //----------------------------------------------
        // Define Expansion
        Exp = MemoryManager<MultiRegions::ContField3D>
                        ::AllocateSharedPtr(vSession,graph3D,vSession->GetVariable(0));
        //----------------------------------------------

        //----------------------------------------------
        // Set up coordinates of mesh for Forcing function evaluation
        coordim = Exp->GetCoordim(0);
        nq      = Exp->GetTotPoints();

        xc0 = Array<OneD,NekDouble>(nq,0.0);
        xc1 = Array<OneD,NekDouble>(nq,0.0);
        xc2 = Array<OneD,NekDouble>(nq,0.0);

        switch(coordim)
        {
        case 3:
            Exp->GetCoords(xc0,xc1,xc2);
            break;
        default:
            ASSERTL0(false,"Coordim not valid");
            break;
        }
        //----------------------------------------------

        //----------------------------------------------
        // Define forcing function for first variable defined in file
        fce = Array<OneD,NekDouble>(nq);
        LibUtilities::EquationSharedPtr ffunc
                                    = vSession->GetFunction("Forcing", 0);

        ffunc->Evaluate(xc0, xc1, xc2, fce);

        //----------------------------------------------

        //----------------------------------------------
        // Setup expansion containing the  forcing function
        Fce = MemoryManager<MultiRegions::ContField3D>::AllocateSharedPtr(*Exp);
        Fce->SetPhys(fce);
        //----------------------------------------------

        //----------------------------------------------
        // Helmholtz solution taking physical forcing
        Exp->HelmSolve(Fce->GetPhys(), Exp->UpdateCoeffs(), flags, factors);
        //----------------------------------------------

        //----------------------------------------------
        // Backward Transform Solution to get solved values at
        Exp->BwdTrans(Exp->GetCoeffs(), Exp->UpdatePhys(), MultiRegions::eGlobal);
        //----------------------------------------------

        //----------------------------------------------
        // See if there is an exact solution, if so
        // evaluate and plot errors
        LibUtilities::EquationSharedPtr ex_sol
                                    = vSession->GetFunction("ExactSolution", 0);

        //-----------------------------------------------
        // Write solution to file
        string   out(vSession->GetSessionName() + ".fld");
        if (vComm->GetSize() > 1)
        {
            out += "." + boost::lexical_cast<string>(vComm->GetRank());
        }
        std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
                                                    = Exp->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

        Exp->GlobalToLocal(Exp->GetCoeffs(),Exp->UpdateCoeffs());
        for(i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back("u");
            Exp->AppendFieldData(FieldDef[i], FieldData[i]);
        }
        graph3D->Write(out, FieldDef, FieldData);
        //-----------------------------------------------

        if(ex_sol)
        {
            //----------------------------------------------
            // evaluate exact solution

            ex_sol->Evaluate(xc0, xc1, xc2, fce);

            //----------------------------------------------

            //--------------------------------------------
            // Calculate L_inf error
            Fce->SetPhys(fce);
            Fce->SetPhysState(true);


            //--------------------------------------------
            // Calculate errors
            NekDouble vLinfError = Exp->Linf(Fce->GetPhys());
            NekDouble vL2Error   = Exp->L2(Fce->GetPhys());
            NekDouble vH1Error   = Exp->H1(Fce->GetPhys());
            if (vComm->GetRank() == 0)
            {
                cout << "L infinity error: " << vLinfError << endl;
                cout << "L 2 error:        " << vL2Error << endl;
                cout << "H 1 error:        " << vH1Error << endl;

                cout << "Time in ExpEval:  " << ex_sol->GetTime() << endl;
            }
            //--------------------------------------------
        }
        //----------------------------------------------
    }
    catch (const std::runtime_error&)
    {
        cout << "Caught an error" << endl;
        return 1;
    }

    vComm->Finalise();

    return 0;
}




/**
 * Performs a case-insensitive string comparison (from web).
 * @param   s1          First string to compare.
 * @param   s2          Second string to compare.
 * @returns             0 if the strings match.
 */
int NoCaseStringCompare(const string & s1, const string& s2)
{
    string::const_iterator it1=s1.begin();
    string::const_iterator it2=s2.begin();

    //stop when either string's end has been reached
    while ( (it1!=s1.end()) && (it2!=s2.end()) )
    {
        if(::toupper(*it1) != ::toupper(*it2)) //letters differ?
        {
            // return -1 to indicate smaller than, 1 otherwise
            return (::toupper(*it1)  < ::toupper(*it2)) ? -1 : 1;
        }

        //proceed to the next character in each string
        ++it1;
        ++it2;
    }

    size_t size1=s1.size();
    size_t size2=s2.size();// cache lengths

    //return -1,0 or 1 according to strings' lengths
    if (size1==size2)
    {
        return 0;
    }

    return (size1 < size2) ? -1 : 1;
}

