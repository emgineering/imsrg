#ifndef ManagerBase_hh
#define ManagerBase_hh 1

#include "Parameters.hh"
#include "ReadWrite.hh"
#include "IMSRGSolver.hh"
#include "HFMBPT.hh"

namespace OpParseUtil
{
    struct OpFromFile {
            std::string file2name,file3name,opname;
            int j,p,t,r; // J rank, parity, dTz, particle rank
        } ;
    std::vector<OpFromFile> GetOpsFromFile(std::vector<std::string> tags);
    std::vector<std::string> ExpandOperatorShorthand(std::vector<std::string> opnames, ModelSpace& modelspace);
};

class ManagerBase 
{
    public: 
        Parameters parameters;
        ReadWrite rw;

        std::string method;
        
        std::vector<OpParseUtil::OpFromFile> opsfromfile_unpacked;
        std::vector<std::string> opnames;

        std::string flowfile;
        std::string intfile;

        ManagerBase(Parameters parameters);
        int Solve();

        void ValidateReadWriteFiles();
        void PrintPerturbativeEstimates(Operator& H);

    private:

        ModelSpace modelspace;
        ModelSpace modelspace_imsrg;
        ModelSpace ms2;
        HFMBPT* hf_p;
        IMSRGSolver imsrgsolver;
        std::vector<Operator> ops;

        bool renormal_order = false;
         

        void TestScratch();
        void ParseParameterShorthand();

        std::vector<Operator>& TransformOperatorsConcurrent(Operator& HNO);

        ModelSpace ConfigureInputModelSpace();
        
        Operator SetupHamiltonian(ModelSpace& modelspace);
        Operator GetHNOInBasis(std::string basis, Operator& Hbare);

        int FCIMethodStuff(Operator& HNO);

        void TruncateModelSpace(Operator& HNO);
        void SolveIMSRG(Operator& HNO);

        // void ReadTwoBody(std::string inputtbme, std::string fmt2, Operator &Hbare);
        // void ReadThreeBody(std::string input3bme, std::string input3bme_type, Operator &Hbare);


        // Operator ManagerBase::TransformToBasis();
        // void ComputeBasisTransform(Operator& H, std::string basis);
        // Operator GetNormalOrderedH(std::string basis);

        void RenormalOrder(Operator& HNO);
        void WriteOutput();
        void TransformAndWriteOperatorsSequential();
        void WriteOmega();

};

#endif