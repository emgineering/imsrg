#ifndef ManagerBase_hh
#define ManagerBase_hh 1

#include "Parameters.hh"
#include "ReadWrite.hh"

class ManagerBase 
{
    public: 
        struct OpFromFile {
            std::string file2name,file3name,opname;
            int j,p,t,r; // J rank, parity, dTz, particle rank
        } ;

        Parameters parameters;
        ReadWrite rw;

        ManagerBase(Parameters parameters);
        int Solve();

        // TODO: Put in ReadWrite
        void TestScratch();
        void EnsureReadable(std::string filename);

        // TODO: consider delegating operator parsing
        std::vector<OpFromFile> GetOpsFromFile();

        // Drafts
        void InitializeH();
        void ParseOperators();
        void SolveHF();
        void SolveIMSRG();
        void TransformOperators();
};

#endif