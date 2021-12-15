#ifndef ManagerBase_hh
#define ManagerBase_hh 1

#include "Parameters.hh"
#include "ReadWrite.hh"
#include "ModelSpace.hh"

class ManagerBase 
{
    public: 
        struct OpFromFile {
            std::string file2name,file3name,opname;
            int j,p,t,r; // J rank, parity, dTz, particle rank
        } ;

        Parameters parameters;
        ReadWrite rw;

        std::string method;
        std::string physical_system;

        std::string flowfile;
        std::string intfile;

        ManagerBase(Parameters parameters);
        int Solve();


    private:
        // TODO: Put in ReadWrite
        void TestScratch();
        void EnsureReadable(std::string filename);

        // TODO: consider delegating operator parsing
        std::vector<OpFromFile> GetOpsFromFile();

        ModelSpace ConfigureInputModelSpace();
        std::vector<std::string> ExpandOperatorShorthand(std::vector<std::string> opnames, ModelSpace modelspace);

        void ReadTwoBody(std::string inputtbme, std::string fmt2, Operator &Hbare);
        void ReadThreeBody(std::string input3bme, std::string input3bme_type, Operator &Hbare)


        // Drafts
        void InitializeH();
        void ParseOperators();
        void SolveHF();
        void SolveIMSRG();
        void TransformOperators();
};

#endif