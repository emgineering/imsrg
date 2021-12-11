#ifndef ManagerBase_hh
#define ManagerBase_hh 1

#include "Parameters.hh"
#include "ReadWrite.hh"

class ManagerBase 
{
    public: 
        Parameters parameters;
        ReadWrite rw;

        ManagerBase(Parameters parameters);
        int Solve();

        // Drafts
        void InitializeH();
        void ParseOperators();
        void SolveHF();
        void SolveIMSRG();
        void TransformOperators();
};

#endif