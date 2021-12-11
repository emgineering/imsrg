#ifndef ManagerBase_hh
#define ManagerBase_hh 1

#include "Parameters.hh"

class ManagerBase 
{
    public: 
        Parameters parameters;

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