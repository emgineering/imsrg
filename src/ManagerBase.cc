#include "ManagerBase.hh"
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <omp.h>
#include "IMSRG.hh"
#include "Parameters.hh"
#include "PhysicalConstants.hh"
#include "ModelSpace.hh"
#include "IMSRGSolver.hh"
#include "HFMBPT.hh"
#include "version.hh"

ManagerBase::ManagerBase(Parameters parameters) : 
  parameters(parameters), 
  method(parameters.s("method")),
  opnames(parameters.v("Operators")),
  flowfile(parameters.s("flowfile")),
  intfile(parameters.s("intfile"))
{
    // initialize ReadWrite
    std::string LECs = parameters.s("LECs");
    std::string scratch = parameters.s("scratch");
    std::string fmt3 = parameters.s("fmt3");

    rw.SetLECs_preset(LECs);
    rw.SetScratchDir(scratch);
    rw.Set3NFormat( fmt3 );
}
  

int ManagerBase::Solve() 
{
  std::string basis = parameters.s("basis");
  bool only_2b_omega = parameters.s("only_2b_omega") == "true";
  bool write_omega = parameters.s("write_omega") == "true";
  std::vector<std::string> spwf = parameters.v("SPWF");

  std::vector<Operator> ops;

  ParseParameterShorthand();
  ValidateReadWriteFiles();

  ModelSpace modelspace = ConfigureInputModelSpace();

  // For both dagger operators and single particle wave functions, it's convenient to
  // just get every orbit in the valence space. So if SPWF="valence" ,  we append all valence orbits
  if ( std::find( spwf.begin(), spwf.end(), "valence" ) != spwf.end() )
  {
    // this erase/remove idiom is needed because remove just shuffles things around rather than actually removing it.
    spwf.erase( std::remove( spwf.begin(), spwf.end(), "valence" ), std::end(spwf) );
    for ( auto v : modelspace.valence )
    {
      spwf.push_back( modelspace.Index2String(v) );
    }
  }

  opnames = OpParseUtil::ExpandOperatorShorthand(opnames, modelspace);

  Operator Hbare = SetupHamiltonian(modelspace);
  Operator& HNO = Hbare;

  // Performs basis transformations to HF, NAT if necessary, and does normal ordering
  HNO = GetHNOInBasis(basis, Hbare);

  // Log output
  imsrg_util::WriteSPWaveFunctions( spwf, *hf_p, intfile);
  if (method != "HF")
    PrintPerturbativeEstimates(HNO);
  std::cout << "done with pert stuff, method = " << method << std::endl;


  if ( method != "magnus" )
  {
    // Calculate all the desired operators. If we're using magnus, we'll do this after the flow is over
    TransformOperatorsConcurrent(HNO);
  }

  if (basis=="HF" or basis=="NAT")
  {
    std::cout << basis << " Single particle energies and wave functions:" << std::endl;
    (*hf_p).PrintSPEandWF();
    std::cout << std::endl;
  }

  // TODO: avoid early returns if possible. Maybe deal with this via subclassing instead of if statements
  // since methods seem to give significantly different behaviour in contained tasks
  if ( method == "HF" )
  {
    HNO.PrintTimes();
    return 0;
  }
  else if ( method == "FCI")
  {
    return FCIMethodStuff(HNO);
  }

  // TODO: where does this belong?
   if (only_2b_omega)
  {
    std::cout << " Restricting the Magnus operator Omega to be 2b." << std::endl;
    Commutator::SetOnly2bOmega(only_2b_omega);
  }

  TruncateModelSpace(HNO);

  if ( method == "MP3" )
  {
    HNO.PrintTimes();
    return 0;
  }

  SolveIMSRG(HNO);


  // if (IMSRG3)
  // {
  //   std::cout << "Norm of 3-body = " << imsrgsolver.GetH_s().ThreeBodyNorm() << std::endl;
  // }
  
  
  if (write_omega) WriteOmega();
  
  HNO.PrintTimes();
  return EXIT_SUCCESS;
}



void ManagerBase::ValidateReadWriteFiles()
{
  std::string inputtbme = parameters.s("2bme");
  std::string input3bme = parameters.s("3bme");
  std::string fmt2 = parameters.s("fmt2");


  // test 2bme file
  if (inputtbme != "none" and fmt2.find("oakridge")==std::string::npos and fmt2 != "schematic" )
  {
    if( not std::ifstream(inputtbme).good() )
    {
      std::cout << "trouble reading " << inputtbme << "  fmt2 = " << fmt2 << "   exiting. " << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  // test 3bme file
  if (input3bme != "none")
  {
    if( not std::ifstream(input3bme).good() )
    {
      std::cout << "trouble reading " << input3bme << " exiting. " << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  std::vector<std::string> opsfromfile = parameters.v("OperatorsFromFile");
  opsfromfile_unpacked = OpParseUtil::GetOpsFromFile(opsfromfile);

  TestScratch();
}

std::vector<OpParseUtil::OpFromFile> OpParseUtil::GetOpsFromFile(std::vector<std::string> tags)
{
  // unpack the awkward input format for reading an operator from file, and put it into a struct.
  // the format should look like OpName^j_t_p_r^/path/to/2bfile^/path/to/3bfile  if particle rank of Op is 2-body, then 3bfile is not needed.

  std::vector<OpFromFile> unpacked;

  // If we're reading in other operators, make sure those are ok too
  for (auto& tag : tags)
  {
     std::istringstream ss(tag);
     std::string opname,qnumbers,f2name,f3name="";

     OpFromFile opff;
  
     getline(ss,opname,'^');
     getline(ss,qnumbers,'^');
     getline(ss,f2name,'^');
     if ( not ss.eof() )  getline(ss,f3name,'^');
     opff.opname = opname;
     opff.file2name = f2name;
     opff.file3name = f3name;

      ss.str(qnumbers);
      ss.clear();
      std::string tmp;
      getline(ss,tmp,'_');
      std::istringstream(tmp) >> opff.j;
      getline(ss,tmp,'_');
      std::istringstream(tmp) >> opff.t;
      getline(ss,tmp,'_');
      std::istringstream(tmp) >> opff.p;
      getline(ss,tmp,'_');
      std::istringstream(tmp) >> opff.r;
      
      std::cout << "Parsed tag. opname = " << opff.opname << "  " << opff.j << " " << opff.t << " " << opff.p << " " << opff.r << "   file2 = " << opff.file2name   << "    file3 = " << opff.file3name << std::endl;

      // now make sure the files exist before we add them to the list.

//     if( not std::ifstream(f2name).good() )
     if( not std::ifstream(opff.file2name).good() )
     {
//       std::cout << "trouble reading " << f2name << " exiting. " << std::endl;
       std::cout << "trouble reading " << opff.file2name << " exiting. " << std::endl;
       exit(EXIT_FAILURE);
     }

     if ( opff.file3name != "") // is there a 3-body file too?
     {
//       getline(ss,f3name,'^');
//       if( not std::ifstream(f3name).good() )
       if( not std::ifstream(opff.file3name).good() )
       {
         std::cout << "trouble reading " << opff.file3name << " exiting. " << std::endl;
//         std::cout << "trouble reading " << f3name << " exiting. " << std::endl;
         exit(EXIT_FAILURE);
       }
     }
     // if the files look good, then add it to the list
     unpacked.push_back( opff );
  }
  return unpacked;
}

// TODO:
// - check if directly setting parameters here even works, or
// - parse these imsrg-related parameters at point of use, or
// - make all these member variables
void ManagerBase::ParseParameterShorthand()
{
  // deal with some short-hand method names
  if (method == "NSmagnus") // "No split" magnus
  {
    parameters.double_par["omega_norm_max"] = 50000;
    // omega_norm_max=50000;
    method = "magnus";
  }
  if (method.find("brueckner") != std::string::npos)
  {
    if (method=="brueckner2") 
    {
      parameters.string_par["brueckner_restart"] = "true";
      // brueckner_restart=true;
    }
    if (method=="brueckner1step")
    {
      parameters.int_par["nsteps"] = 1;
      parameters.string_par["core_generator"] = parameters.s("valence_generator");
      //  nsteps = 1;
      //  core_generator = valence_generator;
    }
    parameters.double_par["omega_norm_max"] = 500;
    parameters.string_par["use_brueckner_bch"] = "true";
    // use_brueckner_bch = true;
    // omega_norm_max=500;
    method = "magnus";
  }
}

void ManagerBase::TestScratch()
{
  std::vector<std::string> opnames = parameters.v("Operators");
  std::string scratch = rw.GetScratchDir();

  // Test whether the scratch directory exists and we can write to it.
  // This is necessary because otherwise you get garbage for transformed operators and it's
  // not obvious what went wrong.
  if ( (method == "magnus") and  ( (opnames.size() + opsfromfile_unpacked.size()) > 0 )  )
  {
    if ( scratch=="/dev/null" or scratch=="/dev/null/")
    {
      std::cout << "ERROR!!! using Magnus with scratch = " << scratch << " but you're also trying to transform some operators. Dying now. " << std::endl;
      exit(EXIT_FAILURE);
    }
    else if ( scratch != "" )
    {
      std::string testfilename = scratch + "/_this_is_a_test_delete_me";
      std::ofstream testout(testfilename);
      testout << "PASSED" << std::endl;
      testout.close();

      // now read it back.
      std::ifstream testin(testfilename);
      std::string checkpassed;
      testin >> checkpassed;
      if ( (checkpassed != "PASSED") or ( not testout.good() ) or ( not testin.good() ) )
      {
        std::cout << "ERROR in " << __FILE__ <<  " failed test write to scratch directory " << scratch << " that's bad. Dying now." << std::endl;
        exit(EXIT_FAILURE);
      }

    }
  }
}

ModelSpace ManagerBase::ConfigureInputModelSpace()
{
  std::string reference = parameters.s("reference");
  std::string valence_space = parameters.s("valence_space");
  std::string custom_valence_space = parameters.s("custom_valence_space");
  std::string physical_system = parameters.s("physical_system");
  std::string basis = parameters.s("basis");
  std::string occ_file = parameters.s("occ_file");

  int eMax = parameters.i("emax");
  int lmax = parameters.i("lmax"); // so far I only use this with atomic systems.
  int E3max = parameters.i("e3max");
  int emax_unocc = parameters.i("emax_unocc");
  int nsteps = parameters.i("nsteps");
  int lmax3 = parameters.i("lmax3");
  int targetMass = parameters.i("A");

  double hw = parameters.d("hw");

  if (custom_valence_space!="") // if a custom space is defined, the input valence_space is just used as a name
  {
    if (valence_space=="") // if no name is given, then just name it "custom"
    {
      parameters.string_par["valence_space"] = "custom";
      flowfile = parameters.DefaultFlowFile();
      intfile = parameters.DefaultIntFile();
    }
    valence_space = custom_valence_space;
  }


  modelspace = ( reference=="default" ? ModelSpace(eMax,valence_space) : ModelSpace(eMax,reference,valence_space) );

//  std::cout << __LINE__ << "  constructed modelspace " << std::endl;
  modelspace.SetE3max(E3max);
  modelspace.SetLmax(lmax);
//  std::cout << __LINE__ << "  done setting E3max and lmax " << std::endl;


  if (emax_unocc>0)
  {
    modelspace.SetEmaxUnocc(emax_unocc);
  }

  if (physical_system == "atomic")
  {
    modelspace.InitSingleSpecies(eMax, reference, valence_space);
  }

  if (occ_file != "none" and occ_file != "" )
  {
    modelspace.Init_occ_from_file(eMax,valence_space,occ_file);
  }


  if (nsteps < 0) // default to 1 step for single ref, 2 steps for valence decoupling
    nsteps = modelspace.valence.size()>0 ? 2 : 1;


  modelspace.SetHbarOmega(hw);
  if (targetMass>0)
     modelspace.SetTargetMass(targetMass);
  if (lmax3>0)
     modelspace.SetLmax3(lmax3);

  return modelspace;
}

std::vector<std::string> OpParseUtil::ExpandOperatorShorthand(std::vector<std::string> opnames, ModelSpace& modelspace)
{
  if ( std::find( opnames.begin(), opnames.end(), "rhop_all") != opnames.end() )
  {
    opnames.erase( std::remove( opnames.begin(), opnames.end(), "rhop_all"), std::end(opnames) );
    for ( double r=0.0; r<=10.0; r+=0.2 )
    {
       std::ostringstream opn;
       opn << "rhop_" << r;
       opnames.push_back( opn.str() );
    }
  }

  if ( std::find( opnames.begin(), opnames.end(), "rhon_all") != opnames.end() )
  {
    opnames.erase( std::remove( opnames.begin(), opnames.end(), "rhon_all"), std::end(opnames) );
    for ( double r=0.0; r<=10.0; r+=0.2 )
    {
       std::ostringstream opn;
       opn << "rhon_" << r;
       opnames.push_back( opn.str() );
    }
  }

  if ( std::find( opnames.begin(), opnames.end(), "DaggerHF_valence") != opnames.end() )
  {
    opnames.erase( std::remove( opnames.begin(), opnames.end(), "DaggerHF_valence"), std::end(opnames) );
    for ( auto v : modelspace.valence )
    {
      opnames.push_back( "DaggerHF_"+modelspace.Index2String(v) );
    }
    std::cout << "I found DaggerHF_valence, so I'm changing the opnames list to :" << std::endl;
    for ( auto opn : opnames ) std::cout << opn << " ,  ";
    std::cout << std::endl;
  }

  if ( std::find( opnames.begin(), opnames.end(), "DaggerAlln_valence") != opnames.end() )
  {
    opnames.erase( std::remove( opnames.begin(), opnames.end(), "DaggerAlln_valence"), std::end(opnames) );
    for ( auto v : modelspace.valence )
    {
      opnames.push_back( "DaggerAlln_"+modelspace.Index2String(v) );
    }
    std::cout << "I found DaggerAlln_valence, so I'm changing the opnames list to :" << std::endl;
    for ( auto opn : opnames ) std::cout << opn << " ,  ";
    std::cout << std::endl;
  }

  return opnames;
}


Operator ManagerBase::SetupHamiltonian(ModelSpace& modelspace)
{
  std::string inputtbme = parameters.s("2bme");
  std::string input3bme = parameters.s("3bme");
  std::string fmt2 = parameters.s("fmt2");
  std::string input3bme_type = parameters.s("3bme_type");
  std::string no2b_precision = parameters.s("no2b_precision");
  std::string physical_system = parameters.s("physical_system");

  bool goose_tank = parameters.s("goose_tank") == "true";
  bool store_3bme_pn = (parameters.s("store_3bme_pn")=="true");
  bool nucleon_mass_correction = parameters.s("nucleon_mass_correction") == "true";
  bool relativistic_correction = parameters.s("relativistic_correction") == "true";

  int atomicZ = parameters.i("atomicZ");

  double hw_trap = parameters.d("hw_trap");
  double threebody_threshold = parameters.d("threebody_threshold");
  double BetaCM = parameters.d("BetaCM");
  double hwBetaCM = parameters.d("hwBetaCM");
  double BetaISO = parameters.d("BetaISO");

  int file2e1max = parameters.i("file2e1max");
  int file2e2max = parameters.i("file2e2max");
  int file2lmax = parameters.i("file2lmax");
  int file3e1max = parameters.i("file3e1max");
  int file3e2max = parameters.i("file3e2max");
  int file3e3max = parameters.i("file3e3max");

  //  std::cout << "Making the Hamiltonian..." << std::endl;
  int particle_rank = input3bme=="none" ? 2 : 3;
  Operator Hbare = Operator(modelspace,0,0,0,particle_rank);
  Hbare.SetHermitian();


  Commutator::SetUseGooseTank(goose_tank);
  Commutator::SetThreebodyThreshold(threebody_threshold);

  std::cout << "Reading interactions..." << std::endl;


  if (inputtbme != "none")
  {
    if (fmt2 == "me2j")
      rw.ReadBareTBME_Darmstadt(inputtbme, Hbare,file2e1max,file2e2max,file2lmax);
    else if (fmt2 == "navratil" or fmt2 == "Navratil")
      rw.ReadBareTBME_Navratil(inputtbme, Hbare);
    else if (fmt2 == "oslo" )
      rw.ReadTBME_Oslo(inputtbme, Hbare);
    else if (fmt2.find("oakridge") != std::string::npos )
    { // input format should be: singleparticle.dat,vnn.dat
      size_t comma_pos = inputtbme.find_first_of(",");
      if ( fmt2.find("bin") != std::string::npos )
        rw.ReadTBME_OakRidge( inputtbme.substr(0,comma_pos),  inputtbme.substr( comma_pos+1 ), Hbare, "binary");
      else
        rw.ReadTBME_OakRidge( inputtbme.substr(0,comma_pos),  inputtbme.substr( comma_pos+1 ), Hbare, "ascii");
    }
    else if (fmt2 == "takayuki" )
      rw.ReadTwoBody_Takayuki( inputtbme, Hbare);
    else if (fmt2 == "nushellx" )
      rw.ReadNuShellX_int( Hbare, inputtbme );
    else if (fmt2 == "schematic" )
    {
      std::cout << "using schematic potential " << inputtbme << std::endl;
      if ( inputtbme == "Minnesota") Hbare += imsrg_util::MinnesotaPotential( modelspace );
    }

    std::cout << "done reading 2N" << std::endl;
  }

  // Read in the 3-body file
  if (Hbare.particle_rank >=3)
  {
    if(input3bme_type == "full")
    {
      rw.Read_Darmstadt_3body(input3bme, Hbare, file3e1max,file3e2max,file3e3max);
    }
    if(input3bme_type == "no2b")
    {

      Hbare.ThreeBody.SetMode("no2b");
      if (no2b_precision == "half")  Hbare.ThreeBody.SetMode("no2bhalf");

      Hbare.ThreeBody.ReadFile( {input3bme}, {file3e1max, file3e2max, file3e3max, file3e1max} );
      rw.File3N = input3bme;

    }
    else if(input3bme_type == "mono")
    {
      Hbare.ThreeBody.SetMode("mono");
      Hbare.ThreeBody.ReadFile( {input3bme}, {file3e1max, file3e2max, file3e3max, file3e1max} );
      rw.File3N = input3bme;
    }
    std::cout << "done reading 3N" << std::endl;
  }

  if (store_3bme_pn)
  {
    Hbare.ThreeBody.TransformToPN();
  }




  if (inputtbme == "none" and physical_system == "atomic")
  {

    using PhysConst::M_ELECTRON;
    using PhysConst::M_NUCLEON;
    int Z = (atomicZ>=0) ?  atomicZ : modelspace.GetTargetZ() ;
    Hbare -= Z*imsrg_util::VCentralCoulomb_Op(modelspace, modelspace.Lmax) * sqrt((M_ELECTRON*1e6)/M_NUCLEON ) ;
    Hbare += imsrg_util::VCoulomb_Op(modelspace, modelspace.Lmax) * sqrt((M_ELECTRON*1e6)/M_NUCLEON ) ;  // convert oscillator length from fm with nucleon mass to nm with electon mass (in eV).
    Hbare += imsrg_util::KineticEnergy_Op(modelspace); // Don't need to rescale this, because it's related to the oscillator frequency, which we input.
    Hbare /= PhysConst::HARTREE; // Convert to Hartree
  }

  if (fmt2 != "nushellx" and physical_system != "atomic" and hw_trap < 0)  // Don't need to add kinetic energy if we read a shell model interaction
  {
    Hbare += imsrg_util::Trel_Op(modelspace);
    if (Hbare.OneBody.has_nan())
    {
       std::cout << "  Looks like the Trel op is hosed from the get go. Dying." << std::endl;
       std::exit(EXIT_FAILURE);
    }
  }

  // Add an external harmonic trap
  if ( hw_trap > 0 )
  {
    Hbare += 0.5 * (PhysConst::M_NUCLEON * hw_trap * hw_trap)/(PhysConst::HBARC*PhysConst::HBARC) * imsrg_util::RSquaredOp(modelspace); 
    Hbare += imsrg_util::KineticEnergy_Op(modelspace); // use lab-frame kinetic energy
  }

  // correction to kinetic energy because M_proton != M_neutron
  if ( nucleon_mass_correction)
  {
    Hbare += imsrg_util::Trel_Masscorrection_Op(modelspace);
  }

  if ( relativistic_correction)
  {
    Hbare += imsrg_util::KineticEnergy_RelativisticCorr(modelspace);
  }




  // Add a Lawson center of mass term. If hwBetaCM is specified, use that frequency, otherwise use the basis frequency
  if (std::abs(BetaCM)>1e-3)
  {
    if (hwBetaCM < 0) hwBetaCM = modelspace.GetHbarOmega();
    std::ostringstream hcm_opname;
    hcm_opname << "HCM_" << hwBetaCM;
    Hbare += BetaCM * imsrg_util::OperatorFromString( modelspace, hcm_opname.str());
  }


  // Add an isospin term (similar in behaviour to Lawson; to correct for HF breaking isospin https://linkinghub.elsevier.com/retrieve/pii/037594747090117X )
  if (std::abs(BetaISO)>1e-3)
   {
    Hbare += BetaISO * imsrg_util::OperatorFromString( modelspace, "Iso2");
  }

  return Hbare;
}

Operator ManagerBase::GetHNOInBasis(std::string basis, Operator& Hbare)
{
  bool freeze_occupations = parameters.s("freeze_occupations")=="true";
  bool discard_no2b_from_3n = parameters.s("discard_no2b_from_3n")=="true";

  bool IMSRG3 = parameters.s("IMSRG3") == "true";
  bool discard_residual_input3N = parameters.s("discard_residual_input3N")=="true";
  std::string input3bme_type = parameters.s("3bme_type");
  bool perturbative_triples = (parameters.s("perturbative_triples")=="true");

  bool order_NAT_by_energy = (parameters.s("order_NAT_by_energy")=="true") ? true : false;
  bool use_NAT_occupations = (parameters.s("use_NAT_occupations")=="true") ? true : false;
  std::string NAT_order = parameters.s("NAT_order");

  double OccNat3Cut = parameters.d("OccNat3Cut");
  double dE3max = parameters.d("dE3max");


  std::cout << "Creating HF" << std::endl;
  HFMBPT hf(Hbare); // HFMBPT inherits from HartreeFock, so this works for HF and NAT bases.
  hf_p = &hf;

  if (not freeze_occupations )  hf.UnFreezeOccupations();
  if ( discard_no2b_from_3n) hf.DiscardNO2Bfrom3N();

  if (basis!="oscillator")
  {
    std::cout << "Solving" << std::endl;
    hf.Solve();
  }

  // decide what to keep after normal ordering
  int hno_particle_rank = 2;
  if ((IMSRG3) and (Hbare.ThreeBodyNorm() > 1e-5))  hno_particle_rank = 3;
  if (discard_residual_input3N) hno_particle_rank = 2;
  if (input3bme_type=="no2b") hno_particle_rank = 2;

  Operator& HNO = Hbare; // The reference & means we overwrite Hbare and save some memory
  if (basis == "HF" and method !="HF")
  {
    HNO = hf.GetNormalOrderedH( hno_particle_rank );
    if ((IMSRG3 or perturbative_triples) and OccNat3Cut>0 ) hf.GetNaturalOrbitals();
  }
  else if (basis == "NAT") // we want to use the natural orbital basis
  {
    // for backwards compatibility: order_NAT_by_energy overrides NAT_order
    if (order_NAT_by_energy) NAT_order = "energy";

    hf.UseNATOccupations( use_NAT_occupations );
    hf.OrderNATBy( NAT_order );

//  GetNaturalOrbitals() calls GetDensityMatrix(), which computes the 1b density matrix up to MBPT2
//  using the NO2B Hamiltonian in the HF basis, obtained with GetNormalOrderedH().
//  Then it calls DiagonalizeRho() which diagonalizes the density matrix, yielding the natural orbital basis.
    hf.GetNaturalOrbitals();
    HNO = hf.GetNormalOrderedHNAT( hno_particle_rank );

//  SRS: I'm commenting this out because this is not reasonably-expected default behavior
//    // For now, even if we use the NAT occupations, we switch back to naive occupations after the normal ordering
//    // This should be investigated in more detail.
//    if (use_NAT_occupations)
//    {
//      hf.FillLowestOrbits();
//      std::cout << "Undoing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << std::endl;
//      HNO = HNO.UndoNormalOrdering();
//      hf.UpdateReference();
//      modelspace.SetReference(modelspace.core); // change the reference
//      std::cout << "Doing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << std::endl;
//      HNO = HNO.DoNormalOrdering();
//    }

  }
  else if (basis == "oscillator")
  {
    HNO = Hbare.DoNormalOrdering();
  }

  // TODO: ask about this section
  if (perturbative_triples)
  {
    modelspace.SetdE3max(dE3max);
    modelspace.SetOccNat3Cut(OccNat3Cut);
    std::array<size_t,2> nstates = modelspace.CountThreeBodyStatesInsideCut();
    std::cout << "We will compute perturbative triples corrections" << std::endl;
    std::cout << "Truncations: dE3max = " << dE3max << "   OccNat3Cut = " << std::scientific << OccNat3Cut << "  ->  number of 3-body states kept:  " << nstates[0] << " out of " << nstates[1] << std::endl << std::fixed;
  }

  if (IMSRG3  )
  {
    modelspace.SetdE3max(dE3max);
    modelspace.SetOccNat3Cut(OccNat3Cut);
    std::array<size_t,2> nstates = modelspace.CountThreeBodyStatesInsideCut();
    std::cout << "You have chosen IMSRG3. good luck..." << std::endl;
    std::cout << "Truncations: dE3max = " << dE3max << "   OccNat3Cut = " << std::scientific << OccNat3Cut << "  ->  number of 3-body states kept:  " << nstates[0] << " out of " << nstates[1] << std::endl;

    if (hno_particle_rank<3 ) // if we're doing IMSRG3, we need a 3 body operator
    {
      Operator H3(modelspace,0,0,0,3);
      std::cout << "Constructed H3" << std::endl;
      H3.ZeroBody = HNO.ZeroBody;
      H3.OneBody = HNO.OneBody;
      H3.TwoBody = HNO.TwoBody;
      HNO = H3;
      std::cout << "Replacing HNO" << std::endl;
      std::cout << "Hbare Three Body Norm is " << Hbare.ThreeBodyNorm() << std::endl;
      HNO.ThreeBody.SwitchToPN_and_discard();
    }
  }

  // TODO: move these corrections?
  double BetaCM = parameters.d("BetaCM");
  double hwBetaCM = parameters.d("hwBetaCM");
  double BetaISO = parameters.d("BetaISO");

  HNO -= BetaCM * 1.5*hwBetaCM; // This is just the zero-body piece. The other stuff was added earlier.

  double Tz = 0.5 * (modelspace.GetAref() - 2 * modelspace.GetZref());
  double T = abs(Tz);
  HNO -= BetaISO * T * (T+1);

  std::cout << "Hbare 0b = " << std::setprecision(8) << HNO.ZeroBody << std::endl;


  return HNO;
}


void ManagerBase::PrintPerturbativeEstimates(Operator& H){

  std::cout << "Perturbative estimates of gs energy:" << std::endl;
  double EMP2 = H.GetMP2_Energy();
  double EMP2_3B = H.GetMP2_3BEnergy();
  std::cout << "EMP2 = " << EMP2 << std::endl;
  std::cout << "EMP2_3B = " << EMP2_3B << std::endl;
  std::array<double,3> Emp_3 = H.GetMP3_Energy();
  double EMP3 = Emp_3[0]+Emp_3[1]+Emp_3[2];
  std::cout << "E3_pp = " << Emp_3[0] << "  E3_hh = " << Emp_3[1] << " E3_ph = " << Emp_3[2] << "   EMP3 = " << EMP3 << std::endl;
  std::cout << "To 3rd order, E = " << H.ZeroBody + EMP2 + EMP3 + EMP2_3B << std::endl;
}

std::vector<Operator>& ManagerBase::TransformOperatorsConcurrent(Operator& HNO)
{
  std::vector<std::string> opnamesPT1 = parameters.v("OperatorsPT1");
  std::vector<std::string> opnamesRPA = parameters.v("OperatorsRPA");
  std::vector<std::string> opnamesTDA = parameters.v("OperatorsTDA");

  int file3e1max = parameters.i("file3e1max");
  int file3e2max = parameters.i("file3e2max");
  int file3e3max = parameters.i("file3e3max");

  std::string input_op_fmt = parameters.s("input_op_fmt");
  std::string basis = parameters.s("basis");

  using PhysConst::PROTON_RCH2;
  using PhysConst::NEUTRON_RCH2;
  using PhysConst::DARWIN_FOLDY;

  HFMBPT hf = *hf_p;

  std::vector<Operator> ops;

    for (auto& opname : opnames)
    {
        ops.emplace_back( imsrg_util::OperatorFromString(modelspace, opname) );
    }

    // Calculate first order perturbative correction to some operators, if that's what we asked for.
    // Strictly speaking, it doesn't make much sense to do this and then proceed with the IMSRG calculation,
    // but I'm not here to tell people what to do...
    for (auto& opnamept1 : opnamesPT1 )
    {
      ops.emplace_back( imsrg_util::FirstOrderCorr_1b( imsrg_util::OperatorFromString(modelspace,opnamept1)   , HNO ) );
      opnames.push_back( opnamept1+"PT1" );
    }
    for (auto& opnametda : opnamesTDA )
    {  // passing the argument "TDA" just sets the phhp and hpph blocks to zero in the RPA calculation
      ops.emplace_back( imsrg_util::RPA_resummed_1b( imsrg_util::OperatorFromString(modelspace,opnametda)   , HNO, "TDA" ) );
      opnames.push_back( opnametda+"TDA" );
    }
    for (auto& opnamerpa : opnamesRPA )
    {
      ops.emplace_back( imsrg_util::RPA_resummed_1b( imsrg_util::OperatorFromString(modelspace,opnamerpa)   , HNO, "RPA" ) );
      opnames.push_back( opnamerpa+"RPA" );
    }
  


    for ( auto& opff : opsfromfile_unpacked)
    {
      Operator op(modelspace, opff.j, opff.t, opff.p, opff.r );
      if (opff.r>2) op.ThreeBody.Allocate();
      if ( input_op_fmt == "navratil" )
      {
        rw.Read2bCurrent_Navratil( opff.file2name, op );
      }
      else if ( input_op_fmt == "miyagi" )
      {
        if (opff.file2name != "")
        {   
            Operator optmp = rw.ReadOperator2b_Miyagi( opff.file2name, modelspace );
            op.TwoBody = optmp.TwoBody;
        }
        if ( opff.r>2 and opff.file3name != "")  rw.Read_Darmstadt_3body( opff.file3name, op,  file3e1max,file3e2max,file3e3max);
      }
      ops.push_back( op );
      opnames.push_back( opff.opname );
    }


   if (ops.size()>0)
   {
     std::cout << "operators to transform: " << std::endl;
     for ( auto& opn : opnames ) std::cout << opn << " ";
     std::cout << std::endl;
   }

//  for (auto& op : ops)
   for (size_t i=0;i<ops.size();++i)
   {
//     std::cout << "Before transforming  " << opnames[i] << " has 3b norm " << ops[i].ThreeBodyNorm() << std::endl;
      // We don't transform a DaggerHF, because we want the a^dagger to already refer to the HF basis.
     if ((basis == "HF") and (opnames[i].find("DaggerHF") == std::string::npos)  )
     {
       ops[i] = hf.TransformToHFBasis(ops[i]);
     }
     else if ((basis == "NAT") and (opnames[i].find("DaggerHF") == std::string::npos)  )
     {
       ops[i] = hf.TransformHOToNATBasis(ops[i]);
     }
//     std::cout << "After transforming  " << opnames[i] << " has 3b norm " << ops[i].ThreeBodyNorm() << std::endl;
     ops[i] = ops[i].DoNormalOrdering();
//     std::cout << "Before normal ordering  " << opnames[i] << " has 3b norm " << ops[i].ThreeBodyNorm() << std::endl;
       std::cout << basis << " expectation value  " << opnames[i] << "  " << ops[i].ZeroBody << std::endl;
     if (method == "MP3")
     {
       double dop = ops[i].MP1_Eval( HNO );
       std::cout << "Operator 1st order correction  " << dop << "  ->  " << ops[i].ZeroBody + dop << std::endl;
     }
    if ( opnames[i] == "Rp2" )
    {
      double Rp2 = ops[i].ZeroBody;
      int Z = modelspace.GetTargetZ();
      int A = modelspace.GetTargetMass();
      std::cout << " HF point proton radius = " << sqrt( Rp2 ) << std::endl;
      std::cout << " HF charge radius = " << ( abs(Rp2)<1e-6 ? 0.0 : sqrt( Rp2 + PROTON_RCH2 + NEUTRON_RCH2*(A-Z)/Z + DARWIN_FOLDY) ) << std::endl;
    }
   }// for ops.size


   return ops;
}

// TODO: deal with this
// TODO #2: carefully trace where modelspace is passed by reference/value.
// Changed modelspace to HNO.modelspace in this method and not sure if this is valid
int ManagerBase::FCIMethodStuff(Operator& HNO)
{
  std::string valence_file_format = parameters.s("valence_file_format");

   if ( valence_file_format == "tokyo" )
   {
      HNO = HNO.UndoNormalOrdering();
      for (size_t i=0; i<ops.size();i++)
      {
         ops[i] = ops[i].UndoNormalOrdering();
      }

      (*HNO.modelspace).SetReference("vacuum");
      rw.WriteTokyo(HNO,intfile+".snt", "");
      // Haven't yet implemented FCI operators for Tokyo format. I should do this...
      for (size_t i=0; i<ops.size();i++)
      {
         if (ops[i].GetJRank()==0 and ops[i].GetTRank()==0 )
         {
           rw.WriteTokyo(ops[i], intfile + "_" + opnames[i] + ".snt","");
         }
         else
         {
          rw.WriteTensorTokyo(intfile+opnames[i]+".snt",ops[i]);
         }
      }
   }
   else // Write in NuShellX Format
   {
     // we want the 1b piece to be diagonal in the vacuum NO representation
      HNO = HNO.UndoNormalOrdering();
      double previous_zero_body = HNO.ZeroBody;
      (*HNO.modelspace).SetReference("vacuum");
      HartreeFock hfvac(HNO);
      hfvac.Solve();
  
  //    Operator Hvac = hfvac.GetNormalOrderedH();
      HNO = hfvac.GetNormalOrderedH();
      std::cout << "HNO had zero body = " << HNO.ZeroBody << "  and I add " << previous_zero_body << " to it. " << std::endl;
      HNO.ZeroBody += previous_zero_body;
  
      rw.WriteNuShellX_int(HNO,intfile+".int");
      rw.WriteNuShellX_sps(HNO,intfile+".sp");
  
      std::cout << "NO wrt vacuum. One Body term is hopfully still diagonal?" << std::endl << HNO.OneBody << std::endl;
  
      for (index_t i=0;i<ops.size();++i)
      {
        ops[i] = ops[i].UndoNormalOrdering();
        if ((ops[i].GetJRank()+ops[i].GetTRank()+ops[i].GetParity())<1)
        {
          rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int");
        }
        else
        {
          rw.WriteTensorOneBody(intfile+opnames[i]+"_1b.op",ops[i],opnames[i]);
          rw.WriteTensorTwoBody(intfile+opnames[i]+"_2b.op",ops[i],opnames[i]);
        }
      }
    }
    HNO.PrintTimes();
    return 0;
}

  // We may want to use a smaller model space for the IMSRG evolution than we used for the HF step.
  // This is most effective when using natural orbitals.
void ManagerBase::TruncateModelSpace(Operator& HNO)
{
  int eMax_imsrg = parameters.i("emax_imsrg");
  int e2Max_imsrg = parameters.i("e2max_imsrg");
  int e3Max_imsrg = parameters.i("e3max_imsrg");
  int E3max = parameters.i("e3max");
  if (e2Max_imsrg==-1 and eMax_imsrg != -1) e2Max_imsrg = 2*eMax_imsrg;
  if (e3Max_imsrg==-1 and eMax_imsrg != -1) e3Max_imsrg = std::min(E3max, 3*eMax_imsrg);

  std::string physical_system = parameters.s("physical_system");
  std::string occ_file = parameters.s("occ_file");
  std::string reference = parameters.s("reference");
  std::string valence_space = parameters.s("valence_space");


  modelspace_imsrg = modelspace;

  if ( (eMax_imsrg == -1) and (e2Max_imsrg == -1) and (e3Max_imsrg) == -1)
  {
    HNO.SetModelSpace(modelspace_imsrg);
    return;
  }
  else
  {
     std::cout << "Truncating modelspace for IMSRG calculation: emax e2max e3max  ->  " << eMax_imsrg << " " << e2Max_imsrg << " " << e3Max_imsrg << std::endl;
     modelspace_imsrg.SetEmax( eMax_imsrg);
     modelspace_imsrg.SetE2max( e2Max_imsrg);
     modelspace_imsrg.SetE3max( e3Max_imsrg);
     modelspace_imsrg.Init( eMax_imsrg, reference, valence_space);
   //  if (emax_unocc>0) modelspace_imsrg.SetEmaxUnocc(emax_unocc);
     if (physical_system == "atomic") modelspace_imsrg.InitSingleSpecies(eMax_imsrg, reference, valence_space);
     if (occ_file != "none" and occ_file != "" ) modelspace_imsrg.Init_occ_from_file(eMax_imsrg,valence_space,occ_file);
//     if (physical_system == "atomic") modelspace_imsrg.InitSingleSpecies(eMax_imsrg, eMax_imsrg, e3Max_imsrg, reference, valence_space);
//     if (occ_file != "none" and occ_file != "" ) modelspace_imsrg.Init_occ_from_file(eMax_imsrg,e2Max_imsrg,e3Max_imsrg,valence_space,occ_file);


     // If the occupations in modelspace were different from the naive filling, we want to keep those.
     std::map<index_t,double> hole_map;
     for ( auto& i_new : modelspace_imsrg.all_orbits )
     {
        Orbit& oi_new = modelspace_imsrg.GetOrbit(i_new);
        index_t i_old = modelspace.GetOrbitIndex( oi_new.n, oi_new.l, oi_new.j2, oi_new.tz2 );
        Orbit& oi_old = modelspace.GetOrbit(i_old);
        hole_map[i_new] = oi_old.occ;
     }
     modelspace_imsrg.SetReference( hole_map );


     HNO = HNO.Truncate(modelspace_imsrg);

    // After truncating, get the perturbative energies again to see how much things changed.
     PrintPerturbativeEstimates(HNO);
  }
}

void ManagerBase::SolveIMSRG(Operator& HNO)
{
  std::string denominator_partitioning = parameters.s("denominator_partitioning");

  double eta_criterion = parameters.d("eta_criterion");
  double smax = parameters.d("smax");
  double ds_0 = parameters.d("ds_0");
  double dsmax = parameters.d("dsmax");
  double denominator_delta = parameters.d("denominator_delta");
  double domega = parameters.d("domega");
  double omega_norm_max = parameters.d("omega_norm_max");
  double OccNat3Cut = parameters.d("OccNat3Cut");
  double ode_tolerance = parameters.d("ode_tolerance");
  double threebody_threshold = parameters.d("threebody_threshold");
  double dE3max = parameters.d("dE3max");

  bool only_2b_eta = (parameters.s("only_2b_eta")=="true");
  bool hunter_gatherer = parameters.s("hunter_gatherer") == "true";
  bool perturbative_triples = (parameters.s("perturbative_triples")=="true");
  bool use_brueckner_bch = parameters.s("use_brueckner_bch") == "true";
  bool IMSRG3 = parameters.s("IMSRG3") == "true";
  bool imsrg3_n7 = parameters.s("imsrg3_n7") == "true";
  bool imsrg3_at_end = parameters.s("imsrg3_at_end") == "true";
  bool brueckner_restart = false;

  int nsteps = parameters.i("nsteps");

  std::string valence_generator = parameters.s("valence_generator");
  std::string denominator_delta_orbit = parameters.s("denominator_delta_orbit");
  std::string core_generator = parameters.s("core_generator");
  std::string reference = parameters.s("reference");
  std::string valence_space = parameters.s("valence_space");

  // need to set this post-initialization if using imsrgsolver as a member variable
  imsrgsolver.SetHin(HNO);

  imsrgsolver.SetReadWrite(rw);
  imsrgsolver.SetMethod(method);
  imsrgsolver.SetDenominatorPartitioning(denominator_partitioning);
  imsrgsolver.SetEtaCriterion(eta_criterion);
  imsrgsolver.GetGenerator().SetOnly2bEta(only_2b_eta);
  imsrgsolver.max_omega_written = 500;
  imsrgsolver.SetHunterGatherer( hunter_gatherer );
  imsrgsolver.SetPerturbativeTriples(perturbative_triples);
  imsrgsolver.SetSmax(smax);
  imsrgsolver.SetFlowFile(flowfile);
  imsrgsolver.SetDs(ds_0);
  imsrgsolver.SetDsmax(dsmax);
  imsrgsolver.SetDenominatorDelta(denominator_delta);
  imsrgsolver.SetdOmega(domega);
  imsrgsolver.SetOmegaNormMax(omega_norm_max);
  imsrgsolver.SetODETolerance(ode_tolerance);
  if (denominator_delta_orbit != "none")
    imsrgsolver.SetDenominatorDeltaOrbit(denominator_delta_orbit);

  Commutator::SetUseBruecknerBCH(use_brueckner_bch);
  Commutator::SetUseIMSRG3(IMSRG3);
  Commutator::SetUseIMSRG3N7(imsrg3_n7);
  if (use_brueckner_bch)
  {
    std::cout << "Using Brueckner flavor of BCH" << std::endl;
  }
  if (IMSRG3)
  {
    std::cout << "Using IMSRG(3) commutators. This will probably be slow..." << std::endl;
  }
  if (imsrg3_n7)
  {
    std::cout << "  only including IMSRG3 commutator terms that scale up to n7" << std::endl;
  }
  if ( threebody_threshold > 1e-12 )
  {
    std::cout << "skipping IMSRG(3) commutator terms if norm of either operator is below " << threebody_threshold << std::endl;
  }


  // Here's where we need to have the operators
//  std::cout << "MADE IT TO LINE " << __LINE__ << std::endl;


  if (method == "flow" or method == "flow_RK4" )
  {
    for (auto& op : ops )  imsrgsolver.AddOperator( op );
    std::cout << " Added ops. FlowingOps.size = " << imsrgsolver.FlowingOps.size() << std::endl;
  }


  imsrgsolver.SetGenerator(core_generator);
  if (core_generator.find("imaginary")!=std::string::npos or core_generator.find("wegner")!=std::string::npos )
  {
   if (ds_0>1e-2)
   {
     ds_0 = 1e-4;
     dsmax = 1e-2;
     imsrgsolver.SetDs(ds_0);
     imsrgsolver.SetDsmax(dsmax);
   }
  }

  imsrgsolver.Solve();

  if (IMSRG3)
  {
    std::cout << "Norm of 3-body = " << imsrgsolver.GetH_s().ThreeBodyNorm() << std::endl;
  }
  if ( perturbative_triples and method=="magnus" )
  {
//    modelspace.SetdE3max(dE3max);
//    modelspace.SetOccNat3Cut(OccNat3Cut);
//    size_t nstates_kept = modelspace.CountThreeBodyStatesInsideCut();
//    std::array<size_t,2> nstates = modelspace.CountThreeBodyStatesInsideCut();
//    std::cout << "Truncations: dE3max = " << dE3max << "   OccNat3Cut = " << std::scientific << OccNat3Cut << "  ->  number of 3-body states kept:  " << nstates[0] << " out of " << nstates[1] << std::endl << std::fixed;
    double dE_triples = imsrgsolver.GetPerturbativeTriples();
    std::cout << "Perturbative triples:  " << std::setw(16) << std::setprecision(8) << dE_triples << " -> " << imsrgsolver.GetH_s().ZeroBody + dE_triples << std::endl;
  }



  if (brueckner_restart)
  {
     arma::mat newC = (*hf_p).C * arma::expmat( -imsrgsolver.GetOmega(0).OneBody  );
//     if (input3bme != "none") Hbare.SetParticleRank(3);
     HNO = (*hf_p).GetNormalOrderedH(newC);
     imsrgsolver.SetHin(HNO);
     imsrgsolver.s = 0;
     imsrgsolver.Solve();
  }

  if (nsteps > 1 and valence_space != reference) // two-step decoupling, do core first
  {
    if (method == "magnus") smax *= 2;

    imsrgsolver.SetGenerator(valence_generator);
    std::cout << "Setting generator to " << valence_generator << std::endl;
//    modelspace.ResetFirstPass();
    modelspace_imsrg.ResetFirstPass();
    if (valence_generator.find("imaginary")!=std::string::npos or valence_generator.find("wegner")!=std::string::npos)
    {
     if (ds_0>1e-2)
     {
       ds_0 = 1e-4;
       dsmax = 1e-2;
       imsrgsolver.SetDs(ds_0);
       imsrgsolver.SetDsmax(dsmax);
     }
    }
    imsrgsolver.SetSmax(smax);
    imsrgsolver.Solve();
  }


  if ( imsrg3_at_end )
  {
    if ( method.find("magnus") != std::string::npos )
    {
      std::cout << "Performing final BCH transformation at the IMSRG(3) level" << std::endl;


//      modelspace.SetdE3max(dE3max);
//      modelspace.SetOccNat3Cut(OccNat3Cut);
//      int new_E3max = std::min(modelspace.GetE3max(), int( std::ceil(3*std::max( modelspace.GetEFermi()[-1], modelspace.GetEFermi()[+1])+dE3max)));
      modelspace_imsrg.SetdE3max(dE3max);
      modelspace_imsrg.SetOccNat3Cut(OccNat3Cut);
      int new_E3max = std::min(modelspace_imsrg.GetE3max(), int( std::ceil(3*std::max( modelspace_imsrg.GetEFermi()[-1], modelspace_imsrg.GetEFermi()[+1])+dE3max)));
      std::cout << "Setting new E3max = " << new_E3max << std::endl;
//      modelspace.SetE3max(  new_E3max);
      modelspace_imsrg.SetE3max(  new_E3max);

//      Operator H3(modelspace,0,0,0,3);
      Operator H3(modelspace_imsrg,0,0,0,3);
      std::cout << "Constructed H3" << std::endl;
      H3.ZeroBody = HNO.ZeroBody;
      H3.OneBody = HNO.OneBody;
      H3.TwoBody = HNO.TwoBody;
      HNO = H3;
      std::cout << "Replacing HNO" << std::endl;
      // TODO: why Hbare here?
      // std::cout << "Hbare Three Body Norm is " << Hbare.ThreeBodyNorm() << std::endl;
      HNO.ThreeBody.SwitchToPN_and_discard();


      Commutator::SetUseIMSRG3(true);
      Commutator::SetUseIMSRG3N7(imsrg3_n7);
//      Operator H_with_3 = imsrgsolver.Transform(  *(imsrgsolver.H_0) );
      Operator H_with_3 = imsrgsolver.Transform( HNO );

      // Now throw away the residual 3-body so we don't need to keep it after re-normal ordering
      H_with_3.SetNumberLegs(4);
      H_with_3.SetParticleRank(2);
      imsrgsolver.FlowingOps[0] = H_with_3;
    }
    else
    {
      std::cout << "selected imsrg3_at_end, but method != magnus, so I don't know what to do. Ignoring." << std::endl;
    }

  }
}


// If we're doing targeted/ensemble normal ordering
// we now re-normal order wrt to the core
// and do any remaining flow.
void ManagerBase::RenormalOrder(Operator& HNO)
{
  bool IMSRG3 = parameters.s("IMSRG3") == "true";

  

//  ModelSpace ms2(modelspace);
  // ModelSpace ms2(modelspace_imsrg);
  ms2 = ModelSpace(modelspace_imsrg);
  ms2.SetReference(ms2.core); // change the reference
//  if (modelspace.valence.size() > 0 )
  if (modelspace_imsrg.valence.size() > 0 )
//  if (modelspace.valence.size() > 0 or basis=="NAT")
  {
//    renormal_order = modelspace.holes.size() != modelspace.core.size();
    renormal_order = modelspace.holes.size() != modelspace_imsrg.core.size();
    if (not renormal_order)
    {
//      for (auto c : modelspace.core)
      for (auto c : modelspace_imsrg.core)
      {
//         if ( (find( modelspace.holes.begin(), modelspace.holes.end(), c) == modelspace.holes.end()) or (std::abs(1-modelspace.GetOrbit(c).occ)>1e-6))
         if ( (find( modelspace_imsrg.holes.begin(), modelspace_imsrg.holes.end(), c) == modelspace_imsrg.holes.end()) or (std::abs(1-modelspace_imsrg.GetOrbit(c).occ)>1e-6))
         {
           renormal_order = true;
           break;
         }
      }
    }
  }
  if ( renormal_order )
  {

    HNO = imsrgsolver.GetH_s();

//    int nOmega = imsrgsolver.GetOmegaSize() + imsrgsolver.GetNOmegaWritten();
//    std::cout << "Undoing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << std::endl;
    std::cout << "Undoing NO wrt A=" << modelspace_imsrg.GetAref() << " Z=" << modelspace_imsrg.GetZref() << std::endl;
    std::cout << "Before doing so, the spes are " << std::endl;
//    for ( auto i : modelspace.all_orbits ) std::cout << "  " << i << " : " << HNO.OneBody(i,i) << std::endl;
    for ( auto i : modelspace_imsrg.all_orbits ) std::cout << "  " << i << " : " << HNO.OneBody(i,i) << std::endl;
    if (IMSRG3)
    {
      std::cout << "Re-normal-ordering wrt the core. For now, we just throw away the 3N at this step." << std::endl;
      HNO.SetNumberLegs(4);
      HNO.SetParticleRank(2);
    }

    HNO = HNO.UndoNormalOrdering();
    HNO.SetModelSpace(ms2);
    std::cout << "Doing NO wrt A=" << ms2.GetAref() << " Z=" << ms2.GetZref() << "  norbits = " << ms2.GetNumberOrbits() << std::endl;
    HNO = HNO.DoNormalOrdering();

    imsrgsolver.FlowingOps[0] = HNO;

// More flowing is unnecessary, since things should stay decoupled.
//    imsrgsolver.SetHin(HNO);
//    imsrgsolver.SetEtaCriterion(1e-4);
//    imsrgsolver.Solve();
    // Change operators to the new basis, then apply the rest of the transformation
//    std::cout << "Final transformation on the operators..." << std::endl;
//    int iop = 0;
//    for (auto& op : ops)
//    {
//      std::cout << opnames[iop++] << std::endl;
//      op = op.UndoNormalOrdering();
//      op.SetModelSpace(ms2);
//      op = op.DoNormalOrdering();
//      // transform using the remaining omegas
//      op = imsrgsolver.Transform_Partial(op,nOmega);
//    }
  }
}


// FIX
  // Write the output
void ManagerBase::WriteOutput()
{
  std::string valence_file_format = parameters.s("valence_file_format");

  using PhysConst::PROTON_RCH2;
  using PhysConst::NEUTRON_RCH2;
  using PhysConst::DARWIN_FOLDY;


  // If we're doing a shell model interaction, write the
  // interaction files to disk.
//  if (modelspace.valence.size() > 0)
  if (modelspace_imsrg.valence.size() > 0)
  {
    if (valence_file_format == "antoine") // this is still being tested...
    {
      rw.WriteAntoine_int(imsrgsolver.GetH_s(),intfile+".ant");
//      rw.WriteAntoine_input(imsrgsolver.GetH_s(),intfile+".inp",modelspace.GetAref(),modelspace.GetZref());
      rw.WriteAntoine_input(imsrgsolver.GetH_s(),intfile+".inp",modelspace_imsrg.GetAref(),modelspace_imsrg.GetZref());
    }
    std::cout << "Writing files: " << intfile << std::endl;
    if (valence_file_format == "tokyo")
    {
     rw.WriteTokyo(imsrgsolver.GetH_s(),intfile+".snt", "");
    }
    else
    {
      rw.WriteNuShellX_int(imsrgsolver.GetH_s(),intfile+".int");
      rw.WriteNuShellX_sps(imsrgsolver.GetH_s(),intfile+".sp");
    }

//    if (method == "magnus" or method=="flow_RK4")
//    {
//       for (index_t i=0;i<ops.size();++i)
//       {
//          if ( ((ops[i].GetJRank()+ops[i].GetTRank()+ops[i].GetParity())<1) and (ops[i].GetNumberLegs()%2==0) )
//          {
//            if (valence_file_format == "tokyo")
//            {
//              rw.WriteTokyo(ops[i],intfile+opnames[i]+".snt", "op");
//            }
//            else
//            {
//              rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int");
//            }
//          }
//          else if ( ops[i].GetNumberLegs()%2==1) // odd number of legs -> this is a dagger operator
//          {
////            rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int"); // do this for now. later make a *.dag format.
//            rw.WriteDaggerOperator( ops[i], intfile+opnames[i]+".dag",opnames[i]);
//          }
//          else
//          {
//            if (valence_file_format == "tokyo")
//            {
//              rw.WriteTensorTokyo(intfile+opnames[i]+"_2b.snt",ops[i]);
//            }
//            else
//            {
//              rw.WriteTensorOneBody(intfile+opnames[i]+"_1b.op",ops[i],opnames[i]);
//              rw.WriteTensorTwoBody(intfile+opnames[i]+"_2b.op",ops[i],opnames[i]);
//            }
//          }
//       }
//    }
  }
  else // single ref. just print the zero body pieces out. (maybe check if its magnus?)
  {
    std::cout << "Core Energy = " << std::setprecision(6) << imsrgsolver.GetH_s().ZeroBody << std::endl;
    if ( method != "magnus")
    {
      for (index_t i=0;i<ops.size();++i)
      {
//        Operator& op = ops[i];
        Operator& op = imsrgsolver.FlowingOps[i+1]; // the first operator is the Hamiltonian
        std::cout << opnames[i] << " = " << op.ZeroBody << std::endl;
        if ( opnames[i] == "Rp2" )
        {
//           int Z = modelspace.GetTargetZ();
//           int A = modelspace.GetTargetMass();
           int Z = modelspace_imsrg.GetTargetZ();
           int A = modelspace_imsrg.GetTargetMass();
           std::cout << " IMSRG point proton radius = " << sqrt( op.ZeroBody ) << std::endl;
           std::cout << " IMSRG charge radius = " << sqrt( op.ZeroBody + PROTON_RCH2 + NEUTRON_RCH2*(A-Z)/Z + DARWIN_FOLDY) << std::endl;
        }
        if ((op.GetJRank()>0) or (op.GetTRank()>0)) // if it's a tensor, you probably want the full operator
        {
          std::cout << "Writing operator to " << intfile+opnames[i]+".op" << std::endl;
          rw.WriteOperatorHuman(op,intfile+opnames[i]+".op");
        }
      }
    }
  }

}


/////////////////////
/// Transform operators and write them

// FIX
void ManagerBase::TransformAndWriteOperatorsSequential()
{

  std::string valence_file_format = parameters.s("valence_file_format");
  std::string input_op_fmt = parameters.s("input_op_fmt");
  std::string basis = parameters.s("basis");

  int file3e1max = parameters.i("file3e1max");
  int file3e2max = parameters.i("file3e2max");
  int file3e3max = parameters.i("file3e3max");

  HFMBPT hf = *hf_p;

  int E3max = parameters.i("e3max");
  int eMax_imsrg = parameters.i("emax_imsrg");
  int e2Max_imsrg = parameters.i("e2max_imsrg");
  int e3Max_imsrg = parameters.i("e3max_imsrg");
  if (e2Max_imsrg==-1 and eMax_imsrg != -1) e2Max_imsrg = 2*eMax_imsrg;
  if (e3Max_imsrg==-1 and eMax_imsrg != -1) e3Max_imsrg = std::min(E3max, 3*eMax_imsrg);



  if (method == "magnus")
  {

    /// if method is magnus, we didn't do this already. So we need to unpack any operators from file.

    for ( auto& opff : opsfromfile_unpacked)
    {
      opnames.push_back( opff.opname + "_FROMFILE");
    }

    int count_from_file =0;


    if (opnames.size()>0) std::cout << "transforming operators" << std::endl;

    for (size_t i=0;i<opnames.size();++i)
    {
      auto opname = opnames[i];
      std::cout << i << ": " << opname << " " << std::endl;

      Operator op;

      if ( opname.find("_FROMFILE") != std::string::npos)
      {
        OpParseUtil::OpFromFile& opff = opsfromfile_unpacked[count_from_file];
         std::cout << "reading " << opff.opname << " with " << opff.j << " " << opff.t << " " << opff.p << " " << opff.r << "  from file " << opff.file2name << std::endl;
        op = Operator(modelspace, opff.j, opff.t, opff.p, opff.r );
        if (opff.r>2) op.ThreeBody.Allocate();
        if ( input_op_fmt == "navratil" )
        {
          rw.Read2bCurrent_Navratil( opff.file2name, op );
        }
        else if ( input_op_fmt == "miyagi" )
        {
          if (opff.file2name != "")
          {   
              Operator optmp = rw.ReadOperator2b_Miyagi( opff.file2name, modelspace );
              op.TwoBody = optmp.TwoBody;
          }
          if ( opff.r>2 and opff.file3name != "")  rw.Read_Darmstadt_3body( opff.file3name, op,  file3e1max,file3e2max,file3e3max);
        }
        count_from_file++;
        opname = opff.opname; // Get rid of the _FROMFILE bit.
      }
      else
      {
         op = imsrg_util::OperatorFromString( modelspace, opname );
      }
//      Operator op = imsrg_util::OperatorFromString( modelspace, opname );

      if ( basis == "HF")
      {
        op = hf.TransformToHFBasis(op).DoNormalOrdering();
      }
      if ( basis == "NAT")
      {
        op = hf.TransformHOToNATBasis(op).DoNormalOrdering();
      }
      std::cout << "   HF: " << op.ZeroBody << std::endl;

      if ( (eMax_imsrg != -1) or (e2Max_imsrg != -1) or (e3Max_imsrg) != -1)
      {
//     ModelSpace modelspace_imsrg = modelspace;
        std::cout << "Truncating modelspace for IMSRG calculation: emax e2max e3max  ->  " << eMax_imsrg << " " << e2Max_imsrg << " " << e3Max_imsrg << std::endl;
        op = op.Truncate(modelspace_imsrg);
      }



      op = imsrgsolver.Transform(op);

      std::cout << "Before renormal ordering Op(5,4) is " << std::setprecision(10) << op.OneBody(5,4) << std::endl;
      if (renormal_order) 
      {
        op = op.UndoNormalOrdering();
        op.SetModelSpace(ms2);
        op = op.DoNormalOrdering();
      }
//      std::cout << " (" << ops[i].ZeroBody << " ) " << std::endl;
      std::cout << "   IMSRG: " << op.ZeroBody << std::endl;
//      rw.WriteOperatorHuman(ops[i],intfile+opnames[i]+"_step2.op");
//      std::cout << "After renormal ordering Op(5,4) is " << std::setprecision(10) << op.OneBody(5,4) << std::endl;



    std::cout << "      " << op.GetJRank() << " " << op.GetTRank() << " " << op.GetParity() << "   " << op.GetNumberLegs() << std::endl;
    if ( ((op.GetJRank()+op.GetTRank()+op.GetParity())<1) and (op.GetNumberLegs()%2==0) )
    {
       std::cout << "writing scalar files " << std::endl;
      if (valence_file_format == "tokyo")
      {
        rw.WriteTokyo(op,intfile+opname+".snt", "op");
      }
      else
      {
        rw.WriteNuShellX_op(op,intfile+opname+".int");
      }
    }
    else if ( op.GetNumberLegs()%2==1) // odd number of legs -> this is a dagger operator
    {
//      rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int"); // do this for now. later make a *.dag format.
      rw.WriteDaggerOperator( op, intfile+opname+".dag",opname);
    }
    else
    {
       std::cout << "writing tensor files " << std::endl;
      if (valence_file_format == "tokyo")
      {
        rw.WriteTensorTokyo(intfile+opname+"_2b.snt",op);
      }
      else
      {
        rw.WriteTensorOneBody(intfile+opname+"_1b.op",op,opname);
        rw.WriteTensorTwoBody(intfile+opname+"_2b.op",op,opname);
      }
    }

    }// for opnames

  }// if method == "magnus"



  if (method == "flow" or method == "flow_RK4" )
  {
    for (size_t i=0;i<ops.size();++i)
    {
      auto op = imsrgsolver.GetOperator(i+1);  // the zero-th operator is the Hamiltonian
      auto opname = opnames[i];

      if (renormal_order) 
      {
        op = op.UndoNormalOrdering();
        op.SetModelSpace(ms2);
        op = op.DoNormalOrdering();
      }
//      std::cout << " (" << ops[i].ZeroBody << " ) " << std::endl;
      std::cout << "   IMSRG: " << op.ZeroBody << std::endl;
//      rw.WriteOperatorHuman(ops[i],intfile+opnames[i]+"_step2.op");



      std::cout << "      " << op.GetJRank() << " " << op.GetTRank() << " " << op.GetParity() << "   " << op.GetNumberLegs() << std::endl;
      if ( ((op.GetJRank()+op.GetTRank()+op.GetParity())<1) and (op.GetNumberLegs()%2==0) )
      {
         std::cout << "writing scalar files " << std::endl;
        if (valence_file_format == "tokyo")
        {
          rw.WriteTokyo(op,intfile+opname+".snt", "op");
        }
        else
        {
          rw.WriteNuShellX_op(op,intfile+opname+".int");
        }
      }
      else if ( op.GetNumberLegs()%2==1) // odd number of legs -> this is a dagger operator
      {
  //      rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int"); // do this for now. later make a *.dag format.
        rw.WriteDaggerOperator( op, intfile+opname+".dag",opname);
      }
      else
      {
         std::cout << "writing tensor files " << std::endl;
        if (valence_file_format == "tokyo")
        {
          rw.WriteTensorTokyo(intfile+opname+"_2b.snt",op);
        }
        else
        {
          rw.WriteTensorOneBody(intfile+opname+"_1b.op",op,opname);
          rw.WriteTensorTwoBody(intfile+opname+"_2b.op",op,opname);
        }
      }
    }
  }

}



void ManagerBase::WriteOmega()
{

// ReadWrite method
//  std::cout << "Made it here and write_omega is " << write_omega << std::endl;

  HFMBPT hf = *hf_p;
  std::string scratch = rw.GetScratchDir();

  imsrgsolver.FlushOmegaToScratch();
  for (int i=0; i < imsrgsolver.GetNOmegaWritten() ; i++)
  {
      std::ostringstream inputfile,outputfile;
      inputfile << scratch << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << i;
      outputfile << intfile << "_Omega_" << i;
      rw.CopyFile( inputfile.str(), outputfile.str() );
  }
//    rw.WriteOmega(intfile,scratch, imsrgsolver.n_omega_written);
  bool filesucess = hf.C.save(intfile+"C.mat");
  if (filesucess == false)
  {
    std::cout<<"Couldn't save HF coefficient matrix."<<std::endl;
  }
  // std::cout << "writing Omega to " << intfile << "_omega.op" << std::endl;
  // rw.WriteOperatorHuman(imsrgsolver.Omega.back(),intfile+"_omega.op");
}




