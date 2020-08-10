/* vim:set noexpandtab tabstop=4 wrap */
#include <sstream>  // for capturing stderr
//#include <fstream>  //
//#include <iomanip>  // setprecision
//#include <cstdlib>  // for rand()
//#include <cmath>    // ?

#include "show_neutrons.hh"
#include "skroot.h"

// NOTE: the following warnings are expected and safe to ignore
// Warning in <TTree::Bronch>: Using split mode on a class: TLorentzVector with a custom Streamer

int MIN_EVENTS=10;      // must process at least this many events or we'll abort without doing analysis
int MAX_EVENTS=100;     // stop after processing this number of events
int WRITE_FREQUENCY=10; // write output events every N events.

// TODO: for calculating energy from momentum
// can we pull/replace this with TParticleTable or something
std::map<int,double> pdg_to_mass;

// some hack functions to try to add error handling to SK functions that don't have them built in
std::streambuf* start_stderr_capture();
std::string end_stderr_capture(std::streambuf* previous_buff);

int main(int argc, char* argv[]){
	// check arguments
	if(argc < 3){
		std::cerr<< "Usage: "<<basename(argv[0])<<" <output_root_name> <input skroot files...>"<< std::endl;
		exit(1);
	}
	
	// retrieve the output file name
	std::string outname(argv[1]);
	
	// parse the list of input files
	std::vector<std::string> file_list;
	for(int i=2; i<argc; i++){
		std::string fname_in = argv[i];
		file_list.push_back(fname_in);
	}
	
	// Generate a unique ID for the TreeManager we'll create.
	// How do we acquire a unique ID? Guess and hope for the best!
	int TreeManager_ID = 10;
	
	// Construct a ShowNeutrons class object.
	// We're using a class to encapsulate the variables we need to pass around
	// while allowing us to break code up into bite-size functions
	ShowNeutrons myNeutrons = ShowNeutrons(TreeManager_ID);
	myNeutrons.SetVerbosity(10);
	
	// Create a TreeManager and load it the input files
	int get_ok = myNeutrons.LoadSkrootFile(file_list);
	if(not get_ok){
		std::cerr<<"Failed to initialize input files for reading. Terminating"<<std::endl;
		return -1;
	}
	
	// Open output file
	get_ok = myNeutrons.CreateOutputFile(outname);
	if(not get_ok){
		std::cerr<<"Failed to initialize output file "<<outname<<" for writing. Terminating"<<std::endl;
		return -1;
	}
	
	// Process the data - here be the main event processing loop
	std::cout<<"Will process "<<myNeutrons.GetTotalEntries()<<" entries"<<std::endl;
	for(int entry_i=0; entry_i<std::min(MAX_EVENTS,myNeutrons.GetTotalEntries()); ++entry_i){
		get_ok = myNeutrons.ProcessEntry(entry_i);
		if(get_ok==1) break;
	}
	if((get_ok==1)||(myNeutrons.GetEntryNum()==(MAX_EVENTS-1))){
		std::cout<<"Successfully processed input files"<<std::endl;
	} else {
		// uhhh what happened (TODO better error checking in ProcessEntry)
		std::cerr<<"Error during event processing! Processed "
			 <<myNeutrons.GetEntryNum()<<"/"<<myNeutrons.GetTotalEntries()<<" events";
		if(myNeutrons.GetEntryNum()>MIN_EVENTS){
			std::cerr<<", finalizing analysis for events processed"<<std::endl;
		} else {
			std::cerr<<". Terminating"<<std::endl;
			return -1;
		}
	}
	
	// Write the output TTree
	// this is called intermittently to save events as we process them, but also
	// call it it at the end of processing to ensure everything is written out
	get_ok = myNeutrons.WriteTree();
	if(not get_ok){
		std::cerr<<"Error writing output TTree. Terminating"<<std::endl;
		return -1;
	}
	
	// Generate and add analysis histograms to the output file
	get_ok = myNeutrons.GenerateHistograms();
	if(not get_ok){
		std::cerr<<"Error generating output histograms."<<std::endl;
	}
	
	// close input skroot files, delete the TTreeManager
	skroot_close(&TreeManager_ID);
	// delete the SuperManager
	skroot_end();
	
	return 0;
}

int ShowNeutrons::LoadSkrootFile(std::vector<std::string> input_files){
	// open one or more skroot files
	// =============================
	
	// skroot_open_* invokes the singleton SuperManager class to create a TTreeManager
	// The TTreeManager handles file interaction, and will be associated with a given unique ID number.
	
	// skroot_open_ calls the TreeManager constructor with mode = 0;
	//  mode = 0: root 2 root (single or multi file input)
	// skroot_open_ arguments are a unique ID, ROOT output file path, and the length of that path
	// (which must be <1024). This calls CloneTree to copy events from an input SKROOT file
	// to the output SKROOT file (with some branches optionally omitted)
	// where does the input ROOT file come from? Who knows. TODO.
	
	// skroot_open_write_ calls the TreeManager constructor with mode =1;
	//  mode = 1: zbs 2 root (single or multi zbs file input)
	// skroot_open_write_ arguments are a unique ID, a ROOT output file path, and the length of that path
	// (which must be <1024). This will call 'InitZBS2Root(fname)' in the TreeManager constructor
	// where does the input ZBS file come from? Who knows. TODO
	
	// skroot_open_read_ calls the TreeManager constructor with mode = 2;
	//  mode = 2: read root -> no output written, just access root contents
	// arguments are just the unique ID.
	
	// The SuperManager keeps a map of ID->TreeManager objects,
	// BUT you need to determine the unique ID yourself!
	// If the ID is already allocated, skroot_open_* will print an error message to stderr
	// but there is no native way to tell programmatically that the allocation failed!!
	
	// FIXME Right now in the event of failure all it does is print to stderr
	// The right solution would be to modify the TreeManager+SuperManager classes
	// so that they return something to indicate success or failure.
	// but for now, as a quick band-aid, let's just try to catch the stderr output
	// and check if we got an error message.
	
	// set up capturing of stderr
	std::streambuf* nominal_stderr_bufferp = start_stderr_capture();
	// obtain the TreeManager, catching the error if it fails
	skroot_open_read(&TreeManager_ID);
	// release stderr and retrieve the captured message
	// (note this also sends the captured output to stderr as well)
	std::string error_message = end_stderr_capture(nominal_stderr_bufferp);
	// check for errors
	if(error_message!=""){
		std::cerr<<"failed to open TreeManager ID "<<TreeManager_ID<<", with message: "<<error_message;
		return false;
	} else {
		// we can safely continue
		std::cout<<"TreeManager opened successfully"<<std::endl;
	}
	
	// Set the list of input files. We can call skroot_set_input_file multiple times
	// (before initialize) to add multiple files to the chain of files to process.
	// this is just a wrapper around TreeManager::SetInputFile(const char* filename)
	for(int filei=0; filei<input_files.size(); ++filei){
		std::string fname_in = input_files.at(filei);
		skroot_set_input_file(&TreeManager_ID, fname_in.c_str(), fname_in.size());
	}
	
	// Add the input files to the TChain to read, and initialize branch addresses
	// This MUST NOT be called for mode 1 (ZBS2ROOT) (doing so will terminate the application)
	// For mode 0 this also calls TTree::CloneTree to make the output clone tree
	// this is just a wrapper around TreeManager::Initialize()
	skroot_initialize(&TreeManager_ID);
	
	// get total number of events in all loaded files
	skroot_get_entries_(&TreeManager_ID, &total_entries);
	
	// retrieve pointers to the SKROOT class objects
	// these will be updated with each TreeManager::GetEntry() call
	TreeManager* mgr = skroot_get_mgr(&TreeManager_ID);
	file_header = mgr->GetHEAD();   // NOTE! TreeManager has both GetHead() and GetHEAD() methods!!
	// MC truth particles
	mc_info = mgr->GetMC();
	
//	// LoweInfo - bonsai, clusfit, spallation likelihood, etc. Also some MC stuff, apparently.
//	lowe_info = mgr->GetLOWE();
//	// information about reconstruction of the primary muon
//	muon_info   = mgr->GetMU();
//	// what else is in the RunInfo as opposed to HEAD?
//	run_info    = mgr->GetRINFO();
//	// time to previous event? Relevant? (for analysis or MC?) Populated?
//	prev_t0     = mgr->GetPREVT0();
//	// triggering... what is in each?
//	hw_trig_list = mgr->GetHWTRGLIST();    // relevant for MC?
//	sw_trig_list = mgr->GetSoftTrgList();  // why do we have both a list of sw trigs...
//	sw_trig      = mgr->GetSTrigger();     // and a single sw trig...?
//	// event-wise info
//	event_header = mgr->GetHead();         // whats in here?
//	event_trailer = mgr->GetTrail();       // whats in here?
//	// hit information
//	tq_real_info = mgr->GetTQREALINFO();   // whats the difference between TQREAL...
//	tq_areal_info = mgr->GetTQAREALINFO(); // ...and TQAREAL?
//	// anything else useful???
	
	// List of object Getters in SKROOT files
	//////////////////////////////////////////////////////////////
	//  TTree*           GetTree()        { return theTree; }
	//  TTree*           GetOTree()       { return theOTree; }
	//  Header*          GetHEAD()        { return HEAD; }
	//  TQReal*          GetTQREALINFO()  { return TQREALINFO; }
	//  TQReal*          GetTQAREALINFO() { return TQAREALINFO; }
	//  LoweInfo*        GetLOWE()        { return LOWE; }
	//  AtmpdInfo*       GetATMPD()       { return ATMPD; }
	//  UpmuInfo*        GetUPMU()        { return UPMU; }
	//  MuInfo*          GetMU()          { return MU; }
	//  SLEInfo*         GetSLE()         { return SLE; }
	//  SoftTrgList*     GetSoftTrgList() { return SWTRGLIST; }
	//  MCInfo*          GetMC()          { return MC; }
	//  TClonesArray*    GetTQLIST()      { return TQLIST; }
	//  TClonesArray*    GetODTQLIST()    { return ODTQLIST; }
	//  TClonesArray*    GetHWTRGLIST()   { return HWTRGLIST; }
	//  Pedestal*        GetPEDESTAL()    { return PEDESTAL; }
	//  EventHeader*     GetHead()        { return Head; }
	//  EventTrailer*    GetTrail()       { return Trail; }
	//  SoftwareTrigger* GetSTrigger()    { return STrigger; }
	//  QBeeStatus*      GetQBSTATUS()    { return QBSTATUS; }
	//  QBeeStatus*      GetDBSTATUS()    { return DBSTATUS; }
	//  Spacer*          GetSPACER()      { return SPACER; }
	//  PrevT0*          GetPREVT0()      { return PREVT0; }
	//  MismatchedHit*   GetMIS()         { return MIS; }
	//  TClonesArray*    GetGPSList()     { return GPSList; }
	//  TClonesArray*    T2KGetGPSList()  { return T2KGPSList; }
	//  SlowControl*     GetSLWCTRL()     { return SLWCTRL; }
	//  RunInfo*         GetRINFO()       { return RINFO; }
	//  IDOD_Xtlk*       GetIDODXTLK()    { return IDODXTLK; }
	///////////////////////////////////////////////////////////////
	
	// Get one-off information from the TreeManager
	// --------------------------------------------
	
	// Calculate the time difference between first and last events
	run_duration_seconds = GetTdiffFirstLastEvent(start_time, stop_time);
	
	// get the MC version and TBA (table of top-bottom asymmetry of water transparency) version
	GetMCVersion();
	
	// OK, done.
	return 1;
}

int ShowNeutrons::GetMCVersion(){
	// Get MC and absorption table version from first MCInfo entry
	// ===========================================================
	// Get the TTree
	TTree* tree = skroot_get_tree(&TreeManager_ID);
	if(not tree){
		std::cerr<<"No tree from skroot_get_tree?!"<<std::endl;
		return 0;
	}
	// Get the MCInfo branch
	TBranch* MCbranch = tree->GetBranch("MC");
	if(not MCbranch){
		std::cerr<<"No MC branch in data tree?!"<<std::endl;
		return 0;
	}
	
	// read MCInfo for 1st event
	get_ok = MCbranch->GetEntry(0);
	if(not mc_info){
		std::cout<<"MCInfo still doesn't exist?!"<<std::endl;
		return 0;
	}
	
	// note the MC and absorption table versions
	skdetsim_version = mc_info->ivmcp;    // version of monte carlo program (SK_GEOMETRY+1000)
	tba_table_version = mc_info->ivabl;   // version of absorbtion length
	
	return true;
}

int ShowNeutrons::GetEntryNum(){
	return current_entry_num;
};
int ShowNeutrons::GetTotalEntries(){
	return total_entries;
};
int ShowNeutrons::WriteTree(){
	if(verbosity>2) std::cout<<"writing TTree"<<std::endl;
	outfile->cd();
	// TObject::Write returns the total number of bytes written to the file.
	// It returns 0 if the object cannot be written.
	int bytes = outtree->Write("",TObject::kOverwrite);
	if(bytes<=0){
		std::cerr<<"ShowNeutrons Error writing TTree!"<<std::endl;
	} else if(verbosity>2){
		std::cout<<"Wrote "<<get_ok<<" bytes"<<std::endl;
	}
	return bytes;
};
void ShowNeutrons::CloseFile(){
	outtree->ResetBranchAddresses();
	outfile->Write("*",TObject::kOverwrite);
	outfile->Close();
	delete outfile;
	outfile = nullptr;
};
void ShowNeutrons::SetVerbosity(int verb){
	verbosity=verb;
}

double ShowNeutrons::GetTdiffFirstLastEvent(time_t &first_event_time, time_t &last_event_time){
	// Get event(?) time of first and last TTree entry and convert to timestamps
	// ==========================================================================
	// Get the TTree
	TTree* tree = skroot_get_tree(&TreeManager_ID);
	// Get the HEAD branch
	TBranch* HEADbranch = tree->GetBranch("HEADER");
	
	// read HEADER for 1st event
	HEADbranch->GetEntry(0);
	
	// temporary structs to hold the time variables
	struct tm start_tm, stop_tm;
	
	// get time and date
	start_tm.tm_year = file_header->ndaysk[0];
	start_tm.tm_mon  = file_header->ndaysk[1] - 1; // * see below
	start_tm.tm_mday = file_header->ndaysk[2];
	start_tm.tm_hour = file_header->ntimsk[0];
	start_tm.tm_min  = file_header->ntimsk[1];
	start_tm.tm_sec  = file_header->ntimsk[2];
	// build timestamp (seconds since unix epoch)
	first_event_time = mktime(&start_tm);
	
	/* seems like SK file_header stores months 1-12 and days of month 1-31,
	   but c++ time struct stores months 0-11 and days of month 1-31 */
	
	// read only HEADER branch (last event)
	HEADbranch->GetEntry(total_entries-1);
	
	// get time and date
	stop_tm.tm_year = file_header->ndaysk[0];     // year ... since? gregorian?
	stop_tm.tm_mon  = file_header->ndaysk[1] - 1; // calendar month (offset by 1 to make january = 0)
	stop_tm.tm_mday = file_header->ndaysk[2];     // day of month
	stop_tm.tm_hour = file_header->ntimsk[0];     // hour of day (24hr)
	stop_tm.tm_min  = file_header->ntimsk[1];     // minutes
	stop_tm.tm_sec  = file_header->ntimsk[2];     // seconds
	// build timestamp
	last_event_time = mktime(&stop_tm);    // 
	
	// return run duration in seconds
	return (stop_time - start_time);
}

// a couple of last helper functions before the main business end of this tool
std::streambuf* start_stderr_capture(){
	// redirect stderr to a custom stringbuf for capturing return from subsequent function calls
	// =========================================================================================
	// declare a new stringstream to catch cerr and point stderr at it
	// capture the return, which notes where it was previously
	return std::cerr.rdbuf(new std::stringbuf);
}

std::string end_stderr_capture(std::streambuf* previous_buff){
	// stop capturing stderr, retrieve output, and return normal behaviour
	// ===================================================================
	// restore stderr to its previous state, and retreieve our own streambuf
	std::stringbuf* capture_buffer = static_cast<std::stringbuf*>(std::cerr.rdbuf(previous_buff));
	// for backwards compatibility, print whatever was returned via the normal channel
	std::cerr<<capture_buffer->str();
	// get the captured contents
	std::string capture = capture_buffer->str();
	// free the temporary buffer
	delete capture_buffer;
	return capture;
}

int ShowNeutrons::CreateOutputFile(std::string filename){
	// create the output ROOT file and TTree for writing
	// =================================================
	outfile = new TFile(filename.c_str(), "RECREATE");
	outtree = new TTree("eventtree", "Events with Neutron Captures");
	
	// Set up all the output TTree branches
	// Each TTree entry will correspond to one event, which may have multiple neutrons
	// Each neutron capture may have multiple gammas
	// Each gamma may have multiple PMT hits
	
/*
	// create pointers to use for branch addresses 
	// XXX ONLY USE POINTERS FOR READING TTREES XXX
	// --------------------------------------------
	std::string* filenamep = &filename;
	std::map<int,std::string>* neutron_process_mapp = &neutron_process_map;
	
	// primary particles
//	std::vector<int>* primary_idp = &primary_id;
	std::vector<int>* primary_pdgp = &primary_pdg;
	std::vector<double>* primary_energyp = &primary_energy; // [MeV]
	std::vector<TLorentzVector>* primary_start_posp = &primary_start_pos; // [cm]
	std::vector<TLorentzVector>* primary_end_posp = &primary_end_pos; // [cm]
	
	// parent nuclide
//	std::vector<int>* parent_primary_idp = &parent_primary_id;
//	std::vector<int>* nuclide_idp = &nuclide_id;
	std::vector<int>* nuclide_pdgp = &nuclide_pdg;
	std::vector<TLorentzVector>* nuclide_creation_posp = &nuclide_creation_pos;
	std::vector<TLorentzVector>* nuclide_decay_posp = &nuclide_decay_pos;
	
	// neutrons
//	std::vector<std::vector<int> >* parent_nuclide_idp = &parent_nuclide_id;
//	std::vector<std::vector<int> >* neutron_idp = &neutron_id;
	std::vector<std::vector<TLorentzVector> >* neutron_start_posp = &neutron_start_pos;
	std::vector<std::vector<TLorentzVector> >* neutron_end_posp = &neutron_end_pos;
	std::vector<std::vector<double> >* neutron_start_energyp = &neutron_start_energy;
	std::vector<std::vector<double> >* neutron_end_energyp = &neutron_end_energy;
	std::vector<std::vector<int> >* neutron_end_processp = &neutron_end_process;
	
	// gammas
//	std::vector<std::vector<std::vector<int> > >* parent_neutron_idp = &parent_neutron_id;
	std::vector<std::vector<std::vector<double> > >* gamma_energyp = &gamma_energy;
*/
	
	// create branches
	// ---------------
	// file level
	outtree->Branch("filename",&filename);
	outtree->Branch("skdetsim_version",&skdetsim_version);
	outtree->Branch("tba_table_version",&tba_table_version);
	
	// reference constant - this is built as we go so always use the LAST TTree entry
	// TODO put it in the TTree metadata or something
	outtree->Branch("neutron_process_map",&neutron_process_map);
	
	// event level
	outtree->Branch("entry_number",&entry_number);
	
	// primary particle
//	outtree->Branch("primary_id",&primary_id);
	outtree->Branch("primary_pdg",&primary_pdg);
	outtree->Branch("primary_energy",&primary_energy);
	outtree->Branch("primary_start_pos",&primary_start_pos);
	outtree->Branch("primary_end_pos",&primary_end_pos);
	
	// parent nuclide - one for each neutron
//	outtree->Branch("parent_primary_id",&parent_primary_id);
//	outtree->Branch("nuclide_id",&nuclide_id);
	outtree->Branch("nuclide_pdg",&nuclide_pdg);
	outtree->Branch("nuclide_creation_pos",&nuclide_creation_pos);
	outtree->Branch("nuclide_decay_pos",&nuclide_decay_pos);
	
	// neutron
//	outtree->Branch("parent_nuclide_id",&parent_nuclide_id);
//	outtree->Branch("neutron_id",&neutron_id);
	outtree->Branch("neutron_start_pos",&neutron_start_pos);
	outtree->Branch("neutron_end_pos",&neutron_end_pos);
	outtree->Branch("neutron_start_energy",&neutron_start_energy);
	outtree->Branch("neutron_end_process",&neutron_end_process);
	
	// gamma
//	outtree->Branch("parent_neutron_id",&parent_neutron_id);
	outtree->Branch("gamma_energy",&gamma_energy);
	
	return 1;
}

int ShowNeutrons::ProcessEntry(int entry_i){
	current_entry_num = entry_i;
	if(verbosity>2) std::cout<<"Processing entry "<<current_entry_num<<std::endl;
	// Load event
	// ---------------
	// for some reason the TreeManager requires two calls to load the entry.
	// the fortran interface mimicks (wraps) the same calls, so unless we work with
	// the TTree directly, this is necessarily a two-step process.
	// The first step (skroot_next_entry or skroot_jump_entry) sets the event number
	// and clears the branch objects by invoking TreeManager::Clear()
	// XXX if not using the TreeManager, be sure to call Clear() on all objects!!! XXX
	//skroot_next_entry_(&TreeManager_ID,&get_ok);
	skroot_jump_entry_(&TreeManager_ID,&current_entry_num,&get_ok);
	// It gets set to TRUE when we run off the end of the TTree
	if(get_ok){
		std::cout<<"End of processing entries, returning"<<std::endl;
		return 1;
	}
	// 2. now load the objects with the information from this event.
	if(verbosity>2) std::cout<<"skroot get next entry"<<std::endl;
	skroot_get_entry_(&TreeManager_ID);  // TODO wrap steps 1. & 2. in a single call.
	
	// clear all the output vectors
	if(verbosity>2) std::cout<<"clearing output vectors"<<std::endl;
	ClearOutputTreeBranches();
	
	// file information
	if(verbosity>2) std::cout<<"getting input filename"<<std::endl;
	TTree* tree = skroot_get_tree(&TreeManager_ID);
	filename = std::string(tree->GetCurrentFile()->GetName());
	
	// event information
	if(verbosity>2) std::cout<<"retrieving event information"<<std::endl;
//	event_number = file_header->nevsk;
//	swtrg_id = file_header->swtrg_id;
//	run_number = 
//	subrun_number = 
//	// build a timestamp and calculate time since last event?
//	// using PrevT0?
	
	entry_number = current_entry_num;
	
	// MC truth information is accessed from the MCInfo struct
	// by directly accessing its members - no getters here
	
	// alright, this script is a bust because the MCInfo only stores
	// NEUT primaries - no secondaries, so no neutrons, no gammas, etc.
	// 
	// If we had the information, first we'd scan through and record neutrons.
	// If no neutrons, continue to next event.
	// Otherwise, re-scan and record parent and decay gamma information.
	
	// As it is, we can't do this here, so i guess just make human readable
	// what information is present in the MCInfo struct, and maybe others later.
	
	// ========================================
	
	// note: we only have primary particles here
	int num_primaries = mc_info->nvc;
	// for some reason 'vertices' are independent of 'particles'
	// in practice it just means using some index juggling when retrieving information
	int num_vertices = mc_info->nvtxvc;
	
	// a typical IBD event will have the electron anti-neutrino (pdg -12)
	// the positron (pdg -11) and the neutron (pdg 2112)
	if(verbosity>2) std::cout<<"looping over "<<num_primaries<<" primary particles"<<std::endl;
	for(int primary_i=0; primary_i<num_primaries; ++primary_i){
		
		primary_pdg.push_back(mc_info->ipvc[primary_i]);
		primary_energy.push_back(mc_info->energy[primary_i]); // fixme
		// energy might not be populated, calculate it from momentum:
		TVector3 momentum_vector(mc_info->pvc[primary_i][0],
					 mc_info->pvc[primary_i][1],
					 mc_info->pvc[primary_i][2]);
		// we'll need to look up the mass
		double mass=0;
		if(pdg_to_mass.count(primary_pdg.back())){
			mass = pdg_to_mass.at(primary_pdg.back());
		}
		double energy_calc = sqrt(momentum_vector.Mag2()+pow(mass,2));
		
		double pos_x=-999, pos_y=-999, pos_z=-999, pos_t=-999;
		// positions and times need to be looked up from the vertex array
		// - only the vertex index is stored in the particle array
		// i *think* the indices are fortran-style indexes (start from 1)
		// note also the number of vertices is specified separately from
		// the number of particles, so check we're within bounds
		int start_vertex_index = mc_info->ivtivc[primary_i]-1;
		if(start_vertex_index>=0 && start_vertex_index<num_vertices){
			pos_x = mc_info->pvtxvc[start_vertex_index][0];
			pos_y = mc_info->pvtxvc[start_vertex_index][1];
			pos_z = mc_info->pvtxvc[start_vertex_index][2];
			pos_t = mc_info->timvvc[start_vertex_index];
		}
		primary_start_pos.push_back(TLorentzVector(pos_x,pos_y,pos_z,pos_t));
		
		// repeat for the end vertex
		pos_x=-999, pos_y=-999, pos_z=-999, pos_t=-999;
		int end_vertex_index = mc_info->ivtfvc[primary_i]-1;
		if(end_vertex_index>=0 && end_vertex_index<num_vertices){
			pos_x = mc_info->pvtxvc[end_vertex_index][0];
			pos_y = mc_info->pvtxvc[end_vertex_index][1];
			pos_z = mc_info->pvtxvc[end_vertex_index][2];
			pos_t = mc_info->timvvc[end_vertex_index];
		}
		primary_end_pos.push_back(TLorentzVector(pos_x,pos_y,pos_z,pos_t));
		
//		// other available info: for vertices (todo: bounds check on index)
//		// 1. "kind of vertex" ???
//		int start_vertex_kind = mc_info->iflvvc[start_vertex_index];
//		// 2. "parent particle" ??? << index? of what? pdg code?
//		int parent_particle_something = mc_info->iparvc[start_vertex_index];
//		
//		// other available info: for particles
//		// 1. ID of origin particle parent particle ????
//		int id_of_something = mc_info->iorgvc[primary_i];
//		// 2. Cherenkov flag ...? int or bool? what is it flagging? "is it producing cherenkov"?
//		int cherenkov_flag = mc_info->icrnvc[primary_i];
//		// 3. final state flag ...? int or bool? what about the final state is this flagging?
//		int final_state_flag = mc_info->iflgvc[primary_i];
	}
	
//	// departure process
//	n_end_process_name = "some_meaningful_string";
//	
//	// on encountering a new process, add it to the map
//	int process_keynum=-1;
//	for(auto&& aprocess : neutron_process_map){
//		if(aprocess.second==n_end_process_name){
//			process_keynum=aprocess.first;
//			break;
//		}
//	}
//	if(process_keynum<0){
//		// add this process to the map
//		process_keynum=neutron_process_map.size();
//		neutron_process_map.emplace(process_keynum,n_end_process_name);
//	}
	
	if(verbosity>2) PrintBranches();
	
	// Fill the tree
	if(verbosity>2) std::cout<<"filling output TTree entry"<<std::endl;
	outtree->Fill();
	
	// update the output file
	if((current_entry_num%WRITE_FREQUENCY)==0) WriteTree();
	
	// indicate end of tree
	return (current_entry_num==(total_entries-1));
}

int ShowNeutrons::GenerateHistograms(){
	// here be where we make the interesting plots
	// XXX TODO XXX
	return 1;
}

void ShowNeutrons::ClearOutputTreeBranches(){
//	primary_id.clear();
	primary_pdg.clear();
	primary_energy.clear();
	primary_start_pos.clear();
	primary_end_pos.clear();
	
	// parent nuclide
//	parent_primary_id.clear();
//	nuclide_id.clear();
	nuclide_pdg.clear();
	nuclide_creation_pos.clear();
	nuclide_decay_pos.clear();
	
	// neutrons
//	parent_nuclide_id.clear();
//	neutron_id.clear();
	neutron_start_pos.clear();
	neutron_end_pos.clear();
	neutron_start_energy.clear();
	neutron_end_energy.clear();
	neutron_end_process.clear();
	
	// gammas
//	parent_neutron_id.clear();
	gamma_energy.clear();
	
	return;
}

ShowNeutrons::~ShowNeutrons(){
	// Reset output TTree branch addresses. Probably redundant, but no harm in doing so.
	CloseFile();
}

void ShowNeutrons::PrintBranches(){
	std::cout<<"==========================================================="<<std::endl;
	std::cout<<"PRINTING EVENT"<<std::endl;
	std::cout<<"filename: "<<filename<<std::endl
		 <<"skdetsim_version: "<<skdetsim_version<<std::endl
		 <<"tba_table_version: "<<tba_table_version<<std::endl
		 <<"neutron_process_map: {";
	for(std::map<int,std::string>::const_iterator aprocess=neutron_process_map.begin();
	    aprocess!=neutron_process_map.end(); ++aprocess){
		std::cout<<"["<<aprocess->first<<"="<<aprocess->second<<"], ";
	}
	std::cout<<"\b\b}"<<std::endl;
	std::cout<<"entry_number: "<<entry_number<<std::endl
		 <<"num primaries: "<<primary_pdg.size()<<std::endl;
		 
	// print primaries
	for(int primary_i=0; primary_i<primary_pdg.size(); ++primary_i){
		std::cout<<"\tprimary ("<<primary_i<<"): "<<std::endl
			 <<"\t\tprimary pdg: "<<primary_pdg.at(primary_i)<<std::endl
			 <<"\t\tprimary energy: "<<primary_energy.at(primary_i)<<std::endl
			 <<"\t\tprimary start pos:"
			 <<" ("<<primary_start_pos.at(primary_i).X()
			 <<", "<<primary_start_pos.at(primary_i).Y()
			 <<", "<<primary_start_pos.at(primary_i).Z()<<")"<<std::endl
			 <<"\t\tprimary end pos:"
			 <<" ("<<primary_end_pos.at(primary_i).X()
			 <<", "<<primary_end_pos.at(primary_i).Y()
			 <<", "<<primary_end_pos.at(primary_i).Z()<<")"<<std::endl;
	}
	
	int total_neutrons=0;
	int total_gammas=0;
	// print nuclides
	std::cout<<"num nuclides: "<<nuclide_pdg.size()<<std::endl;
	for(int nuclide_i=0; nuclide_i<nuclide_pdg.size(); ++nuclide_i){
		std::cout<<"\tnuclide ("<<nuclide_i<<"): "<<std::endl
			 <<"\t\tnuclide pdg: "<<nuclide_pdg.at(nuclide_i)<<std::endl
			 <<"\t\tnuclide energy: "<<primary_energy.at(nuclide_i)<<std::endl
			 <<"\t\tnuclide start pos:"
			 <<" ("<<nuclide_creation_pos.at(nuclide_i).X()
			 <<", "<<nuclide_creation_pos.at(nuclide_i).Y()
			 <<", "<<nuclide_creation_pos.at(nuclide_i).Z()<<")"<<std::endl
			 <<"\t\tnuclide end pos:"
			 <<" ("<<nuclide_decay_pos.at(nuclide_i).X()
			 <<", "<<nuclide_decay_pos.at(nuclide_i).Y()
			 <<", "<<nuclide_decay_pos.at(nuclide_i).Z()<<")"<<std::endl;
		
		// each nuclide may have multiple neutrons
		// (every neutron should have a nuclide ... backgrounds?)
		std::cout<<"\tnum neutrons from this nuclide: "<<neutron_start_energy.at(nuclide_i).size()<<std::endl;
		for(int neutron_i=0; neutron_i<neutron_start_energy.at(nuclide_i).size(); ++neutron_i){
			std::cout<<"\t\t("<<neutron_i<<"): "<<std::endl
				 <<"\t\t\tneutron end process: "<<neutron_end_process.at(nuclide_i).at(neutron_i)<<std::endl
				 <<"\t\t\tneutron start energy: "<<neutron_start_energy.at(nuclide_i).at(neutron_i)<<std::endl
				 <<"\t\t\tneutron start pos:"
				 <<" ("<<neutron_start_pos.at(nuclide_i).at(neutron_i).X()
				 <<", "<<neutron_start_pos.at(nuclide_i).at(neutron_i).Y()
				 <<", "<<neutron_start_pos.at(nuclide_i).at(neutron_i).Z()<<")"<<std::endl
				 <<"\t\t\tneutron end pos:"
				 <<" ("<<neutron_end_pos.at(nuclide_i).at(neutron_i).X()
				 <<", "<<neutron_end_pos.at(nuclide_i).at(neutron_i).Y()
				 <<", "<<neutron_end_pos.at(nuclide_i).at(neutron_i).Z()<<")"<<std::endl;
			total_neutrons++;
			
			// each neutron may have multiple decay gammas
			std::cout<<"\t\tnum gammas from this neutron: "
				 <<gamma_energy.at(nuclide_i).at(neutron_i).size()<<std::endl;
			for(int gamma_i=0; gamma_i<gamma_energy.at(nuclide_i).at(neutron_i).size(); ++gamma_i){
				std::cout<<"\t\t\t("<<gamma_i<<"): "<<std::endl
				 <<"\t\t\t\tgamma energy: "<<gamma_energy.at(nuclide_i).at(neutron_i).at(gamma_i)<<std::endl;
				total_gammas++;
			}
		}
	}
	std::cout<<"total neutrons in event: "<<total_neutrons<<std::endl
		 <<"total gammas in event: "<<total_gammas<<std::endl;
	// shallow sanity checks
	std::cout<<"consistency check: num neutrons outer vs num nuclides: "
		 <<(neutron_start_energy.size()==nuclide_pdg.size())<<std::endl
		 <<"consistency check: num gammas outer vs num nuclides: "
		 <<(gamma_energy.size()==nuclide_pdg.size())<<std::endl;
	std::cout<<"==========================================================="<<std::endl;
}
