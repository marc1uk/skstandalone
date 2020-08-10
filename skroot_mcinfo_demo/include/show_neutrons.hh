/* vim:set noexpandtab tabstop=4 wrap */
#ifndef SHOW_NEUTRONS_H
#define SHOW_NEUTRONS_H

#include "skroot.h"

#include <iostream> // std::cout, cerr
#include <string>
#include <vector>
#include <utility>  // std::pair
#include <map>
#include <time.h>   // for time_t, for determining run periods

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

class ShowNeutrons{
	public:
	ShowNeutrons(int tree_man_id) : TreeManager_ID(tree_man_id){};
	~ShowNeutrons();
	
	// FUNCTIONS
	// =========
	int LoadSkrootFile(std::vector<std::string> input_files);
	double GetTdiffFirstLastEvent(time_t &first_event_time, time_t &last_event_time);
	int CreateOutputFile(std::string filename);
	int ProcessEntry(int entry_i);      // the main event loop
	int GenerateHistograms();  // make and write out the analysis plots
	void ClearOutputTreeBranches();
	void PrintBranches();
	
	// Getters and Setters
	int GetMCVersion();
	int GetEntryNum();
	int GetTotalEntries();
	int WriteTree();
	void CloseFile();
	void SetVerbosity(int verb);
	
	// VARIABLES
	// =========
	private:
	
	// a unique file handle identifier for this TreeManager
	int TreeManager_ID;
	TreeManager* tree_manager=nullptr;
	
	// input file parameters
	int total_entries;   // in all loaded files
	int current_entry_num = 0;
	time_t  start_time, stop_time;
	double run_duration_seconds;
	Header   *file_header = nullptr;
	MCInfo   *mc_info     = nullptr;
//	LoweInfo *lowe_info   = nullptr;
//	MuInfo   *muon_info   = nullptr;
//	RunInfo* run_info    = nullptr;
//	PrevT0   *prev_t0     = nullptr;
//	TClonesArray* hw_trig_list = nullptr;
//	SoftTrgList*sw_trig_list = nullptr;
//	SoftwareTrigger* sw_trig      = nullptr;
//	EventHeader* event_header = nullptr;
//	EventTrailer* event_trailer = nullptr;
//	TQReal* tq_real_info = nullptr;
//	TQReal* tq_areal_info = nullptr;
	
	// output file parameters
	TFile *outfile = nullptr;
	TTree *outtree = nullptr;
	// see bottom for output variables, since we'll have a lot of them.
	
	
	// verbosity and error handling
	int verbosity=10;
	int get_ok;
	
	// OUTPUTS
	// =======
	// variables we wish to be in the output tree
	
	// one-off map of neutron end processes from int to string
	std::map<int,std::string> neutron_process_map;
	
	// file-wise
	std::string filename;
	int skdetsim_version;
	int tba_table_version;
	
	// event-wise
	int entry_number;   // TTree entry number to be able to identify the source event
	
	// primary particle - an event can have many primary particles
//	std::vector<int> primary_id;
	std::vector<int> primary_pdg;
	std::vector<double> primary_energy; // [MeV]
	std::vector<TLorentzVector> primary_start_pos; // [cm, ns]
	std::vector<TLorentzVector> primary_end_pos; // [cm, ns]
	
	// parent nuclide - potentially many per event
//	std::vector<int> parent_primary_id;
//	std::vector<int> nuclide_id; // so we can pair neutrons to their parent
	std::vector<int> nuclide_pdg;
	std::vector<TLorentzVector> nuclide_creation_pos; // [cm, ns]  - likely to not go far but just in case
	std::vector<TLorentzVector> nuclide_decay_pos; // [cm, ns]
	
	// neutrons  - potentially many per nuclide? or just one?
//	std::vector<std::vector<int> > parent_nuclide_id;
//	std::vector<std::vector<int> > neutron_id; // so we can pair gammas to their neutrons
	std::vector<std::vector<TLorentzVector> > neutron_start_pos; // [cm, ns]
	std::vector<std::vector<TLorentzVector> > neutron_end_pos; // [cm, ns]
	std::vector<std::vector<double> > neutron_start_energy; // [MeV]
	std::vector<std::vector<double> > neutron_end_energy; // [MeV] - on capture/decay
	std::vector<std::vector<int> > neutron_end_process; // enum to map to e.g. capture, decay, escape world...
	
	// gammas - potentially many per neutron
//	std::vector<std::vector<std::vector<int> > > parent_neutron_id;
	std::vector<std::vector<std::vector<double> > > gamma_energy; // [MeV]
	
//	// probably don't actually care about these
//	int run_number;
//	int subrun_number;
//	int swtrg_id;
//	int event_number;
	
//	// what are these?
//	mdrnsk = file_header->mdrnsk;
//	idtgsk = file_header->idtgsk;
//	ifevsk = file_header->ifevsk;
	
	// detector information
	// total charge? time distribution of hits?
	
};
#endif

/* OK plots we want to be able to make:

Maybe
-----
+ neutron multiplicity?
	- vs neutron creation energy
	- categorise by source nuclide?
	- vs primary particle pdg (probably a muon)
	- vs primary particle energy

Later
-----
+ time to parent nuclide decay
	- parent nuclide time of creation
	- parent nuclide time of decay
- distance to parent nuclide decay
	- parent nuclide creation position
	- parent nuclide decay position

Now
---
+ neutron energy at creation
	- categorise by source: parent nuclei? spallation from primary muon? decay from a spallation-induced isotope?
		- primary particle pdg?
		- primary particle energy?
+ neutron energy at capture
	- categorise by capture nucleus
+ distance to neutron capture
	- neutron position at creation
	- neutron position at capture
	- vs neutron energy
	- split by x, y, z, where z is neutron initial direction, or primary particle initial direction
+ time to neutron capture
	- neutron time of creation
	- neutron time of capture
	- vs neutron energy
+ gamma energy
	- categorise by capture nucleus (H, Gd 157, 158...)
	- vs time since capture
+ gamma multiplicity
	- categorise by capture nucleus
	- vs gamma energy

*/
