/* vim:set noexpandtab tabstop=4 wrap */
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

int main(int argc, const char* argv[]){
	/* demo script to check we can read the resulting TTree */
	
	// check inputs
	if(argc==1){
		std::cout<<"usage: "<<argv[0]<<" <input_file>"<<std::endl;
		return 0;
	}
	
	// ###################################################################
	// OPEN FILE, GET TTREE
	// ###################################################################
	
	// get input file and retrieve TTree
	TFile* f = TFile::Open(argv[1]);
	TTree* t = (TTree*)f->Get("eventtree");
	
	// ###################################################################
	// DECLARE TTREE VARIABLES
	// ###################################################################
	
	// this stuff is copied from show_neutrons.hh
	// misc - retreieve from the last entry, supposedly
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
	
	// ###################################################################
	// GET VARIABLE POINTERS
	// ###################################################################
	
	// create pointers to use for branch addresses 
	// XXX ONLY USE ADDRESSES TO POINTERS WITH TTree::SetBranchAddress XXX
	// XXX        USE ADDRESSES TO OBJECTS WITH TTree::Branch
	// -------------------------------------------------------------------
	std::map<int,std::string>* neutron_process_mapp = &neutron_process_map;
	std::string* filenamep = &filename;
	
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
	
	// ###################################################################
	// SETUP TREE FOR READING
	// ###################################################################
	
	t->SetBranchAddress("filename",&filenamep);
	t->SetBranchAddress("skdetsim_version",&skdetsim_version);
	t->SetBranchAddress("tba_table_version",&tba_table_version);
	t->SetBranchAddress("neutron_process_map",&neutron_process_mapp);
	t->SetBranchAddress("entry_number",&entry_number);
//	t->SetBranchAddress("primary_id",&primary_idp);
	t->SetBranchAddress("primary_pdg",&primary_pdgp);
	t->SetBranchAddress("primary_energy",&primary_energyp);
	t->SetBranchAddress("primary_start_pos",&primary_start_posp);
	t->SetBranchAddress("primary_end_pos",&primary_end_posp);
//	t->SetBranchAddress("parent_primary_id",&parent_primary_idp);
//	t->SetBranchAddress("nuclide_id",&nuclide_idp);
	t->SetBranchAddress("nuclide_pdg",&nuclide_pdgp);
	t->SetBranchAddress("nuclide_creation_pos",&nuclide_creation_posp);
	t->SetBranchAddress("nuclide_decay_pos",&nuclide_decay_posp);
//	t->SetBranchAddress("parent_nuclide_id",&parent_nuclide_idp);
//	t->SetBranchAddress("neutron_id",&neutron_idp);
	t->SetBranchAddress("neutron_start_pos",&neutron_start_posp);
	t->SetBranchAddress("neutron_end_pos",&neutron_end_posp);
	t->SetBranchAddress("neutron_start_energy",&neutron_start_energyp);
	t->SetBranchAddress("neutron_end_process",&neutron_end_processp);
//	t->SetBranchAddress("parent_neutron_id",&parent_neutron_idp);
	t->SetBranchAddress("gamma_energy",&gamma_energyp);
	
	// ###################################################################
	// PRINT EVENTS
	// ###################################################################
	
	for(int entry_i=0; entry_i<t->GetEntries(); ++entry_i){
		t->GetEntry(entry_i);
		
		// literally a rip of ShowNeutrons::PrintBranches()
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
	
	// ###################################################################
	// CLEANUP
	// ###################################################################
	
	std::cout<<"end of tree"<<std::endl;
	t->ResetBranchAddresses();
	f->Close();
	
	// ###################################################################
	// C++11 test
	// ###################################################################
	// just checking c++11 stuff is working
	std::map<int,std::string> somevals;
	somevals.emplace(5,"potatoes");
	somevals.emplace(6,"courgettes");
	for(auto&& apair : somevals){
		std::cout<<"next pair: "<<apair.first<<"="<<apair.second<<std::endl;
	}
	std::cout<<"we have "<<somevals.count(5)<<" entries with a key of 5"<<std::endl;
	return 1;
}
