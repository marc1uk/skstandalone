#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include "TVector3.h"

#include "skroot.h"

double calcAbsTime(unsigned int t0, unsigned int cnt32){
	double abst0 = (double)t0*0.5208333*1e-9;
	double abst1 = ( double )( cnt32 >> 17 ) 
               * 0.5208333 * 32768. * 131072. * 1.e-9;
	return (abst0 + abst1);
}

int main(int argc, char* argv[])
{
  // check arguments
  if(argc < 3){
    cerr << "Usage: " << basename(argv[0]) << "  <output_root_name> <input skroot files...>" << endl;
    exit(1);
  }

  string outname(argv[1]);

  // open output file
  TFile *outfile = new TFile(outname.c_str(), "recreate");
  TTree *outtree = new TTree("outtree", "tree");

  int nrunsk;
  int nsubsk;
  int nevsk;
  int ndaysk0;
  int ndaysk1;
  int ndaysk2;
  int ntimsk0;
  int ntimsk1;
  int ntimsk2;
  int mdrnsk;
  int idtgsk;
  int ifevsk;
  int swtrg_id;

  double bsgood;
  double bsdirks;
  Double_t bsenergy;
  double bsvertex0;
  double bsvertex1;
  double bsvertex2;
  double bsdir0;
  double bsdir1;
  double bsdir2;
  double bsresult;
  int bsn50;
  double bscossun;
  double bswallsk;
  double ltimediff;
  double lnsratio;
  int spaevnum;
  double spaloglike;
  double spadt;
  double spadll;
  double spadlt;
  double spamuyn;
  double spamugdn;

  double muqismsk;
  double muinfo[200];
  double muboy_dedx[200];

  outtree->Branch("nrunsk",&nrunsk,"nrunsk/I");
  outtree->Branch("nsubsk",&nsubsk,"nsubsk/I");
  outtree->Branch("nevsk",&nevsk,"nevsk/I");
  outtree->Branch("ndaysk0",&ndaysk0,"ndaysk0/I");
  outtree->Branch("ndaysk1",&ndaysk1,"ndaysk1/I");
  outtree->Branch("ndaysk2",&ndaysk2,"ndaysk2/I");
  outtree->Branch("ntimsk0",&ntimsk0,"ntimsk0/I");
  outtree->Branch("ntimsk1",&ntimsk1,"ntimsk1/I");
  outtree->Branch("ntimsk2",&ntimsk2,"ntimsk2/I");
  outtree->Branch("mdrnsk",&mdrnsk,"mdrnsk/I");
  outtree->Branch("idtgsk",&idtgsk,"idtgsk/I");
  outtree->Branch("ifevsk",&ifevsk,"ifevsk/I");
  outtree->Branch("swtrg_id",&swtrg_id,"swtrg_id/I");

  outtree->Branch("bsenergy", &bsenergy, "bsenergy/D");
  outtree->Branch("bsvertex0",&bsvertex0,"bsvertex0/D");
  outtree->Branch("bsvertex1",&bsvertex1,"bsvertex1/D");
  outtree->Branch("bsvertex2",&bsvertex2,"bsvertex2/D");
  outtree->Branch("bsdir0",&bsdir0,"bsdir0/D");
  outtree->Branch("bsdir1",&bsdir1,"bsdir1/D");
  outtree->Branch("bsdir2",&bsdir2,"bsdir2/D");
  outtree->Branch("bsresult",&bsresult,"bsresult/D");
  outtree->Branch("bsgood",&bsgood,"bsgood/D");
  outtree->Branch("bsdirks",&bsdirks,"bsdirks/D");
  outtree->Branch("bsn50",&bsn50,"bsn50/I");
  outtree->Branch("bscossun",&bscossun,"bscossun/D");
  outtree->Branch("bswallsk",&bswallsk,"bswallsk/D");
  outtree->Branch("litimediff",&ltimediff,"ltimediff/D");
  outtree->Branch("lnsratio",&lnsratio,"lnsratio/D");
  outtree->Branch("spaevnum",&spaevnum,"spaevnum/I");
  outtree->Branch("spaloglike",&spaloglike,"spaloglike/D");
  outtree->Branch("spadt",&spadt,"spadt/D");
  outtree->Branch("spadll",&spadll,"spadll/D");
  outtree->Branch("spadlt",&spadlt,"spadlt/D");
  outtree->Branch("spamuyn",&spamuyn,"spamuyn/D");
  outtree->Branch("spamugdn",&spamugdn,"spamugdn/D");

  outtree->Branch("muqismsk",&muqismsk,"muqismsk/D");

  // open skroot file (read only)
  int id = 10;
  skroot_open_read(&id);
  for(int i=2; i<argc; i++){
    string fname_in = argv[i];
    skroot_set_input_file(&id, fname_in.c_str(), fname_in.size());
  }
  skroot_initialize(&id);

  // obtain tree
  TreeManager* mgr = skroot_get_mgr(&id);
  TTree* tree = mgr->GetTree();

  // obtain off-line branches
  Header   *HEAD = mgr->GetHEAD();
  LoweInfo *LOWE = mgr->GetLOWE();
  MuInfo   *MU   = mgr->GetMU();

  // total number of events
  int nentries = tree->GetEntries();

  std::cout << "nentries = " << nentries << std::endl;
  // main loop
  for(int iev=0; iev<nentries; iev++){
    // read all branches
    tree->GetEntry(iev);

    //int eventtime_s = LOWE->linfo[61];
    //int eventtime_ns = LOWE->linfo[62];

    //double time = calcAbsTime(HEAD->t0,HEAD->counter_32);

    nrunsk = HEAD->nrunsk;
    nsubsk = HEAD->nsubsk;
    nevsk = HEAD->nevsk;
    ndaysk0 = HEAD->ndaysk[0];
    ndaysk1 = HEAD->ndaysk[1];
    ndaysk2 = HEAD->ndaysk[2];
    ntimsk0 = HEAD->ntimsk[0];
    ntimsk1 = HEAD->ntimsk[1];
    ntimsk2 = HEAD->ntimsk[2];
    mdrnsk = HEAD->mdrnsk;
    idtgsk = HEAD->idtgsk;
    ifevsk = HEAD->ifevsk;
    swtrg_id = HEAD->swtrg_id;

    bsenergy = LOWE->bsenergy;
    bsgood = LOWE->bsgood[1];
    bsdirks = LOWE->bsdirks;
    bsvertex0 = LOWE->bsvertex[0];
    bsvertex1 = LOWE->bsvertex[1];
    bsvertex2 = LOWE->bsvertex[2];
    bsdir0 = LOWE->bsdir[0];
    bsdir1 = LOWE->bsdir[1];
    bsdir2 = LOWE->bsdir[2];
    bsresult = LOWE->bsresult[0];
    bsn50 = LOWE->bsn50;
    bscossun = LOWE->bscossun;
    bswallsk = *reinterpret_cast<float*>(&LOWE->linfo[9]);
    ltimediff = LOWE->ltimediff;
    lnsratio = LOWE->lnsratio;
    spaevnum = LOWE->spaevnum;
    spaloglike = LOWE->spaloglike;
    spadt = LOWE->spadt;
    spadll = LOWE->spadll;
    spadlt = LOWE->spadlt;
    spamuyn = LOWE->spamuyn;
    spamugdn = LOWE->spamugdn;

    int ith = 10; 
    cout << *reinterpret_cast<float*>(&MU->muinfo[ith]) 
         << " " 
	 << MU->muboy_dedx[ith]
	 << endl; 

    muqismsk = MU->muqismsk;

    outtree->Fill();
  }

  outfile->cd();

  // close output file
  outtree->Write();
  outfile->Close();

  // close skroot file
  skroot_close(&id);
  skroot_end();

  return 0;
}
