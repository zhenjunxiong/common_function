#include "histo_read.h"

using namespace std;

void test(){
TFile* file = new TFile("/home/zjx/workdir/run3/test_code/bk.root","read");
string hashlistname = "table-maker/output";
string listname = "TrackBarrel_jpsiO2MCdebugCuts5";
TList* list = get_List(file,hashlistname,listname);

TH1F* Ele_pin = (TH1F*)list->FindObject("TPCnSigEle_pIN");
TF1* Ele_pin_fit_low = new TF1("Ele_pin_fit_low","1/([0]+[1]*x)+[2]+[3]*x-3",0.3,15);
TF1* Ele_pin_fit_high = new TF1("Ele_pin_fit_high","1/([0]+[1]*x)+[2]+[3]*x+3",0.3,15);
Ele_pin_fit_low->SetParameters(0.65,-5.37,-0.86,0.16);
Ele_pin_fit_high->SetParameters(0.65,-5.37,-0.86,0.16);

Ele_pin->Draw();
Ele_pin_fit_low->Draw("same");
Ele_pin_fit_high->Draw("same");
}