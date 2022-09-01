#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
using namespace std;
TTree* getbranch4polarization(TTree* input_tree);

//void getbranch4palorization();
float pairpt,paireta,pairphi,trk1pt,trk2pt,trk1eta,trk2eta,trk1phi,trk2phi;
namespace polarization
{
    double get_PhiHelicityFrame(TLorentzVector electronPositive, TLorentzVector electronNegative, TLorentzVector possibleJPsi );
    double get_CosThetaHelicityFrame(TLorentzVector electronPositive, TLorentzVector electronNegative, TLorentzVector possibleJPsi);
    double get_PhiCollinsSoper(TLorentzVector electronPositive, TLorentzVector electronNegative, TLorentzVector possibleJPsi);
    double get_CosThetaCollinsSoper(TLorentzVector electronPositive, TLorentzVector electronNegative, TLorentzVector possibleJPsi);
};

TTree* getbranch4polarization(TTree* input_tree)
{   
    // TFile *mc17 = new TFile("/home/zjx/workdir/polarization/input_test.root","read");
	// TTree *input_tree =  (TTree*) mc17->Get("fTreePairsSPM");
    double PhiHelicityFrame,CosThetaHelicityFrame,PhiCollinsSoper,CosThetaCollinsSoper;
    // TFile *output_root = new TFile("/home/zjx/workdir/polarization/output_test.root","recreate");
    TTree* output_tree = input_tree->CloneTree(0);
    output_tree->Branch("PhiHelicityFrame",&PhiHelicityFrame,"PhiHelicityFrame/D");
    output_tree->Branch("CosThetaHelicityFrame",&CosThetaHelicityFrame,"CosThetaHelicityFrame/D");
    output_tree->Branch("PhiCollinsSoper",&PhiCollinsSoper,"PhiCollinsSoper/D");
    output_tree->Branch("CosThetaCollinsSoper",&CosThetaCollinsSoper,"CosThetaCollinsSoper/D");

    input_tree->SetBranchAddress("pairPt",&pairpt);
    input_tree->SetBranchAddress("pairEta",&paireta);
    input_tree->SetBranchAddress("pairPhi",&pairphi);
    input_tree->SetBranchAddress("trk1Pt",&trk1pt);
    input_tree->SetBranchAddress("trk1Eta",&trk1eta);
    input_tree->SetBranchAddress("trk1Phi",&trk1phi);
    input_tree->SetBranchAddress("trk2Pt",&trk2pt);
    input_tree->SetBranchAddress("trk2Eta",&trk2eta);
    input_tree->SetBranchAddress("trk2Phi",&trk2phi);
    int nentries = input_tree->GetEntries();
    cout<<"nentries = "<<nentries<<endl;
    for (int ientry = 0; ientry < nentries; ientry++)
    {
        input_tree->GetEntry(ientry);
        TLorentzVector electronPositive;
        electronPositive.SetPtEtaPhiM(trk1pt,trk1eta,trk1phi,0.511e-3);
        TLorentzVector electronNegative;
        electronNegative.SetPtEtaPhiM(trk2pt,trk2eta,trk2phi,0.511e-3);
        TLorentzVector possibleJPsi;
        possibleJPsi.SetPtEtaPhiM(pairpt,paireta,pairphi,3.096);
        PhiHelicityFrame = polarization::get_PhiHelicityFrame(electronPositive,electronNegative,possibleJPsi);
        CosThetaHelicityFrame = polarization::get_CosThetaHelicityFrame(electronPositive,electronNegative,possibleJPsi);
        PhiCollinsSoper = polarization::get_PhiCollinsSoper(electronPositive,electronNegative,possibleJPsi);
        CosThetaCollinsSoper = polarization::get_CosThetaCollinsSoper(electronPositive,electronNegative,possibleJPsi);
        output_tree->Fill();
    }
    // output_tree->Write();
    // output_root->Close();

    return output_tree;

}

//define the function
double polarization::get_CosThetaHelicityFrame(TLorentzVector electronPositive, TLorentzVector electronNegative, TLorentzVector possibleJPsi)
{
    double HalfSqrtSnn   = 6800.;//half of sqrt(S_nn) in GeV
    double MassOfproton = 0.938;//mass of used nuclei
    double MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn - MassOfproton*MassOfproton);//calculate the momentum of the beam
    TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn); // projectile
    TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn); // target
    TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();// beta of the possibleJPsi
    TLorentzVector pelectron1  = electronPositive;// get the positive electron momentum 
    TLorentzVector pelectron2  = electronNegative;// get the negative electron momentum
    TLorentzVector pProjDielectron = pProjCM;// get the projectile beam momentum
    TLorentzVector pTargDielectron = pTargCM;// get the target beam momentum
    // calculate the four-momentum of the dielectron system in the Helicity frame
    pelectron1.Boost(beta);
    pelectron2.Boost(beta);
    pProjDielectron.Boost(beta);
    pTargDielectron.Boost(beta);
    TVector3 zaxis = (possibleJPsi.Vect()).Unit();
    double CosThetaHE = zaxis.Dot((pelectron1.Vect()).Unit());
    return   CosThetaHE;
}

double polarization::get_PhiHelicityFrame(TLorentzVector electronPositive, TLorentzVector electronNegative, TLorentzVector possibleJPsi)
{
    double HalfSqrtSnn   = 6800.;//sqrt(S_nn) in GeV
    double MassOfproton = 0.938;//mass of used nuclei
    double MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn - MassOfproton*MassOfproton );//calculate the momentum of the beam
    TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn); // projectile
    TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn); // target
    TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();// beta of the possibleJPsi
    TLorentzVector pelectron1  = electronPositive;// get the positive electron momentum 
    TLorentzVector pelectron2  = electronNegative;// get the negative electron momentum
    TLorentzVector pProjDielectron = pProjCM;// get the projectile beam momentum
    TLorentzVector pTargDielectron = pTargCM;// get the target beam momentum
    // calculate the four-momentum of the dielectron system in the Helicity frame
    pelectron1.Boost(beta);
    pelectron2.Boost(beta);
    pProjDielectron.Boost(beta);
    pTargDielectron.Boost(beta);
    TVector3 zaxis = (possibleJPsi.Vect()).Unit();
    TVector3 yaxis = ((pProjDielectron.Vect()).Cross(pTargDielectron.Vect())).Unit();
    TVector3 xaxis = (yaxis.Cross(zaxis)).Unit();
    double phiHe = TMath::ATan2((pelectron1.Vect()).Dot(yaxis),(pelectron1.Vect()).Dot(xaxis));
    return   phiHe;
}

double polarization::get_CosThetaCollinsSoper(TLorentzVector electronPositive, TLorentzVector electronNegative, TLorentzVector possibleJPsi)
{
    double HalfSqrtSnn   = 6800.;//sqrt(S_nn) in GeV
    double MassOfproton = 0.938;//mass of used nuclei
    double MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn - MassOfproton*MassOfproton );//calculate the momentum of the beam
    TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn); // projectile
    TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn); // target
    TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();// beta of the possibleJPsi
    TLorentzVector pelectron1  = electronPositive;// get the positive electron momentum 
    TLorentzVector pelectron2  = electronNegative;// get the negative electron momentum
    TLorentzVector pProjDielectron = pProjCM;// get the projectile beam momentum
    TLorentzVector pTargDielectron = pTargCM;// get the target beam momentum
    // calculate the four-momentum of the dielectron system in the Helicity frame
    pelectron1.Boost(beta);
    pelectron2.Boost(beta);
    pProjDielectron.Boost(beta);
    pTargDielectron.Boost(beta);

    TVector3 zaxisCS=(((pProjDielectron.Vect()).Unit())-((pTargDielectron.Vect()).Unit())).Unit();
    Double_t CosThetaCS = zaxisCS.Dot((pelectron1.Vect()).Unit());
    return   CosThetaCS;
}

double polarization::get_PhiCollinsSoper(TLorentzVector electronPositive, TLorentzVector electronNegative, TLorentzVector possibleJPsi)
{
    double HalfSqrtSnn   = 6800.;//sqrt(S_nn) in GeV
    double MassOfproton = 0.938;//mass of used nuclei
    double MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn - MassOfproton*MassOfproton );//calculate the momentum of the beam
    TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn); // projectile
    TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn); // target
    TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();// beta of the possibleJPsi
    TLorentzVector pelectron1  = electronPositive;// get the positive electron momentum 
    TLorentzVector pelectron2  = electronNegative;// get the negative electron momentum
    TLorentzVector pProjDielectron = pProjCM;// get the projectile beam momentum
    TLorentzVector pTargDielectron = pTargCM;// get the target beam momentum
    // calculate the four-momentum of the dielectron system in the Helicity frame
    pelectron1.Boost(beta);
    pelectron2.Boost(beta);
    pProjDielectron.Boost(beta);
    pTargDielectron.Boost(beta);

    TVector3 zaxisCS=(((pProjDielectron.Vect()).Unit())-((pTargDielectron.Vect()).Unit())).Unit();
    TVector3 yaxisCS=(((pProjDielectron.Vect()).Unit()).Cross((pTargDielectron.Vect()).Unit())).Unit();
    TVector3 xaxisCS=(yaxisCS.Cross(zaxisCS)).Unit();
    Double_t phiCS = TMath::ATan2((pelectron1.Vect()).Dot(yaxisCS),(pelectron1.Vect()).Dot(xaxisCS));
    return   phiCS;
}
