#include <TFile.h>
#include <TH1.h>

using namespace std;
TList* get_List(TFile*, string);

TList* get_List(TFile *file, string HashListname,string Listname)
{
    THashList *HashList = (THashList*)file->Get(HashListname.c_str());
    TList *List = (TList*)HashList->FindObject(Listname.c_str());
    return List;
}

