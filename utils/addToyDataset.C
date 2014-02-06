// 
void addToyDataset(const char * workspaceName, const char * toyFileName, const char * toyName, const char* outputFileName) {
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");

  // the workspace with the required POIs
  TFile workspace(workspaceName);
  RooWorkspace * myWS = workspace.Get("w");

  // the file with the toy 
  TFile fileWithToy(toyFileName);
  RooAbsData * toyDataset = (RooAbsData *) fileWithToy.Get(Form("toys/%s",toyName));

  if (!toyDataset) {
    cout << "Error: toyDataset " << Form("toys/%s",toyName) 
         << " not found in file " << toyFileName << endl;
    return;
  }
  toyDataset->SetName(Form("toys/%s",toyName));
  // import (the two for completness
  myWS->import(*toyDataset);
  
  // write on a new workspace file
  myWS->writeToFile(outputFileName);
}
