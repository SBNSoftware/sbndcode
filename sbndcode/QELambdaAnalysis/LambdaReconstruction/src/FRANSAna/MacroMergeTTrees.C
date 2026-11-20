
// run as
//  root
// .L MacroMergeTTrees.C
// MacroMergeTTrees("CCQE/", "BG") (directory with th trees to merge, and the label)
// MacroMergeTTrees("VOLambda/", "S")
void MacroMergeTTrees(std::string inputName, std::string tree_dirname, std::string tree_name = "FRAMSTree"){

  std::vector<std::string> file_names;

  TString inputNameTSt(inputName.c_str());

  if( inputNameTSt.EndsWith(".list") ){
    std::ifstream infile(inputName);
    while(!infile.eof()){
      std::string fname;
      infile>>fname;
      file_names.push_back(fname);
    }
  }
  else{
    const char *ext=".root";

    TSystemDirectory dir(inputName.c_str(), inputName.c_str());
    TList *files = dir.GetListOfFiles();
    if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
        fname = file->GetName();
        if (!file->IsDirectory() && fname.EndsWith(ext)) {
            cout << fname.Data() << endl;
            file_names.push_back( (inputName+fname.Data()) );
        }
      }
    }
  }
  

  std::string outName = "Merged"+tree_name+"_"+tree_dirname+".root";
  TFile *fout = new TFile(outName.c_str(),"RECREATE");
  fout->cd();
  gDirectory->mkdir(tree_dirname.c_str());
  gDirectory->cd(tree_dirname.c_str());

  TTree *mergedtree;

  TList *tree_list = new TList;

  for(size_t k=0; k<file_names.size(); k++){
    std::cout<<k<<" "<<file_names[k]<<std::endl;
    TFile *ft = new TFile(file_names[0].c_str(),"READ");
    TTree *t1 = (TTree *)ft->Get((tree_dirname+"/"+tree_name).c_str());
    //if(k>20) continue;
    tree_list->Add(t1);
  }

  fout->cd();
  gDirectory->cd(tree_dirname.c_str());
  mergedtree = TTree::MergeTrees(tree_list);
  mergedtree->SetName(tree_name.c_str());
  mergedtree->Write();


  fout->Close();

  return 0;
}
