{
  gSystem->Load("libFWCoreFWLite.so"); 
  AutoLibraryLoader::enable();
  gSystem->Load("libDataFormatsFWLite.so");
  gSystem->Load("libDataFormatsPatCandidates.so");
  gROOT->ProcessLine("namespace edm {typedef edm::Wrapper<vector<float> > Wrapper<vector<float,allocator<float> > >; }");
  gROOT->ProcessLine("namespace edm {typedef edm::Wrapper<vector<double> > Wrapper<vector<double,allocator<double> > >; }");

  std::cout << "*****************************" << std::endl;
  std::cout << " now compiling fw light macro - ignore all CompositeCandidate errors " << std::endl;
  gSystem->CompileMacro("frameworklightana.C","k");
  frameworklightana();
}
