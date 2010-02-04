{
  gSystem->Load("libFWCoreFWLite.so"); 
  AutoLibraryLoader::enable();
  gSystem->Load("libDataFormatsFWLite.so");
//   gROOT->ProcessLine("namespace edm {typedef edm::Wrapper<vector<float> > Wrapper<vector<float,allocator<float> > >; }");
//    gROOT->ProcessLine("namespace edm {typedef edm::Wrapper<vector<double> > Wrapper<vector<double,allocator<double> > >; }");

  gSystem->CompileMacro("makeplot_fwlite.C","k");
  gSystem->CompileMacro("makeplot_fwlite.C","k");

  makeplot_fwlite();
  gSystem->Exit(1);
}
