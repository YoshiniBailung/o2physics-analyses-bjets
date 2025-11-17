void bjetlinearizer(){



    TFile *file = TFile::Open("myFile.root");
    TTree *intree = (TTree*)file->Get("bjet-tree-merger/myTree");

    Float_t mjetpt = 0;
    Int_t mntracks = 0, mjetflavor;
    Float_t msip2d[50], msipz[50];
    Float_t msigmasip2d[50];
    Float_t msigmasip2dsign[50], msigmasipzsign[50];
    Float_t mtrackpt[50];


    intree->SetBranchAddress("mJetpT", &mjetpt);
    intree->SetBranchAddress("mNTracks", &mntracks);
    intree->SetBranchAddress("mJetFlavor", &mjetflavor);
    intree->SetBranchAddress("mSignedIP2D", &msip2d);
    intree->SetBranchAddress("mSignedIP2DSign", &msigmasip2dsign);
    intree->SetBranchAddress("mSignedIPz", &msipz);
    intree->SetBranchAddress("mSignedIPzSign", &msigmasipzsign);
    intree->SetBranchAddress("mTrackpT", &mtrackpt);

    Int_t NJETS = intree->GetEntries();
    // for (int i =0; i < NJETS; i++){
    //     intree->GetEntry(i);

    //     Printf("Jet %d has pt %f with %d tracks", i, mjetpt, mntracks);

    //     for (int j = 0; j < 10 ; j ++){
    //         Printf("JET %d, tracks have [%d]pt %f", i, j, mtrackpt[j]);
    //     }
    // }

    Int_t ntracks = 0, jetflavor = 0;
    Float_t jetpt = 0, sip2d1 = 0, sip2d2 = 0, sip2d3 = 0, sigma2d1 = 0, sigma2d2 = 0, sigma2d3 = 0; //2d coordinates
    Float_t trackpt1 = 0, trackpt2 = 0, trackpt3 = 0; 
    Float_t sip2ds1 = 0, sip2ds2 = 0, sip2ds3 = 0; //2d ip significance
    Float_t sip3d1 = 0, sip3d2 = 0, sip3d3 = 0, sigma3d1 = 0, sigma3d2 = 0, sigma3d3 = 0; //3d coordinates

    TFile *outFile = new TFile("myPartialFile.root", "RECREATE");

    TTree *myTree = new TTree("myPartialTree", "my Partial tree");
    
    myTree->Branch("jetpT", &jetpt, "jetpT/F");
    myTree->Branch("nTracks", &ntracks, "nTracks/I");
    //myTree->Branch("JetFlavor", &jetflavor, "JetFlavor/I");
    //2D IPs
    myTree->Branch("sIP2D1", &sip2d1, "sIP2D1/F");
    myTree->Branch("sIP2D2", &sip2d2, "sIP2D2/F");
    myTree->Branch("sIP2D3", &sip2d3, "sIP2D3/F");
    //2D IP sigmas
    myTree->Branch("sigma2D1", &sigma2d1, "sigma2D1/F");
    myTree->Branch("sigma2D2", &sigma2d2, "sigma2D2/F");
    myTree->Branch("sigma2D3", &sigma2d3, "sigma2D3/F");
    //2D IP significance
    myTree->Branch("sIP2Ds1", &sip2ds1, "sIP2Ds1/F");
    myTree->Branch("sIP2Ds2", &sip2ds2, "sIP2Ds2/F");
    myTree->Branch("sIP2Ds3", &sip2ds3, "sIP2Ds3/F");
    //3D IPs
    myTree->Branch("sIP3D1", &sip3d1, "sIP3D1/F");
    myTree->Branch("sIP3D2", &sip3d2, "sIP3D2/F");
    myTree->Branch("sIP3D3", &sip3d3, "sIP3D3/F");
    //3D IP sigmas
    myTree->Branch("sigma3D1", &sigma3d1, "sigma3D1/F");
    myTree->Branch("sigma3D2", &sigma3d2, "sigma3D2/F");
    myTree->Branch("sigma3D3", &sigma3d3, "sigma3D3/F");

    float weight = 100.;

    for (int ijet = 0; ijet < NJETS ; ijet++){

        intree->GetEntry(ijet);

        //cout<<"Getting Jet "<<ijet<<" with "<<mntracks<<" tracks"<<endl;
        
        vector<vector<float>> SIP2D;
        vector<float> Sigma2D, SIP3D, Sigma3D;
        
        if (mntracks < 3) continue;

        jetpt = mjetpt;
        ntracks = mntracks;
        //jetflavor = mjetflavor;


        for (int itrack = 0; itrack < mntracks; itrack++){
            SIP2D.push_back({msip2d[itrack]*weight,static_cast<float>(itrack)});
            Sigma2D.push_back(msigmasip2dsign[itrack]*weight);

            SIP3D.push_back(msipz[itrack]);
            Sigma3D.push_back(msigmasipzsign[itrack]*weight);

            //Printf("%f - %f - %f - %f", msip2d[itrack], msigmasip2dsign[itrack], msipz[itrack], msigmasipzsign[itrack]);
        }
        //cout<<"********************************************************************************"<<endl;

        std::sort(SIP2D.begin(),
          SIP2D.end(),
          [] (const std::vector<float> &a, const std::vector<float> &b)
          {
              return a[0] > b[0];
          });
        
        auto get = [&](std::vector<float>& v, int idx) {
            return (v.size() > idx ? v[idx] : -999.f);
        };


        sip2d1 = SIP2D[0][0];
        sip2d2 = SIP2D[1][0];
        sip2d3 = SIP2D[2][0];

        int index1 = static_cast<int>(SIP2D[0][1]);
        int index2 = static_cast<int>(SIP2D[1][1]);
        int index3 = static_cast<int>(SIP2D[2][1]);

        sigma2d1 = Sigma2D[index1];
        sigma2d2 = Sigma2D[index2];
        sigma2d3 = Sigma2D[index3];

        sip2ds1 = sip2d1/sigma2d1;
        sip2ds2 = sip2d2/sigma2d2;
        sip2ds3 = sip2d3/sigma2d3;

        sip3d1 = SIP3D[index1];
        sip3d2 = SIP3D[index2];
        sip3d3 = SIP3D[index3];

        sigma3d1 = Sigma3D[index1];
        sigma3d2 = Sigma3D[index2];
        sigma3d3 = Sigma3D[index3];

        SIP2D.clear();
        Sigma2D.clear();
        SIP3D.clear();
        Sigma3D.clear();    

        cout<<"Filling Jet "<<ijet<<" with "<<ntracks<<" tracks"<<endl;

        myTree->Fill();

    }


    myTree->Write();
    outFile->Write();


       


}