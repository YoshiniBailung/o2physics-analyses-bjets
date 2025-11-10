struct bjetTree {
    float mJetpT;
    float mJetEta;
    float mJetPhi;
    int mNTracks;
    int mNSV;
    float mJetMass;
    int mJetFlavor;
  
    std::array<float, 10> mTrackpT;
    std::array<float, 10> mTrackEta;
    std::array<float, 10> mDotProdTrackJet;
    std::array<float, 10> mDotProdTrackJetOverJet;
    std::array<float, 10> mDeltaRJetTrack;
    std::array<float, 10> mSignedIP2D;
    std::array<float, 10> mSignedIP2DSign;
    std::array<float, 10> mSignedIPz;
    std::array<float, 10> mSignedIPzSign;
    std::array<float, 10> mMomFraction;
    std::array<float, 10> mDeltaRTrackVertex;
  
    std::array<float, 10> mSVpT;
    std::array<float, 10> mDeltaRSVJet;
    std::array<float, 10> mSVMass;
    std::array<float, 10> mSVfE;
    std::array<float, 10> mIPXY;
    std::array<float, 10> mCPA;
    std::array<float, 10> mChi2PCA;
    std::array<float, 10> mDispersion;
    std::array<float, 10> mDecayLength2D;
    std::array<float, 10> mDecayLength2DError;
    std::array<float, 10> mDecayLength3D;
    std::array<float, 10> mDecayLength3DError;
  };

void bjettreemerger(){

    bjetTree treeWords;


    TFile *file = TFile::Open("../TreeCreatorTask/AO2D.root");
    TList *klist = file->GetListOfKeys();

    TTree *thash;

    const Int_t nTrees = 4;

    TString param_branch_names[nTrees] = {
        "O2bjetconstit",
        "O2bjetparam",
        "O2bjetsvparam",
        "O2bjettracksparam",
    };

    const Int_t maxtracks = 10;

    std::vector<float> trackpT;

    myFile = new TFile("myFile.root", "RECREATE");
    gDirectory->mkdir("bjet-tree-merger");
    gDirectory->cd("bjet-tree-merger");
    myTree = new TTree("myTree", "Test tree");
    
    myTree->Branch("mJetpT", &treeWords.mJetpT, "mJetpT/F");
    myTree->Branch("mJetEta", &treeWords.mJetEta, "mJetEta/F");
    myTree->Branch("mJetPhi", &treeWords.mJetPhi, "mJetPhi/F");
    myTree->Branch("mNTracks", &treeWords.mNTracks, "mNTracks/I");
    myTree->Branch("mNSV", &treeWords.mNSV, "mNSV/I");
    myTree->Branch("mJetMass", &treeWords.mJetMass, "mJetMass/F");
    myTree->Branch("mJetFlavor", &treeWords.mJetFlavor, "mJetFlavor/I");
    // Track info branches
    myTree->Branch("mTrackpT", &treeWords.mTrackpT, "mTrackpT[10]/F");
    myTree->Branch("mTrackEta", &treeWords.mTrackEta, "mTrackEta[10]/F");
    myTree->Branch("mDotProdTrackJet", &treeWords.mDotProdTrackJet, "mDotProdTrackJet[10]/F");
    myTree->Branch("mDotProdTrackJetOverJet", &treeWords.mDotProdTrackJetOverJet, "mDotProdTrackJetOverJet[10]/F");
    myTree->Branch("mDeltaRJetTrack", &treeWords.mDeltaRJetTrack, "mDeltaRJetTrack[10]/F");
    myTree->Branch("mSignedIP2D", &treeWords.mSignedIP2D, "mSignedIP2D[10]/F");
    myTree->Branch("mSignedIP2DSign", &treeWords.mSignedIP2DSign, "mSignedIP2DSign[10]/F");
    myTree->Branch("mSignedIPz", &treeWords.mSignedIPz, "mSignedIPz[10]/F");
    myTree->Branch("mSignedIPzSign", &treeWords.mSignedIPzSign, "mSignedIPzSign[10]/F");
    myTree->Branch("mMomFraction", &treeWords.mMomFraction, "mMomFraction[10]/F");
    myTree->Branch("mDeltaRTrackVertex", &treeWords.mDeltaRTrackVertex, "mDeltaRTrackVertex[10]/F");
    // Secondary vertex info branches
    myTree->Branch("mSVpT", &treeWords.mSVpT, "mSVpT[10]/F");
    myTree->Branch("mDeltaRSVJet", &treeWords.mDeltaRSVJet, "mDeltaRSVJet[10]/F");
    myTree->Branch("mSVMass", &treeWords.mSVMass, "mSVMass[10]/F");
    myTree->Branch("mSVfE", &treeWords.mSVfE, "mSVfE[10]/F");
    myTree->Branch("mIPXY", &treeWords.mIPXY, "mIPXY[10]/F");
    myTree->Branch("mCPA", &treeWords.mCPA, "mCPA[10]/F");
    myTree->Branch("mChi2PCA", &treeWords.mChi2PCA, "mChi2PCA[10]/F");
    myTree->Branch("mDispersion", &treeWords.mDispersion, "mDispersion[10]/F");
    myTree->Branch("mDecayLength2D", &treeWords.mDecayLength2D, "mDecayLength2D[10]/F");
    myTree->Branch("mDecayLength2DError", &treeWords.mDecayLength2DError, "mDecayLength2DError[10]/F");
    myTree->Branch("mDecayLength3D", &treeWords.mDecayLength3D, "mDecayLength3D[10]/F");
    myTree->Branch("mDecayLength3DError", &treeWords.mDecayLength3DError, "mDecayLength3DError[10]/F");



    for (int i = 0; i < 1; i++){ //Directory iterator 

        TKey *key = (TKey*)klist->At(i);

        TDirectoryFile *dummy = (TDirectoryFile*)key->ReadObj(); 
        
        thash = (TTree*)dummy->Get(Form("%s",param_branch_names[0].Data())); 
        Int_t tsize = thash->GetEntries(); 
        Int_t Indexbjet = 0;

        thash->SetBranchAddress("fIndexbjetParams",&Indexbjet);

        for (int ientry = 0; ientry < tsize; ++ientry){ // bjet param tree iterator

            thash->GetEntry(ientry);
            
            Float_t trackpt = 0;
            Int_t trackbjetindex = 0;

            TTree *ttrack = (TTree*)dummy->Get(Form("%s",param_branch_names[3].Data()));
            Int_t jetsize = ttrack->GetEntries();

            ttrack->SetBranchAddress("fTrackpT",&trackpt);
            ttrack->SetBranchAddress("fIndexbjetParams",&trackbjetindex);

            for (int jentry = 0; jentry < jetsize; ++jentry){ // bjet tracks param tree iterator
                
                ttrack->GetEntry(jentry);
                
                if (Indexbjet != trackbjetindex) continue;

                
                
            }
        }



    }
        
    //     for (int ientry = 0; ientry < tsize; ientry ++){
        
    //     for (int ibr = 0; ibr < param_size; ibr ++) thash->GetBranch(Form("%s",param_branch_names[ibr].Data())); 
        
    // }
}