
// Just the order in which the are placed in the raw data.
int G_detector_angles[12]={30, 54, 78, 102, 126, 150, 210, 234, 258, 282, 306, 330};


// Detector timing offsets. Very importnt!
float G_pmt_offsets[12][2] = {{-17.2,-30.4},{-13.9,-24.9},{-18.5,-33.2},{-29.6,-28.0},{-14.9,-25.2},{-18.0,-25.9},{-17.7,-27.9},{-14.5,-27.8},{-1.2,-22.8},{-13.9,-28.2},{-15.4,-27.9},{-17.0,-33.0}};

// Calculate opening angle
TVector3 v1,v2;
float vector_angle(float x1,float y1,float z1,float x2,float y2,float z2){
    v1.SetXYZ(x1,y1,z1);
    v2.SetXYZ(x2,y2,z2);
    return 180/3.1415*v1.Angle(v2);
}

// Belw are lists of run numbers for different targets.
int D2Os[7] = {6541, 6551, 6565, 6586, 6607, 6621, 6642};
//
int DUs[56] = {6531, 6532, 6533, 6534, 6535, 6536, 6537, 6538, 6542, 6543, 6544,
          6545, 6546, 6547, 6548, 6549, 6553, 6554, 6556, 6557, 6558, 6559,
          6560, 6561, 6568, 6570, 6571, 6572, 6573, 6587, 6588, 6591, 6592,
          6600, 6601, 6603, 6608, 6609, 6610, 6611, 6612, 6613, 6614, 6615,
          6616, 6617, 6625, 6626, 6627, 6628, 6629, 6630, 6634, 6635, 6636,
          6640};

int Ths[5] = {6622, 6623, 6624, 6638, 6639};

int Als[7] = {6604, 6552, 6569, 6598, 6620, 6632, 6643};

//int Cosmics[4] = {6562,6585,6631,6641};
int Cosmics[2] = {6550, 6605};

int U233s[5] = {6644, 6645, 6646, 6646, 6647};

int Pus[3] = {6648, 6649, 6650};


// This structre contains all the usful informtion of a single hit.
struct Hit{
    float tof=0;
    float top=0;
    float bot=0;
    float pmt_diff=0;
    int det=0;
    float x=0;
    float y=0;
    float z=0;
    float theta_abs=0;
    float erg=0;
    int ForwardDet = 0;
};

void Do_Something_with_histogram(Hit hit1, Hit hit2, TH1F * hist){
    float opening_angle = vector_angle(hit1.x, hit1.y, hit1.z, hit2.x, hit2.y, hit2.z);
    hist->Fill(opening_angle);
}

void test(){
    int TDC1190[10][129]; // Container for the raw data from the DAQ

    // Make this a correct path for your system
    TFile * file=new TFile("/Volumes/JeffMacEx/2nCorrData/Aug_2017/r6634.root");
    TTree * daq_tree = (TTree*)file->Get("FF");

    TLeaf * br_TDC1190 = daq_tree->GetBranch("TDC")->GetLeaf("TDC1190"); // To set branch address
    br_TDC1190->SetAddress(&TDC1190);

    TH1F * my_hist = new TH1F("","",180,20,180);

    for(int ientry = 1; ientry < daq_tree->GetEntries(); ientry++){
        if(ientry%10000 == 0)cout<<"Entry: "<<ientry<<" outta: "<<daq_tree->GetEntries()<<endl;

        std::vector<Hit> hits = std::vector<Hit>();

        daq_tree->GetEntry(ientry);

        for (int i=1;i<10;i++){
            int top_index = 2*i+1;  // index for the array, e.g TDC1190[1][j]. Odd numbers are top PMT, evens are bottom.
            int bot_index = top_index+1;

            // factor of 0.1 for converting to ns

            float pmt_top = 0.1*TDC1190[1][top_index];
            float pmt_bot = 0.1*TDC1190[1][bot_index];
            pmt_top += TDC1190[1][top_index]!=0   ? G_pmt_offsets[i][0] : 0.0;
            pmt_bot += TDC1190[1][bot_index]!=0   ? G_pmt_offsets[i][1] : 0.0;

            float trig=0.1*TDC1190[1][32];

            if (pmt_top!=0  && pmt_bot!=0){

                if (trig!=0){

                    float tof= 0.5* (pmt_top + pmt_bot)-trig;

                    float pmt_diff= pmt_bot-pmt_top; // This order of subtraction is the most practical.

                    if (abs(pmt_diff)<40){
                        Hit * hit = new Hit(); // deleted in set_all_hits, or in the if statement right below.
                        hit->tof=tof;

                        /*
                        //use this to cut on neutrons
                        if (!(hit->tof>45 && hit->tof<145)){
                            delete hit;
                            continue;
                        }
                        */

                        hit->pmt_diff = pmt_diff;
                        hit->top = pmt_top;
                        hit->bot  = pmt_bot;

                        hit->det  =G_detector_angles[i];

                        hit->x = 125*cos(hit->det/180.*3.1415);
                        hit->y = -(125*sin(hit->det/180.*3.1415)); // Is negative to make R-handed coord sys with +z going upward.
                        hit->z = (hit->bot-hit->top)*6.5;

                        hit->theta_abs = TMath::ACos(hit->x/sqrt(hit->x*hit->x+hit->y*hit->y+hit->z*hit->z));

                        if (abs(hit->tof)>1) hit->erg=8127.0/hit->tof/hit->tof;

                        hits.push_back(*hit);
                        if (hits.size()==2){
                            // The function below does whatever you want with a pair of coincident hits
                            Do_Something_with_histogram(hits[0], hits[1], my_hist);
                            hits.clear();
                    }
                }
            }
        }
    }
}

    my_hist->Draw();
}

//    char path[100]=;
//    _F=
//   daq_tree.Add()

//
//
//    Hit * Pulse::__get_Hit__(int TDC1190[][129], int i, float trig){
//            if (i==0 || i==11){
//                cout<<"Forward detectors not supported in this function"<<endl;
//                return NULL
//            }
//            int top_index = 2*i+1;  // index for the array, e.g TDC1190[1][j]. Odd numbers are top PMT, evens are bottom.
//            int bot_index = top_index+1;
//
//
//            float pmt_top = 0.1*TDC1190[1][top_index];
//            float pmt_bot = 0.1*TDC1190[1][bot_index];
//
//
//            if (pmt_top!=0  && pmt_bot!=0){
//
//                if (trig!=0){
//
//                    float tof= 0.5* (pmt_top + pmt_bot)-trig;
//
//                    float pmt_diff= pmt_bot-pmt_top; // This order of subtraction is the most practical.
//
//                    if (abs(pmt_diff)<40){
//                        Hit * hit = new Hit(); // deleted in set_all_hits, or in the if statement right below.
//                        hit->tof=tof;
//
//                        hit->pmt_diff = pmt_diff;
//                        hit->top = pmt_top;
//                        hit->bot  = pmt_bot;
//
//                        hit->det  =G_detector_angles[i];
//
//                        hit->x = 125*cos(hit->det/180.*3.1415);
//                        hit->y = -(125*sin(hit->det/180.*3.1415)); // Is negative to make R-handed coord sys with +z going upward.
//                        hit->z = (hit->bot-hit->top)*6.5;
//
//                        hit->theta_abs = TMath::ACos(hit->x/sqrt(hit->x*hit->x+hit->y*hit->y+hit->z*hit->z));
//
//                        if (abs(hit->tof)>1) hit->erg=8127.0/hit->tof/hit->tof;
//
//                        return hit;
//                    }
//                }
//            }
//            return NULL;
//    }
