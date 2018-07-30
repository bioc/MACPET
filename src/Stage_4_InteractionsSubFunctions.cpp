#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/MatrixAccessor.hpp>
#include <bigmemory/BigMatrix.h>
// [[Rcpp::interfaces(r, cpp)]]

//####################################################
// Function to decide the new peak summit based on the merges:
//####################################################
//[[Rcpp::export]]
SEXP Get_NewPeakSummit_fun_Rcpp(Rcpp::NumericVector const &queryHits,
                                Rcpp::NumericVector const &subjectHits,
                                Rcpp::NumericVector const &PeakSummit,
                                Rcpp::NumericVector const &FDR,
                                int const &Noverlaps, int const &NPeaksMerged){
    //---------------------
    // Initialize:
    //---------------------
    Rcpp::NumericVector PeakSummitNew(NPeaksMerged);
    //---------------------
    // loop on data:
    //---------------------
    int outer = 0;
    int MergedPeakCount = 0;
    while(outer<Noverlaps){
        int SubHitouter=subjectHits[outer]-1;//will also take self
        double FDRouter=FDR[SubHitouter];
        double outerCounts=0;
        double PeakSummit_opt;
        //---------------------
        // loop inner:
        //---------------------
        for(int inner=outer;inner<Noverlaps;inner++){
            if(queryHits[outer]==queryHits[inner]){
                outerCounts+=1;
                int SubHitinner=subjectHits[inner]-1;//will also take self
                double FDRinner=FDR[SubHitinner];
                if(FDRinner<=FDRouter){
                    FDRouter=FDRinner;
                    PeakSummit_opt=PeakSummit[SubHitinner];
                }
            }else{
                break;
            }
        }
        PeakSummitNew[MergedPeakCount]=PeakSummit_opt;
        MergedPeakCount+=1;
        outer+=outerCounts;
    }
    return Rcpp::wrap(PeakSummitNew);
}
// done
//####################################################
// Funtions for classifing the interaction PETs according to their densities
// And getting counts data
//####################################################
//---------------------
// Function to create the PETsInfoMat which holds the classification of each interaction PET in the data
// into the peaks. PETs which fall into one peak are removed
//---------------------
//[[Rcpp::export]]
SEXP Get_PETsInfoMat_fun_Rcpp(Rcpp::NumericVector const &VEC_query,//query(PET ID)
                              Rcpp::NumericVector const &VEC_Type,//If Anchor 1 or 2 for the PETs
                              Rcpp::NumericVector const &VEC_Tag,//The midpoint of the Anchor (for the PET)
                              Rcpp::NumericVector const &VEC_LID,//The Line ID in PeakData(R index)
                              Rcpp::NumericVector const &VEC_PeakSummit,//the midtpoint of the peak
                              int const &NGlobalInterPETs,//the total PETs to be classified
                              int const &NIntTagsloop){//total interaction tags to classify
    // ---------------------------------------
    // Function declaration:
    // ---------------------------------------
    // function for the data frame substitution based on Types variable:
    void Subset_Query_fun_Rcpp(int const &iq_outer, Rcpp::NumericVector &Type1_LID,
                               Rcpp::NumericVector &Type2_LID,Rcpp::NumericVector &Type1_PeakSummit,
                               Rcpp::NumericVector &Type2_PeakSummit, int &Type1_n, int &Type2_n,
                               double &Type1_Tag, double &Type2_Tag, int const &NIntTagsloop,
                               Rcpp::NumericVector const &VEC_query, Rcpp::NumericVector const &VEC_Type,
                               Rcpp::NumericVector const &VEC_LID, Rcpp::NumericVector const &VEC_PeakSummit,
                               Rcpp::NumericVector const &VEC_Tag);
    // Function to classify the tag to the nearest peak if it overlaps with more than two(never happens)
    int Get_TagPeakSummitDist_fun_Rcpp(int const &Type_n,
                                       double  &Type_Tag, Rcpp::NumericVector  const &Type_PeakSummit,
                                       Rcpp::NumericVector const &Type_LID);
    // function for updating the Tag counts in each PBS:
    void Save_TagClassification_fun_Rcpp(int const &Type1_n, int const &Type2_n,
                                         int const &Type1_LID_max,int const &Type2_LID_max,
                                         Rcpp::NumericMatrix &PETsInfoMat,
                                         int &PETscan, int &NPETspassing);
    // ---------------------------------------
    // --------create the PETsInfoMat
    // ---------------------------------------
    // Col: 1-PBSi,2-PBSj(>i),col-3 1 if passing 0 if not
    Rcpp::NumericMatrix PETsInfoMat(NGlobalInterPETs,3);
    int PETscan=0;//counter for the PETs scanned for the PETsInfoMat entries
    int NPETspassing=0;//counter for the total PETs interacting, for redusing the matrix afterwards
    // ---------------------------------------
    // start looping the interaction PETs for classification:
    // In each loop you will subset the input based on the same query name
    // Then for each query, which means a single TAG, you will classify it.
    // If it overlaps with more than one peaks, choose the nearest
    // ---------------------------------------
    int iq_outer=0;//the outer query index row
    // loop on the main query
    while(iq_outer<NIntTagsloop){
        // ---------------------------------------
        // Get the split of the subset in types 1 and 2:
        // ---------------------------------------
        // Type 1 referes to the left-most tag of a PET, Type2 to the rightmost tag
        Rcpp::NumericVector Type1_LID,Type2_LID,Type1_PeakSummit,Type2_PeakSummit;//LID and midtpoint (peaks) positions
        int Type1_n=0,Type2_n=0;//Hits per interaction PET anchor/Tag
        double Type1_Tag,Type2_Tag;//the tag positions
        Subset_Query_fun_Rcpp(iq_outer, Type1_LID,Type2_LID,Type1_PeakSummit, Type2_PeakSummit,Type1_n, Type2_n,
                              Type1_Tag, Type2_Tag, NIntTagsloop,VEC_query, VEC_Type,
                              VEC_LID, VEC_PeakSummit,VEC_Tag);
        // ---------------------------------------
        // For each Tag, if it overlaps with more than one Peak,
        // choose the nearest one.If it does not overlap with any then
        // ---------------------------------------
        // Anchor/TAG 1 of interaction PET:
        int Type1_LID_max=Get_TagPeakSummitDist_fun_Rcpp(Type1_n,Type1_Tag, Type1_PeakSummit,
                                                         Type1_LID);
        // Anchor/TAG 2 of interaction PET:
        int Type2_LID_max=Get_TagPeakSummitDist_fun_Rcpp(Type2_n,Type2_Tag, Type2_PeakSummit,
                                                         Type2_LID);
        // ---------------------------------------
        // Now check the Peaks (LID_max)
        // To save the counts if any.
        // ---------------------------------------
        Save_TagClassification_fun_Rcpp(Type1_n,Type2_n,
                                        Type1_LID_max,Type2_LID_max,
                                        PETsInfoMat,PETscan,NPETspassing);
        // ---------------------------------------
        // increase the iq_outer counter to go to the next interaction PET:
        // NOTE: it is NOT possible that both Type1_n and Type2_n
        // are zero! Because I know from R when I found
        // which interaction PETs overlap with the Peak streams.
        // So iq_outer is always increasing.
        // ---------------------------------------
        iq_outer+=Type1_n+Type2_n;
        PETscan+=1;
    }
    // ---------------------------------------
    // Reduce the matrix to return only the passing PETs:
    // ---------------------------------------
    if(NPETspassing==0){
        return Rcpp::wrap(R_NilValue);//for error handling
    }
    // create new data:
    Rcpp::NumericMatrix PETsInfoMatPassing(NPETspassing,2);
    int Passscan=0;
    // loop and fill in:
    for(int i=0;i<NGlobalInterPETs;i++){
        if(PETsInfoMat(i,2)==1){
            // then passing:
            PETsInfoMatPassing(Passscan,0)=PETsInfoMat(i,0);
            PETsInfoMatPassing(Passscan,1)=PETsInfoMat(i,1);
            Passscan+=1;
        }
        // early break:
        if(Passscan==NPETspassing) break;
    }
    return Rcpp::wrap(PETsInfoMatPassing);
}
// done
//---------------------
// Function for subseting the data.frame sub into two categories based on their Type vector
//---------------------
void Subset_Query_fun_Rcpp(int const &iq_outer, Rcpp::NumericVector &Type1_LID,
                           Rcpp::NumericVector &Type2_LID,Rcpp::NumericVector &Type1_PeakSummit,
                           Rcpp::NumericVector &Type2_PeakSummit, int &Type1_n, int &Type2_n,
                           double &Type1_Tag, double &Type2_Tag, int const &NIntTagsloop,
                           Rcpp::NumericVector const &VEC_query, Rcpp::NumericVector const &VEC_Type,
                           Rcpp::NumericVector const &VEC_LID, Rcpp::NumericVector const &VEC_PeakSummit,
                           Rcpp::NumericVector const &VEC_Tag){
    // ---------------------------------------
    // loop to subset all rows with the same query:
    // This corresponds to one single interaction PET:
    // ---------------------------------------
    for(int iq_inner=iq_outer;iq_inner<NIntTagsloop;iq_inner++){
        // check queries:
        if(VEC_query[iq_inner]==VEC_query[iq_outer]){
            // then same query/TAG.
            if(VEC_Type[iq_inner]==1){
                // Tag comes from first Anchor of the interaction PET:
                Type1_n+=1;
                Type1_LID.push_back(VEC_LID[iq_inner]);
                Type1_PeakSummit.push_back(VEC_PeakSummit[iq_inner]);
                Type1_Tag=VEC_Tag[iq_inner];
            }else{
                // Tag comes from second Anchor of the interaction PET:
                Type2_n+=1;
                Type2_LID.push_back(VEC_LID[iq_inner]);
                Type2_PeakSummit.push_back(VEC_PeakSummit[iq_inner]);
                Type2_Tag=VEC_Tag[iq_inner];
            }
        }else{
            // entered a different query
            // which means different interaction PET so break
            break;
        }
    }
}
// done
//---------------------
// Function to get the mearest peak to the tag if neccessery
//---------------------
int Get_TagPeakSummitDist_fun_Rcpp(int const &Type_n,
                                   double  &Type_Tag, Rcpp::NumericVector const &Type_PeakSummit,
                                   Rcpp::NumericVector const &Type_LID){
    // ---------------------------------------
    // Find the probabilities:
    // ---------------------------------------
    int Type_LID_max=-1;//in case non is found, that is Type_n=0
    double Type_TagMidtDist_max=INFINITY;
    // get all the Peak probabilities the corresponding tag is classified to and return the max
    for(int t=0;t<Type_n;t++){
        // get the probability:
        double Type_TagMidtDist=std::abs(Type_Tag-Type_PeakSummit[t]);
        // check with min:
        if(Type_TagMidtDist<Type_TagMidtDist_max){
            // replace:
            Type_LID_max=Type_LID[t];
            Type_TagMidtDist_max=Type_TagMidtDist;
        }
    }
    return(Type_LID_max);
}
// done
//---------------------
// Function for classifing the optimal PET and counting tags and saving into the matrix
//---------------------
void Save_TagClassification_fun_Rcpp(int const &Type1_n, int const &Type2_n,
                                     int const &Type1_LID_max,int const &Type2_LID_max,
                                     Rcpp::NumericMatrix &PETsInfoMat,
                                     int &PETscan, int &NPETspassing){
    //---------------------
    //take cases and classify the current PET:
    //---------------------
    if(Type1_n!=0&&Type2_n!=0&&Type1_LID_max!=Type2_LID_max){
        // Then the current PET is a connection PET.
        // The Type1_LID_max!=Type2_LID_max ensures that PET is not connecting the same PBS
        //-----------
        // Then update the interaction data:
        //-----------
        NPETspassing+=1;
        PETsInfoMat(PETscan,2)=1;//mark as passing
        // here you will save in the matrices using the pointers.
        // save only on the upper triangle of each matrix.
        if(Type1_LID_max<Type2_LID_max){
            // the order of saving is correct
            PETsInfoMat(PETscan,0)=Type1_LID_max;
            PETsInfoMat(PETscan,1)=Type2_LID_max;
        }else{
            // you need to change the order of the inputs
            // for saving on the upper diagonal
            PETsInfoMat(PETscan,0)=Type2_LID_max;
            PETsInfoMat(PETscan,1)=Type1_LID_max;
        }
    }//dont save it if it falls on the same peak
    // In any other case, one or both of the PET tags does not overlapp with any Peak
    // Those cases are skipped. The PET entry will be empty for PBS-i/j
}
// done
//####################################################
// Function for initializing the InteractionInfMat and the Network:
//####################################################
//---------------------
//Main function for initializing the matrix InteractionInfMat
//---------------------
// [[Rcpp::export]]
SEXP Initiate_InteractionInfMat_fun_Rcpp(Rcpp::NumericMatrix &InteractionInfMat,// the InteractionInfMat,
                                         Rcpp::NumericMatrix &InteractionInfo,//the InteractionInfo
                                         int &NPeaksInvolved,//total of peaks involved
                                         int &NInteractions){//total potential interactions in data
    // ---------------------------------------
    // initiate output:
    // ---------------------------------------
    Rcpp::NumericVector AllInteIndeces(NInteractions);//for looping p-values
    // Vector for the total tags in each peak:
    Rcpp::NumericVector NiNjMat(NPeaksInvolved);
    // ---------------------------------------
    // scan the InteractionInfo and fill the other matrices
    // ---------------------------------------
    for(int i=0;i<NInteractions;i++){
        // ---------------------------------------
        // Save the InteractionInfMat:
        // ---------------------------------------
        InteractionInfMat(i,0)=InteractionInfo(i,0);//PBS i in R for PeaksData
        InteractionInfMat(i,1)=InteractionInfo(i,1);//PBS j in R for PeaksData
        InteractionInfMat(i,2)=InteractionInfo(i,5);//R index for bi products in AdjMat
        InteractionInfMat(i,3)=InteractionInfo(i,6);//R index for bi products in AdjMat
        InteractionInfMat(i,4)=InteractionInfo(i,5);//R index for AdjMat
        InteractionInfMat(i,5)=InteractionInfo(i,6);//R index for AdjMat
        InteractionInfMat(i,6)=InteractionInfo(i,7);//R index for NiNjMat
        InteractionInfMat(i,7)=InteractionInfo(i,8);//R index for NiNjMat
        InteractionInfMat(i,8)=NA_REAL;//p-value
        InteractionInfMat(i,9)=NA_REAL;//FDR
        InteractionInfMat(i,10)=InteractionInfo(i,2);//nij PETs
        InteractionInfMat(i,11)=NA_REAL;//QCell for Vij
        InteractionInfMat(i,13)=InteractionInfo(i,3);//Chrom12ID
        InteractionInfMat(i,14)=InteractionInfo(i,4);//IntraID
        // ---------------------------------------
        // Take the line ids for NiNjMat
        // ---------------------------------------
        int Node_i=InteractionInfo(i,7)-1;//c++ index
        int Node_j=InteractionInfo(i,8)-1;//c++ index
        // ---------------------------------------
        // save the NiNjMat:
        // ---------------------------------------
        NiNjMat[Node_i]+=InteractionInfo(i,2);
        NiNjMat[Node_j]+=InteractionInfo(i,2);
        // interaction vector indicator:
        AllInteIndeces[i]=i+1;//the LID in the significant matrix (R index)
    }
    Rcpp::List Indeces_Vij=Rcpp::List::create(Rcpp::Named("AllInteIndeces")=AllInteIndeces,
                                              Rcpp::Named("NiNjMat")=NiNjMat);
    return Rcpp::wrap(Indeces_Vij);
}
// done
//---------------------
//Main function for initializing the Network
//---------------------
// [[Rcpp::export]]
SEXP Initiate_GenomeMap_fun_Rcpp(int const &NPeaksInvolved_Net,//total peaks involved in network
                                 Rcpp::NumericVector const &AdjNode_Net,//the Adj node ids
                                 Rcpp::NumericVector const &PBS_Net,//IDS in the PeakSummit
                                 Rcpp::NumericVector const &PeakSummit_Net,//the peak summit
                                 int const &Chrom12ID_Net){//the network ID
    // ---------------------------------------
    // Initiate network list:
    // Each element is the Node_i/j and it is a list
    // with names weights/edges with the other nodes
    // and the distance to the other nodes.
    // By the end sort the edges/weights in increasing weights order
    // ---------------------------------------
    Rcpp::List Network(NPeaksInvolved_Net);
    // ---------------------------------------
    // --------loop and fill in the matrix for every entry columnwise
    // ---------------------------------------
    double PrintThres=1;//threshold for printing
    for(int i=0;i<NPeaksInvolved_Net;i++){
        // printing:
        double PrintIntdex=std::floor((i+1.0)/NPeaksInvolved_Net*100);
        if(PrintIntdex==PrintThres){
            PrintThres+=1;
            Rcout<<"|----Initiating the Genome network <-"<<Chrom12ID_Net<<"->: "<<PrintIntdex<<"% completed ----|\r";
        }
        // ---------------------------------------
        // take info for the i:
        // ---------------------------------------
        int PBS_i=PBS_Net[i]-1;//c++ for peaksummit
        int AdjNode_i=AdjNode_Net[i]-1;//c++ index for the network
        // Initiate positions:
        int ID_left=-1, ID_right=-1;
        double Dist_left=INFINITY, Dist_right=-INFINITY; //left-i-right
        // loop through j:
        for(int j=0;j<NPeaksInvolved_Net;j++){
            // check if skipping same peak:
            if(i==j) continue;
            // ---------------------------------------
            // Take the PBS_j from PeaksData
            // ---------------------------------------
            int PBS_j=PBS_Net[j]-1;//c++ for peaksummit
            int AdjNode_j=AdjNode_Net[j]-1;//c++ index for the network
            // ---------------------------------------
            // get distance:
            // ---------------------------------------
            double Dij=PeakSummit_Net[PBS_i]-PeakSummit_Net[PBS_j];//need to be negative too
            // ---------------------------------------
            // check left:
            // ---------------------------------------
            if(Dij>0&&Dij<Dist_left){
                Dist_left=Dij;
                ID_left=AdjNode_j+1;//R index
            }
            // ---------------------------------------
            // check right:
            // ---------------------------------------
            if(Dij<0&&Dij>Dist_right){
                Dist_right=Dij;
                ID_right=AdjNode_j+1;//R index
            }
        }
        // ---------------------------------------
        // save the i entry in the network:
        // ---------------------------------------
        // Create network list for i:
        Rcpp::List Network_i=Rcpp::List::create(Rcpp::Named("edges")=NA_REAL,
                                                Rcpp::Named("weights")=NA_REAL);
        if(ID_left!=-1&&ID_right!=-1){
            Network_i[0]=Rcpp::NumericVector::create(ID_left,ID_right);//R indeces
            Network_i[1]=Rcpp::NumericVector::create(std::abs(Dist_left),std::abs(Dist_right));
        }else if(ID_left!=-1){
            Network_i[0]=Rcpp::NumericVector::create(ID_left);//R indeces
            Network_i[1]=Rcpp::NumericVector::create(std::abs(Dist_left));
        }else if(ID_right!=-1){
            Network_i[0]=Rcpp::NumericVector::create(ID_right);//R indeces
            Network_i[1]=Rcpp::NumericVector::create(std::abs(Dist_right));
        }
        // save to global:
        Network[AdjNode_i]=Network_i;
    }
    return Rcpp::wrap(Network);
}
// done
//---------------------
//Function for getting the mapping index from a matrix to the vector index
//---------------------
// [[Rcpp::export]]
int Get_VectPosIndex_fun_Rcpp(int &NPeaksInvolved, int &Nadj, int &Adj_i, int &Adj_j){
    // The general formula for N=NPeaksInvolved, i=Adj_i, j=Adj_j is:
    // N(N-1)/2-(N-i)(N-i-1)/2+j-i-1, for i<j
    // which is: Nadj-(N-i)(N-i-1)/2+j-i-1.
    // First transform the i and j to c++ indeces:
    // NOTE: dont swap them because they are passed as refence!
    // compute the output:
    int VectID;
    if(Adj_i<Adj_j){
        // correct order
        VectID=Nadj-(NPeaksInvolved-Adj_i)*(NPeaksInvolved-Adj_i-1)/2+Adj_j-Adj_i-1;
    }else if(Adj_i>Adj_j){
        // change orders
        VectID=Nadj-(NPeaksInvolved-Adj_j)*(NPeaksInvolved-Adj_j-1)/2+Adj_i-Adj_j-1;
    }else{
        VectID=-1;//the diagonal is not in, this never happens
    }
    return VectID;
}
// done
//---------------------
//Vectorized function of the above
//---------------------
// [[Rcpp::export]]
SEXP Get_VectPosIndex_Vectorized_fun_Rcpp(int &NPeaksInvolved, int &Nadj, Rcpp::NumericVector const &Adj_i_vect,
                                          Rcpp::NumericVector const &Adj_j_vect){
    //---------------------
    // Function declaration:
    //---------------------
    int Get_VectPosIndex_fun_Rcpp(int &NPeaksInvolved, int &Nadj, int &Adj_i, int &Adj_j);
    //---------------------
    // Loop in each element and calculate the index:
    //---------------------
    int Size=Adj_i_vect.size();
    Rcpp::NumericVector Adj_ij(Size);
    for(int i=0;i<Size;i++){
        // take indeces and convert to c++
        int Adj_i=Adj_i_vect[i]-1;
        int Adj_j=Adj_j_vect[i]-1;
        // convert and save:
        Adj_ij[i]=Get_VectPosIndex_fun_Rcpp(NPeaksInvolved,Nadj,Adj_i,Adj_j)+1;//convert to R index

    }
    return Rcpp::wrap(Adj_ij);
}
// done
//####################################################
//Functions for the expected PETs under H0
//####################################################
//---------------------
// Main function to call the shortest path for global paths from src(source) node using dijkstra algorithm
// The function uses Set for the next node to check
//---------------------
//[[Rcpp::export]]
SEXP Dijkstra_GSP_fun_Rcpp(int &src,//source node in R index
                           Rcpp::List const &Network,//the network to use
                           int const &NPeaksInvolved){//the total Peaks
    //---------------------
    // Initialize:
    //---------------------
    Rcpp::NumericVector GlobalNodesDist(NPeaksInvolved,INFINITY);//dist of src node to each other node
    // Create a set to store vertices that are being prerocessed:
    // first pair position is the distance from src, second in the node (R index)
    std::set< std::pair<double, int> > SetPairs;
    //---------------------
    // set distance of src to 0 since it is the source node:
    // And also set the src pair in the SetPairs:
    //---------------------
    GlobalNodesDist[src-1]=0;
    SetPairs.insert(std::make_pair(0.0, src));//node in R index
    //---------------------
    // Loop for the SP:
    //---------------------
    while(!SetPairs.empty()){
        //---------------------
        // Take node with the shortest distance:
        // The first vertex in Set is the minimum distance vertex, extract it from set.
        //---------------------
        std::pair<double, int> CurPair = *(SetPairs.begin());//get first pair
        SetPairs.erase(SetPairs.begin());//erase it
        //---------------------
        // Take the node of that pair:
        //---------------------
        // vertex label is stored in second of pair (it
        // has to be done this way to keep the vertices
        // sorted distance (distance must be first item
        // in pair)
        int CurNode = CurPair.second;//c++ index
        //---------------------
        // take CurNode info:
        //---------------------
        Rcpp::List Network_CurNode=Network[CurNode-1];
        Rcpp::NumericVector edges_CurNode=Network_CurNode["edges"];
        Rcpp::NumericVector weights_CurNode=Network_CurNode["weights"];
        int Paths_CurNode_tot=edges_CurNode.size();//total paths from the CurNode
        //---------------------
        //Scan the children of the CurNode:
        //---------------------
        for(int v=0;v<Paths_CurNode_tot;v++){
            int edge_v=edges_CurNode[v];
            double weight_v=weights_CurNode[v];
            // Compute the distance and save if smaller:
            double dist_to_v=GlobalNodesDist[CurNode-1]+weight_v;
            if(dist_to_v<GlobalNodesDist[edge_v-1]){
                //---------------------
                // If distance of v is not INF then it must be in
                // our set, so removing it and inserting again
                // with updated less distance.
                // Note : We extract only those vertices from Set
                // for which distance is finalized. So for them,
                // we would never reach here.
                //---------------------
                if (!std::isinf(GlobalNodesDist[edge_v-1])){
                    SetPairs.erase(SetPairs.find(std::make_pair(GlobalNodesDist[edge_v-1], edge_v)));
                }
                // save to shortest path distance:
                GlobalNodesDist[edge_v-1]=dist_to_v;
                // update set:
                SetPairs.insert(std::make_pair(GlobalNodesDist[edge_v-1], edge_v));
            }
        }
    }
    return Rcpp::wrap(GlobalNodesDist);
}
// done
//---------------------
// Function for saving the Vij into the matrix
//---------------------
//[[Rcpp::export]]
void Save_BigMat_fun_Rcpp(SEXP &BigInfoMatDescInst, Rcpp::NumericVector const &GlobalNodesDist,
                          int &k, int &StartInd, int &EndInd,
                          Rcpp::NumericVector &InteractionPairs,
                          Rcpp::NumericVector const &NiNjIndeces,
                          Rcpp::NumericVector const &NiNjMat){
    // ---------------------------------------
    // Convert to Rcpp::XPtr<BigMatrix> objects:
    // ---------------------------------------
    Rcpp::XPtr<BigMatrix> xpBigInfoMatDescInst=Rcpp::XPtr<BigMatrix>(BigInfoMatDescInst);
    // ---------------------------------------
    // --------create MatrixAccessors
    // ---------------------------------------
    MatrixAccessor<double> maBigInfoMatDescInst(*xpBigInfoMatDescInst);
    // ---------------------------------------
    // initiate the h peak>k:
    // ---------------------------------------
    k-=1;//the peak in c++
    int h=k+1;//the next peak for the GlobalNodesDist
    // ---------------------------------------
    // --------loop
    // ---------------------------------------
    for(int VectID_kh=StartInd;VectID_kh<=EndInd;VectID_kh++){
        // take distance:
        double Dkh=GlobalNodesDist[h];
        // take correct indeces for NiNjMat:
        int NiNjMat_k=NiNjIndeces[k]-1;
        int NiNjMat_h=NiNjIndeces[h]-1;
        // Compute Vij weight:
        double Vij=NiNjMat[NiNjMat_k]*NiNjMat[NiNjMat_h]/Dkh;
        // check if valid:
        if(std::isinf(Vij)){//out, the Dkh is 0 and merged
            Vij=NA_REAL;
        }
        // save:
        maBigInfoMatDescInst[0][VectID_kh]=Vij;
        // count combination:
        if(!std::isnan(Vij)){
            // else not inf, and not na, so intra:
            InteractionPairs[0]+=1;
        }
        // go to next
        h+=1;
    }
}
// done
//---------------------
// Function to get the total PETs in each quantile and assign the QCells into the InfMat
//---------------------
//[[Rcpp::export]]
void Get_QCellPETCounts_fun_Rcpp(Rcpp::NumericVector const &BinsVij,
                                 int const & BinsVijSize,
                                 Rcpp::NumericVector const &ObsVij,
                                 Rcpp::NumericMatrix &InteractionInfMat,
                                 Rcpp::NumericVector const &AllInteIndeces,
                                 Rcpp::NumericVector &QCellPETCountsVij){
    //---------------------
    // Initiate the return:
    //---------------------
    int AllInteIndecesSize=AllInteIndeces.size();
    //---------------------
    // Loop on the leftovers
    //---------------------
    for(int index=0;index<AllInteIndecesSize;index++){
        // take interaction
        int CurInt=AllInteIndeces[index]-1;//Now it is in c++ index
        // take Vij:
        double Vij=ObsVij[index];
        // take nij:
        double nij=InteractionInfMat(CurInt,10);
        //---------------------
        // classify Vij
        //---------------------
        int QPosVij=BinsVijSize-1; //last position
        for(int j=0;j<BinsVijSize;j++){
            if(Vij<=BinsVij[j]){
                QPosVij=j;
                break;
            }
        }
        // save the QCell:
        InteractionInfMat(CurInt,11)=QPosVij;
        // count:
        QCellPETCountsVij[QPosVij]+=nij;
    }
}
// done
//---------------------
// Function to convert get the totals ni each QCell
//---------------------
//[[Rcpp::export]]
void Get_QCellCombCounts_fun_Rcpp(Rcpp::NumericVector const &BinsVij,
                                  int const &BinsVijSize,
                                  SEXP &BigInfoMatDescInst,
                                  Rcpp::NumericVector const &VkhOrder,
                                  Rcpp::NumericVector &QCellCombCountsVij_Net,
                                  int const &StartInd, int const &EndInd){
    // ---------------------------------------
    // Convert to Rcpp::XPtr<BigMatrix> objects:
    // ---------------------------------------
    Rcpp::XPtr<BigMatrix> xpBigInfoMatDescInst=Rcpp::XPtr<BigMatrix>(BigInfoMatDescInst);
    // ---------------------------------------
    // --------create MatrixAccessors
    // ---------------------------------------
    MatrixAccessor<double> maBigInfoMatDescInst(*xpBigInfoMatDescInst);
    //---------------------
    // initialize
    //---------------------
    int QPosVij=0;//the position on the QCellCombCountsVij_Net since VkhOrder
    // gives sorted vij
    int TotVkh=EndInd-StartInd+1;
    //---------------------
    // loop:
    //---------------------
    for(int i=0;i<TotVkh;i++){
        // take the order:
        int Order_kh=VkhOrder[i]-1;//c++ index
        // get the correct position on the Dij:
        double VkhPos=StartInd+Order_kh;//c++ index
        // get the Vij:
        double Vij=maBigInfoMatDescInst[0][int(VkhPos)];
        // check if na and break
        if(std::isnan(Vij)) break;//end of valid combinations
        //---------------------
        // classify Dij:
        //---------------------
        if(Vij>BinsVij[QPosVij]){
            // Find the next cell:
            int QPosLast=BinsVijSize-1;//ensures last cell
            for(int j=QPosVij+1;j<BinsVijSize;j++){
                if(Vij<=BinsVij[j]){
                    QPosLast=j;
                    break;
                }
            }
            QPosVij=QPosLast;
        }
        // count:
        QCellCombCountsVij_Net[QPosVij]+=1;
    }
}
// done
//####################################################
//Functions for running the interaction analysis:
//####################################################
//---------------------
//Main function for assessing the current interaction:
//---------------------
// [[Rcpp::export]]
SEXP Assess_Interaction_fun_Rcpp(int &CurInt,//the current interaction id in R index
                                 Rcpp::NumericMatrix &InteractionInfMat,//the InteractionInfMat
                                 Rcpp::Function &Poiss_fun,
                                 Rcpp::NumericMatrix const &BinMatVij){
    // ---------------------------------------
    // Initiate the results in a vector
    // ---------------------------------------
    Rcpp::NumericVector pValues_round_i(2);
    // ---------------------------------------
    // Initiate the return:
    // A vector with the first element the CurInt corresponding in R index
    // and second the p-value for computing the FDR
    // ---------------------------------------
    pValues_round_i[0]=CurInt;//R index
    CurInt-=1;//C++ index
    // ---------------------------------------
    // Gather information about the current interaction:
    // ---------------------------------------
    double nij=InteractionInfMat(CurInt,10);
    // get QCEll:
    int QCell_Vij=InteractionInfMat(CurInt,11);//already in c++
    // get expected:
    double lij=BinMatVij(QCell_Vij,2);
    // ---------------------------------------
    // Get the p-value of the interaction
    // ---------------------------------------
    double pval=Rcpp::as<double>(Poiss_fun(nij,lij));
    // ---------------------------------------
    // Save
    // ---------------------------------------
    pValues_round_i[1]=pval;
    // ---------------------------------------
    // return:
    // ---------------------------------------
    return Rcpp::wrap(pValues_round_i);
}
// done
//####################################################
//Functions for updating the significant interactions
//####################################################
//---------------------
// Function for updating the indeces of the rest of the interactions to be added
// This is usefull for the network update, if an interaction is already updated
// then the network will not be updated again
//---------------------
//[[Rcpp::export]]
void Update_ToBeAddedInter_fun_Rcpp(Rcpp::NumericMatrix &InteractionInfMat,
                                    int &k, int &h, int &i, int &TR_Si,
                                    Rcpp::NumericVector &LastInteractions,
                                    int &Chrom12ID_i){
    //---------------------
    // NOTE: the interaction updated in R is in i index
    // In c++ the i index is the next interaction, so start there.
    //---------------------
    for(int j=i;j<TR_Si;j++){
        // take ID:
        int ID_j=LastInteractions[j]-1;//because R to c++ index
        // change the Node_i/j if needed (note the chromosome ID has to also be the same):
        if((InteractionInfMat(ID_j,2)==h)&&(InteractionInfMat(ID_j,13)==Chrom12ID_i)) InteractionInfMat(ID_j,2)=k;
        if((InteractionInfMat(ID_j,3)==h)&&(InteractionInfMat(ID_j,13)==Chrom12ID_i)) InteractionInfMat(ID_j,3)=k;
    }
}
// done
//---------------------
// Function for updating the rest of the interactions Nodes and checking
// for bi-products:
//---------------------
//[[Rcpp::export]]
SEXP Check_BiProd_fun_Rcpp(Rcpp::NumericMatrix &InteractionInfMat,
                           int &k, int &h, Rcpp::NumericVector &AllInteIndeces, double &TotBiRem, int &Chrom12ID_i,
                           int &OrdersCount){
    // ---------------------------------------
    // Initiate:
    // ---------------------------------------
    int AllInteIndecesSize=AllInteIndeces.size();
    Rcpp::NumericVector BiProductIDSreject;//R ids to be added and removed from AllInteIndeces
    Rcpp::NumericVector BiProductIDSaccepted;//accepted bi-products
    int TotBiAcc=0;
    //---------------------
    // loop and update each element in the InteractionInfMat
    //---------------------
    for(int j=0;j<AllInteIndecesSize;j++){
        // take ID:
        int ID_j=AllInteIndeces[j]-1;//because R to c++ index
        //---------------------
        // change the Node_i/j if needed:
        //---------------------
        if((InteractionInfMat(ID_j,2)==h)&&(InteractionInfMat(ID_j,13)==Chrom12ID_i)) InteractionInfMat(ID_j,2)=k; //change name
        if((InteractionInfMat(ID_j,3)==h)&&(InteractionInfMat(ID_j,13)==Chrom12ID_i)) InteractionInfMat(ID_j,3)=k;//change name
        // ---------------------------------------
        // check if bi-product
        // ---------------------------------------
        if(InteractionInfMat(ID_j,2)==InteractionInfMat(ID_j,3)){// then this is a biproduct.
            if(InteractionInfMat(ID_j,9)<0.05){
                // accept it
                // order:
                InteractionInfMat(ID_j,12)=OrdersCount;
                // save the biproducts:
                BiProductIDSaccepted.push_back(AllInteIndeces[j]);//R indeces
                TotBiAcc+=1;
            }else{
                // reject it
                // count:
                TotBiRem+=1;
                // order:
                InteractionInfMat(ID_j,12)=NA_REAL;//indicating a bi-product non significant
                // save the biproducts:
                BiProductIDSreject.push_back(AllInteIndeces[j]);//R indeces
            }
        }
    }
    Rcpp::List BiPorductsInfo=Rcpp::List::create(Rcpp::Named("BiProductIDSreject")=BiProductIDSreject,
                                                 Rcpp::Named("TotBiRem")=TotBiRem,
                                                 Rcpp::Named("BiProductIDSaccepted")=BiProductIDSaccepted,
                                                 Rcpp::Named("TotBiAcc")=TotBiAcc);

    return Rcpp::wrap(BiPorductsInfo);
}
// done
//####################################################
//Function For summarizing the interaction information
//####################################################
//[[Rcpp::export]]
SEXP Get_InteractionInfo_fun_Rcpp(Rcpp::NumericMatrix &InteractionInfMat,// the InteractionInfMat
                                  int &NInteractions){//total interactions
    // ---------------------------------------
    // --------Initiate the return:
    // ---------------------------------------
    Rcpp::NumericMatrix InteractionInfo(NInteractions,6);
    // ---------------------------------------
    // --------Loop and save:
    // ---------------------------------------
    for(int i=0;i<NInteractions;i++){
        // PBS from:
        InteractionInfo(i,0)=InteractionInfMat(i,0);
        // PBS to:
        InteractionInfo(i,1)=InteractionInfMat(i,1);
        // p-value:
        InteractionInfo(i,2)=InteractionInfMat(i,8);
        // FDR:
        InteractionInfo(i,3)=InteractionInfMat(i,9);
        // order:
        InteractionInfo(i,4)=InteractionInfMat(i,12);
        //total interaction PETs
        InteractionInfo(i,5)=InteractionInfMat(i,10);
    }
    return Rcpp::wrap(InteractionInfo);
}
// done
//####################################################
//Function used in GetSignInteractions.GenomeMap for subseting significant interactions
//####################################################
//[[Rcpp::export]]
SEXP SubsetSignificantInteractions_fun_Rcpp(int const &NInteractionInfo,
                                            Rcpp::NumericVector const &FDR,
                                            Rcpp::NumericVector const &Order,
                                            double const &threshold){
    int MaxOrder=0;//the maximum order to take, has to start at zero
    // loop and check. Since in same order is the same FDR, break if
    // the FDR gets bigger than the threshold, or iff the orders skips numbers
    // This might happen in case one whole order is out from before, since
    // you have subset the unsignifikant.
    for(int i=0;i<NInteractionInfo;i++){
        if((FDR[i]<threshold)&&((Order[i]==MaxOrder+1)||(Order[i]==MaxOrder))){
            // Then it passes, but only if the order order is the previous plus one.
            MaxOrder = Order[i];
        }else{
            break;
        }
    }
    return Rcpp::wrap(MaxOrder);
}
