#include <Rcpp.h>
using namespace Rcpp;

//####################################################
//Functions for preparing local fit in a region:
//####################################################
// ---------------------
// Global function for getting the region info.
// ---------------------
Rcpp::List PrepLocFit_Rcpp(Rcpp::DataFrame &bppass_x,
                           Rcpp::DataFrame const &ChromInf){
    // Inputs:
    // bppass_x=the data frame of the data with Utag,DTag,Region,Chrom,MainIndex info
    // ChromInf=data frame with chrom, size and tot PETs in each chromosome.
    // Outputs(a list):
    // DFit= matrix for the UTag and DTag
    // MainIndex=vector for the main indeces for each PET
    // Chrom=character for the chromosome
    // Regio=integer for the region
    // N=#PETs
    // ChromSize=integer the size of the chromosome(used in fitting to cut out boundaries)
    // NoisePDF=vector of the noise PDF
    // InParamlist=list of the initial parameters.

    // ------------------------
    // Declare functions:
    // ------------------------
    //For Spp and some other matrices/constants
    void Get_Dk2_Rcpp(int const &N,Rcpp::NumericMatrix const &DFit,
                      Rcpp::NumericVector &Dk2,double &Dk2sum0,
                      double &Dk2sum1, int const &Nhalf,double &TagMax,
                      double &TagMin);
    // For EM Spp:
    void SppEM_Rcpp(int const &N,
                    double const &pi,double const &Dk2sum0,
                    double const &Dk2sum1,int const &Nhalf,
                    Rcpp::NumericVector const &Dk2,
                    Rcpp::NumericVector &FeatureID,int &Nsub);
    // declare function for kernels and local maxima:
    Rcpp::NumericMatrix MatchKernelsPairs_Rcpp(Rcpp::NumericVector const &KernSeq,
                                               Rcpp::NumericMatrix const &DFit,
                                               Rcpp::NumericVector const &FeatureID,
                                               int const &Nsub,
                                               int const &KernSeqLength,
                                               double const &BW,
                                               double const &pi);

    //declare function for initiating:
    Rcpp::List InitializeKernels_Rcpp(Rcpp::NumericMatrix const &KernelPairs,
                                      Rcpp::NumericMatrix const &DFit,int const &N);
    // ------------------------
    // take general region info:
    // ------------------------
    Rcpp::StringVector ChromVect=bppass_x["Chrom"];
    std::string Chrom=Rcpp::as<std::string>(ChromVect[0]);//chromosome of region
    Rcpp::NumericVector RegionVect=bppass_x["Region"];
    int Region=RegionVect[0];//region ID.
    // Rcout<<"Which: "<<Chrom<<"--"<<Region<<std::endl;
    Rcpp::NumericVector MainIndex=bppass_x["MainIndex"];//the MainIndex vector.
    int N=MainIndex.size();//number of observations, #PETs.
    int ChromSize;//initiate to find them.
    // break the ChromInf data frame:
    Rcpp::StringVector ChromInfChrom=ChromInf["Chrom"];
    Rcpp::NumericVector ChromInfsize=ChromInf["size"];
    for(int chr=0;chr<ChromInfChrom.size();chr++){
        if(ChromInfChrom[chr]==Chrom){
            ChromSize=ChromInfsize[chr];
            break;
        }
    }//will find one occurence.
    // create DFit:
    Rcpp::NumericMatrix DFit(N,2);
    DFit(_,0)=Rcpp::as<Rcpp::NumericVector>(bppass_x["UTag"]);
    DFit(_,1)=Rcpp::as<Rcpp::NumericVector>(bppass_x["DTag"]);
    // ------------------------
    // Evaluate Noise PDF. Also in the loop get: Dk2-second nearest distance
    // Different sums for the lambdas for Spp.
    // Max and Min Tag coordinates for generating sequence for the Kernels
    // -----------------------
    // initiate:
    // Rcpp::NumericVector NoisePDF(N);//Noise PDF
    Rcpp::NumericVector Dk2(N);//second nearest distance PET.only dk^2 is needed for SPP.
    double Dk2sum0=0,Dk2sum1=0;//for lambda initials for SPP
    const int Nhalf=std::ceil(N/2.0);//half observations needed in general in some functions
    double TagMax=-INFINITY,TagMin=INFINITY;//for region sizes
    // call void to find variables:
    Get_Dk2_Rcpp(N,DFit,Dk2,Dk2sum0,Dk2sum1,Nhalf,TagMax,TagMin);
    // find noise pdf:
    double NoisePDF=0.4*std::pow(TagMax-TagMin,-2.0);
    // ------------------------
    // Run Spatial poisson process for Noise/Feature estimation.
    // -----------------------
    // initiate:
    const double pi=3.14159265358979323846;
    Rcpp::NumericVector FeatureID;//Vector with the indeces of the features.
    int Nsub;//total number of features.
    // call function(get the feature subset matrix):
    SppEM_Rcpp(N,pi,Dk2sum0,Dk2sum1,Nhalf,Dk2,FeatureID,Nsub);
    // ------------------------
    // Run Kernel Estimation:
    // -----------------------
    // create sequence to evaluate the kernels
    double SeqStep=std::max(std::round((TagMax-TagMin)/std::min(double(2*N),500.0)),1.0);//the step of the sequence, minimum 1
    Rcpp::NumericVector KernSeq(1);//the sequence to be used in Kernels
    // fill values by expanding the sequence:
    KernSeq[0]=TagMin;
    int KernSeqLength=KernSeq.size();//size of sequence to generate
    while(KernSeq[KernSeqLength-1]<TagMax){
        KernSeq.push_back(KernSeq[KernSeqLength-1]+SeqStep);//append element at the end
        KernSeqLength=KernSeq.size();//update sequence length
    }
    // -----------------------
    // call function for kernel local maxima:
    // -----------------------
    Rcpp::NumericMatrix KernelPairs=MatchKernelsPairs_Rcpp(KernSeq,DFit,
                                                           FeatureID,Nsub,
                                                           KernSeqLength,
                                                           50.0,pi);
    // -----------------------
    // find initials:
    // -----------------------
    Rcpp::List InParam=InitializeKernels_Rcpp(KernelPairs,DFit,N);
    // ------------------------
    //Create output:
    // -----------------------
    Rcpp::List Res=Rcpp::List::create(Rcpp::Named("DFit")=DFit,
                                      Rcpp::Named("MainIndex")=MainIndex,
                                      Rcpp::Named("Chrom")=Chrom,
                                      Rcpp::Named("Region")=Region,
                                      Rcpp::Named("N")=N,
                                      Rcpp::Named("ChromSize")=ChromSize,
                                      Rcpp::Named("NoisePDF")=NoisePDF,
                                      Rcpp::Named("InParam")=InParam,
                                      Rcpp::Named("KernSeq")=KernSeq);
    return Res;
}
//Done

// ---------------------
// Function for Evaluating Noise PDF, it also returns other things used in subsequent functions:
// ---------------------
void Get_Dk2_Rcpp(int const &N,Rcpp::NumericMatrix const &DFit,
                  Rcpp::NumericVector &Dk2,double &Dk2sum0,double &Dk2sum1,
                  int const &Nhalf,double &TagMax,double &TagMin){
    // Initiate:
    Rcpp::NumericVector KNear(2);//save the nearest distances
    // loop to fill in
    for(int i=0;i<N;i++){//loop through tag i
        // -----------------------
        // start part of dk2/MaxMin Region:
        // -----------------------
        if(DFit(i,0)>TagMax) TagMax=DFit(i,0);//replace max UTag
        if(DFit(i,0)<TagMin) TagMin=DFit(i,0);//replace min UTag
        if(DFit(i,1)>TagMax) TagMax=DFit(i,1);//replace max DTag
        if(DFit(i,1)<TagMin) TagMin=DFit(i,1);//replace min DTag
        //for dk2, KNear[0]=first nearest, KNear[1]=second(keep that)
        KNear[0]=INFINITY,KNear[1]=INFINITY;
        // -----------------------
        // end part of dk2/MaxMin Region:
        // -----------------------
        for(int j=0;j<N;j++){//loop through tag j(internal)
            // -----------------------
            // start part of dk2/MaxMin Region:
            // -----------------------
            // find euclidean distance for jl:
            // NOTE: i need dkij^2 to save, so no need to take sqrt and pow afterwords.
            // just find dkij^2 imediatelly.
            if(i!=j){//else it is self, so skip it
                double dkij2=std::pow(DFit(i,0)-DFit(j,0),2.0)+std::pow(DFit(i,1)-DFit(j,1),2.0);
                // save first and second nearest
                if(dkij2<KNear[0]){
                    KNear[1]=KNear[0];//move the previous first max to second.
                    KNear[0]=dkij2;//replace first max.
                }else if(dkij2<KNear[1]){
                    KNear[1]=dkij2;//to return
                }
            }
            // -----------------------
            // end part of dk2/MaxMin Region:
            // -----------------------
        }
        // -----------------------
        // start part of dk2/MaxMin Region:
        // -----------------------
        Dk2[i]=KNear[1];//second nearest distance
        //save for lambdas
        if(i<=Nhalf){
            Dk2sum0+=Dk2[i];
        }else{
            Dk2sum1+=Dk2[i];
        }
        // -----------------------
        // end part of dk2/MaxMin Region:
        // -----------------------
    }
}
//Done.

// ---------------------
// Function Spatial Poisson Process-returns the matrix used in Kernels:
// ---------------------
void SppEM_Rcpp(int const &N,
                double const &pi,double const &Dk2sum0,
                double const &Dk2sum1,int const &Nhalf,
                Rcpp::NumericVector const &Dk2,
                Rcpp::NumericVector &FeatureID,int &Nsub){
    // ------------
    // initiate parameters and critical values:
    // ------------
    double CritVal=INFINITY,p_0=0.5,p_1=0.5,lambda_0=2*Nhalf/(pi*Dk2sum0),
        lambda_1=2*Nhalf/(pi*Dk2sum1);
    int Emit=1;//em iterations.
    Rcpp::NumericVector Feature(N);//save feature classes here(0=noise/1=feature)
    // declare function
    void EMstepSpp_Rcpp(double &p_0, double &p_1,double &lambda_0,
                        double &lambda_1,Rcpp::NumericVector const &Dk2,
                        double const &pi, int const &N,
                        Rcpp::NumericVector &Feature,double &CritVal,int &Nsub);
    //-----------
    //Run EM:
    //-----------
    while((Emit<=500)&&(CritVal>1e-6)){
        // EM step:
        EMstepSpp_Rcpp(p_0,p_1,lambda_0,lambda_1,Dk2,pi,N,Feature,CritVal,Nsub);
        Emit+=1;
    }
    //-----------
    //get the FeatureID vector for output:
    //-----------
    if(Nsub<=2){//at least 2 observations.
        Rcpp::NumericVector FeaturesIDout(N);
        for(int i=0;i<N;i++){
            FeaturesIDout[i]=i;
        }
        // replace to return:
        FeatureID=FeaturesIDout;
        Nsub=N;//since the Nsub is the whole N now.
    }else{
        Rcpp::NumericVector FeaturesIDout(Nsub);
        int IDind=0;//for index in FeatureIDout
        for(int i=0;i<N;i++){
            if(Feature[i]==1){
                FeaturesIDout[IDind]=i;
                IDind+=1;
            }
        }
        // replace to return:
        FeatureID=FeaturesIDout;
    }
}
//Done.

//---------------------
//Estep And Mstep Function for Spp:
//---------------------
void EMstepSpp_Rcpp(double &p_0, double &p_1,double &lambda_0,double &lambda_1,
                    Rcpp::NumericVector const &Dk2,double const &pi,
                    int const &N,Rcpp::NumericVector &Feature,
                    double &CritVal,int &Nsub){
    //-----------
    // initate for paramters: save also feature at once and return it
    //-----------
    Rcpp::NumericVector n_g(2), dk2resp(2),Upp_g(2),Uplambda_g(2);
    Nsub=0;//have to set it to zero for each run.
    // loop from 1:N and fill up rows:
    for(int i=0;i<N;i++){
        // find first resp:
        double fij01=std::exp((lambda_0-lambda_1)*pi*Dk2[i]+
                              2.0*std::log(lambda_1/lambda_0)+std::log(p_1/p_0));
        double Respi0=1.0/(1.0+fij01);
        // Note Respi1=1-Respi0, so if Respi0>Respi1=1-Respi0=>
        // Respi0>1/2, and lambda_0>lambda_1, then feature in class 0.
        // The same if Respi0<1/2 and lambda_1>lambda_0 feature in class 1.
        // save feature (on previous param but at convergence is ok)
        if(((lambda_0>lambda_1)&&(Respi0>(1.0-Respi0)))||
           ((lambda_1>lambda_0)&&((1-Respi0)>Respi0))){
            Feature[i]=1;
            Nsub+=1;
        }else{
            Feature[i]=0;//this has to be here because Feture is passed by reference!!!
            // so if i have a value of one before, i need to set it to zero now
        }
        // the parameters:
        n_g[0]+=Respi0;
        n_g[1]+=1.0-Respi0;
        dk2resp[0]+=Dk2[i]*Respi0;//updating lambda0
        dk2resp[1]+=Dk2[i]*(1.0-Respi0);//for updating lambda1
    }
    //-----------
    // fill parameters here:
    //-----------
    Upp_g[0]=n_g[0]/N;// update p_0
    Upp_g[1]=n_g[1]/N;// update p_1
    Uplambda_g[0]=2.0*n_g[0]/(pi*dk2resp[0]);//update lambda0
    Uplambda_g[1]=2.0*n_g[1]/(pi*dk2resp[1]);//update lambda1
    //-----------
    // Update criterio:
    //-----------
    CritVal=std::sqrt(std::pow(Upp_g[0]-p_0,2.0)+std::pow(Upp_g[1]-p_1,2.0)+
        std::pow(Uplambda_g[0]-lambda_0,2.0)+std::pow(Uplambda_g[1]-lambda_1,2.0));
    //-----------
    //create output:
    //-----------
    p_0=Upp_g[0],p_1=Upp_g[1],lambda_0=Uplambda_g[0],lambda_1=Uplambda_g[1];

}
//Done.

// ---------------------
// Function for finding peaks based on Kernel estimation and a given BW.
// ---------------------
Rcpp::NumericMatrix MatchKernelsPairs_Rcpp(Rcpp::NumericVector const &KernSeq,
                                           Rcpp::NumericMatrix const &DFit,
                                           Rcpp::NumericVector const &FeatureID,
                                           int const &Nsub,
                                           int const &KernSeqLength,
                                           double const &BW,
                                           double const &pi){
    //-----------
    // initate
    //-----------
    // DivKern is the Kernel normalizations:
    double DivKern=1.0/(std::sqrt(2.0*pi)*Nsub*BW);
    // Matrix for saving Kernel information:
    Rcpp::NumericMatrix KernelInfo(KernSeqLength-1,8);
    // 0=Upoint,1=Dpoint,2-3=U/Dmax=(1 if local max, 0 if not)
    // 4=UDeps(1 if both dens bigger than Eps, 0 if not)(needed if not pair found)
    // 5-6=U/Ddens densities of the U and D tags, 7 maximum of the densities of the row
    int UDeps=0,UMinInd=0,DMinInd=0,ValidPairs=0;//the total U and D which are higher than eps,
    // the minimum index for the u/d local max, and the total valid pairs.
    Rcpp::NumericVector Uvalid,Dvalid;//the indeces of the valid pairs
    //-----------
    // Run to find Kernels:
    //-----------
    for(int i=0;i<(KernSeqLength-1);i++){
        // take UTag,DTag for i:
        double UKernSum_i=0,DKernSum_i1=0;//sums of current kernel dens.
        // loop thought the rest, use Nsub for the subset of the feature:
        for(int j=0;j<Nsub;j++){
            // take UTag,DTag for j:
            UKernSum_i+=std::exp(-std::pow((KernSeq[i]-DFit(FeatureID[j],0))/BW,2.0)/2.0);//the i kernel seq
            DKernSum_i1+=std::exp(-std::pow((KernSeq[i+1]-DFit(FeatureID[j],1))/BW,2.0)/2.0);//the next i+1 kernel seq
        }
        // Note: The KernelInfo(_,0) is filled from 1 to KernSeqLength-1
        // Note: The KernelInfo(_,1) is filled from 2 to KernSeqLength.
        KernelInfo(i,0)=KernSeq[i];//Utag point
        KernelInfo(i,1)=KernSeq[i+1];//Dtag point
        KernelInfo(i,5)=DivKern*UKernSum_i;//Utag dens
        KernelInfo(i,6)=DivKern*DKernSum_i1;//Dtag dens
        KernelInfo(i,7)=std::max(KernelInfo(i,5),KernelInfo(i,6));//max density
        if(KernelInfo(i,5)>DBL_EPSILON&&KernelInfo(i,6)>DBL_EPSILON) KernelInfo(i,4)=1;
        if(KernelInfo(i,4)==1) UDeps+=1;//count up.
        //-----------
        //update local maxima:
        //-----------
        if(i==1){
            // Utag check: at i==0
            if((KernelInfo(i-1,5)>KernelInfo(i,5)+DBL_EPSILON)&&//bigger than next
               (KernelInfo(i-1,5)>DBL_EPSILON)){//bigger than eps
                // then local maxima on UTag:
                KernelInfo(i-1,2)=1;
            }
            // Dtag check: at i==0
            if((KernelInfo(i-1,6)>KernelInfo(i,6)+DBL_EPSILON)&&//bigger than next
               (KernelInfo(i-1,6)>DBL_EPSILON)){//bigger than eps
                // then local maxima on UTag:
                KernelInfo(i-1,3)=1;
            }
        }else if(i==(KernSeqLength-2)){
            // Utag check: at i
            if((KernelInfo(i,5)>KernelInfo(i-1,5)+DBL_EPSILON)&&//bigger than previous
               (KernelInfo(i,5)>DBL_EPSILON)){//bigger than eps
                // then local maxima on UTag:
                KernelInfo(i,2)=1;
            }
            // Dtag check: at i
            if((KernelInfo(i,6)>KernelInfo(i-1,6)+DBL_EPSILON)&&//bigger than previous
               (KernelInfo(i,6)>DBL_EPSILON)){//bigger than eps
                // then local maxima on UTag:
                KernelInfo(i,3)=1;
            }
        }else{
            // for i==2 to i==(KernSeqLength-3), check i-1 for both tags
            // Utag check:
            if((KernelInfo(i-1,5)>KernelInfo(i,5)+DBL_EPSILON)&&//bigger than next
               (KernelInfo(i-1,5)>KernelInfo(i-2,5)+DBL_EPSILON)&&//bigger than the previous
               (KernelInfo(i-1,5)>DBL_EPSILON)){//bigger than eps
                // then local maxima on UTag:
                KernelInfo(i-1,2)=1;
            }
            // Dtag check:
            if((KernelInfo(i-1,6)>KernelInfo(i,6)+DBL_EPSILON)&&//bigger than next
               (KernelInfo(i-1,6)>KernelInfo(i-2,6)+DBL_EPSILON)&&//bigger than the previous
               (KernelInfo(i-1,6)>DBL_EPSILON)){//bigger than eps
                // then local maxima on UTag:
                KernelInfo(i,3)=1;
            }
        }
        //-----------
        //check for valid pairs so far:
        //-----------
        for(int d=i;d>=DMinInd;d--){//for the Dtag column
            if(KernelInfo(d,3)==1){//local Dtag max
                for(int u=d;u>=UMinInd;u--){//for the Utag column
                    if(KernelInfo(u,2)==1){
                        // then a pair.
                        ValidPairs+=1;//count the pair
                        Uvalid.push_back(u);//save Utag
                        Dvalid.push_back(d);//save Utag
                        UMinInd=u+1;//increase UMinInd to the above UTag
                        DMinInd=d+1;//increase UMinInd to the above DTag
                        break;//break loop
                    }
                }
            }
        }
    }
    //-----------
    // Prepare the return:
    //-----------
    Rcpp::NumericMatrix Res;//initiate output.
    if(ValidPairs!=0){
        // i have pairs
        Rcpp::NumericMatrix KernelMatrix(ValidPairs,3);
        for(int i=0;i<ValidPairs;i++){
            KernelMatrix(i,0)=KernelInfo(Uvalid[i],0);//Utag
            KernelMatrix(i,1)=KernelInfo(Dvalid[i],1);//Dtag
            KernelMatrix(i,2)=std::max(KernelInfo(Uvalid[i],5),KernelInfo(Dvalid[i],6));//max dens
        }
        Res=KernelMatrix;//replace to return
    }else{
        // return one pair:
        Rcpp::NumericMatrix KernelMatrix(1,3);
        KernelMatrix(0,0)=mean(DFit(_,0));//Utag
        KernelMatrix(0,1)=mean(DFit(_,1));//Dtag
        KernelMatrix(0,2)=1;//max dens
        Res=KernelMatrix;//replace to return
    }
    return Res;
}
//Done

// ---------------------
// Function initiating the kernels:
// ---------------------
Rcpp::List InitializeKernels_Rcpp(Rcpp::NumericMatrix const &KernelPairs,
                                  Rcpp::NumericMatrix const &DFit,int const &N){
    // -------------
    // initiate:
    // -------------
    int G=KernelPairs.nrow();
    double pgval_noise=1.0/(double(G)+1.0);//adding noise
    // create the p_g:
    Rcpp::NumericVector p_g_noise(G+1,pgval_noise),
    sdx_g(G),lambdax_g(G),mx_g(G),sdy_g(G),lambday_g(G),my_g(G);
    // -------------
    // loop though clusters and fill the parameters
    // -------------
    for(int g=0;g<G;g++){
        // get peak points:
        double UPointStart_g=KernelPairs(g,0)-250;
        double UPointEnd_g=KernelPairs(g,0)+250;
        double DPointStart_g=KernelPairs(g,1)-250;
        double DPointEnd_g=KernelPairs(g,1)+250;
        // find which tags are in the interval:
        // for var use: E(x^2)-(E(x))^2
        int Ng=0;//the total tags included in the above intervals
        double UtagSum=0,DtagSum=0,UtagSumP2=0,DtagSumP2=0;//sum for means and variances
        for(int i=0;i<N;i++){
            // check if either of the tags are included in the intervals
            bool i_inn=(UPointStart_g<=DFit(i,0)&&DFit(i,0)<=UPointEnd_g)||
                (DPointStart_g<=DFit(i,1)&&DFit(i,1)<=DPointEnd_g);
            if(i_inn){
                // save the info
                Ng+=1;
                UtagSum+=DFit(i,0);
                DtagSum+=DFit(i,1);
                UtagSumP2+=std::pow(DFit(i,0),2.0);
                DtagSumP2+=std::pow(DFit(i,1),2.0);
            }
        }
        // -------------
        // define modes:
        // -------------
        mx_g[g]=KernelPairs(g,0);
        my_g[g]=KernelPairs(g,1);
        // -------------
        // define E(X) and E(x^2):
        // -------------
        double Meanx_g=UtagSum/Ng;
        double Meany_g=DtagSum/Ng;
        double Meanx2_g=UtagSumP2/Ng;
        double Meany2_g=DtagSumP2/Ng;
        // -------------
        // define sd:
        // -------------
        sdx_g[g]=std::sqrt(Meanx2_g-std::pow(Meanx_g,2.0));
        sdy_g[g]=std::sqrt(Meany2_g-std::pow(Meany_g,2.0));
        // -------------
        // define lambdas/skewness:
        // -------------
        lambdax_g[g]=-0.5;
        lambday_g[g]=0.5;
    }
    // -------------
    // create output:
    // -------------
    Rcpp::List Res=Rcpp::List::create(Rcpp::Named("p_g")=p_g_noise,
                                      Rcpp::Named("sdx_g")=sdx_g,
                                      Rcpp::Named("lambdax_g")=lambdax_g,
                                      Rcpp::Named("mx_g")=mx_g,
                                      Rcpp::Named("sdy_g")=sdy_g,
                                      Rcpp::Named("lambday_g")=lambday_g,
                                      Rcpp::Named("my_g")=my_g);
    return Res;
}
//Done

//####################################################
//--------------------- Local Model Fit functions:
//####################################################
//---------------------
//Function for preparing and fitting in each region: (exported in R)
//---------------------
// [[Rcpp::export]]
SEXP FitCallLocal_fun_Rcpp( Rcpp::DataFrame &bppass_x,
                            Rcpp::DataFrame &ChromInf){
    // ---------------------------------------
    // Function declaration:
    // ---------------------------------------
    // function for preparing the local fit
    Rcpp::List PrepLocFit_Rcpp(Rcpp::DataFrame &bppass_x,Rcpp::DataFrame const &ChromInf);
    // function for fitting the local fit.
    Rcpp::List GetLocFit_Rcpp(Rcpp::NumericMatrix const &DFit,Rcpp::List const &InParam,
                              int const &N,int const &Region,
                              int const &ChromSize, double const &NoisePDF,
                              Rcpp::NumericVector const &MainIndex,std::string &Chrom,
                              Rcpp::NumericVector const &KernSeq);
    // ---------------------------------------
    // --------Prepare the data for fitting:
    // ---------------------------------------
    Rcpp::List LocInf=PrepLocFit_Rcpp(bppass_x,ChromInf);
    // break:
    Rcpp::List InParam=LocInf["InParam"];
    Rcpp::NumericVector MainIndex=LocInf["MainIndex"];
    std::string Chrom=LocInf["Chrom"];
    int Region=LocInf["Region"];
    int N=LocInf["N"];
    int ChromSize=LocInf["ChromSize"];
    double NoisePDF=LocInf["NoisePDF"];//noise pdf
    Rcpp::NumericMatrix DFit=LocInf["DFit"];//Fitting data
    Rcpp::NumericVector KernSeq=LocInf["KernSeq"];//Kernel Sequence
    // ---------------------------------------
    // --------Fit the peaks and find optimal model:
    // ---------------------------------------
    Rcpp::List FitReslocal=GetLocFit_Rcpp(DFit,InParam,N,Region,ChromSize,
                                          NoisePDF,MainIndex,Chrom,KernSeq);
    return Rcpp::wrap(FitReslocal);
}

//Done

//---------------------
//Function for fitting the local models and chosing the best one by BIC:
//---------------------
Rcpp::List GetLocFit_Rcpp(Rcpp::NumericMatrix const &DFit,Rcpp::List const &InParam,
                          int const &N,int const &Region,
                          int const &ChromSize, double const &NoisePDF,
                          Rcpp::NumericVector const &MainIndex,std::string &Chrom,
                          Rcpp::NumericVector const &KernSeq){
    //------------
    // declare functions:
    //------------
    // Function for running EM for Laplace Model
    Rcpp::List MainEMLoop_Rcpp(Rcpp::NumericMatrix const &DFit,
                               Rcpp::List const &InParam,int const &N,
                               double const &NoisePDF);
    // declare function for merging:
    void MergeOvPeak_Rcpp(Rcpp::NumericMatrix const &DFit,
                          int const &N,double const &NoisePDF,
                          Rcpp::List &OptParam,Rcpp::NumericVector &OptClass,
                          Rcpp::NumericVector &OptClassTot,
                          Rcpp::NumericVector const &KernSeq);
    // declare function for finding peak info:
    Rcpp::DataFrame GetPeakInf_Rcpp(Rcpp::List const &OptParam,
                                    Rcpp::NumericVector const &OptClassTot,
                                    int const &ChromSize, int const &Region,
                                    std::string const &Chrom,Rcpp::NumericVector &OldPeakNames);
    // declare function for fixing the classification at the end:
    void OptClass_Update_Rcpp(Rcpp::NumericVector &OptClass,Rcpp::NumericVector const &OldPeakNames,
                              int const &N);
    //------------
    // Initiate:
    //------------
    Rcpp::List OptParam;
    Rcpp::NumericVector OptClass,OptClassTot;
    bool ModelPass=FALSE;//to check if model is accepted
    //------------
    // Fit model:
    //------------
    Rcpp::List FitRes=MainEMLoop_Rcpp(DFit,InParam,N,NoisePDF);
    //------------
    // save only if no-noise model:
    //------------
    if(!Rcpp::as<bool>(FitRes["OnlyNoise"])){//if not noise:
        ModelPass=TRUE;
        OptParam=FitRes["Param"];
        OptClass=FitRes["Classification"];
        OptClassTot=FitRes["ClassTot"];
    }
    //------------
    // Create output:
    //------------
    Rcpp::List Res=Rcpp::List::create(Rcpp::Named("Classification")=R_NilValue,
                                      Rcpp::Named("PeakInf")=R_NilValue);
    //------------
    // take the best model: check if no model is fitted
    //------------
    if(ModelPass){// then the model is fitted else not
        //------------
        // Merge overlapping:
        //------------
        MergeOvPeak_Rcpp(DFit,N,NoisePDF,OptParam,OptClass,OptClassTot,KernSeq);
        //------------
        // Find Peak Info:
        //------------
        // get model info output:
        Rcpp::NumericVector OldPeakNames;//vector with the old peak names to fix the classes
        Rcpp::DataFrame PeakInf=GetPeakInf_Rcpp(OptParam,OptClassTot,
                                                ChromSize,Region,Chrom,OldPeakNames);
        //------------
        // Update output:
        //------------
        if(PeakInf.nrow()!=0){
            // fix the classification:
            OptClass_Update_Rcpp(OptClass,OldPeakNames,N);
            //create the matrix of the Classification:
            Rcpp::NumericMatrix Classification(N,2);
            // fill in:
            Classification(_,0)=OptClass;//classes
            Classification(_,1)=MainIndex;//index for complete data
            // Create list:
            Res["Classification"]=Classification;
            Res["PeakInf"]=PeakInf;
        }

    }

    return Res;

}
//Done


//---------------------
//Function for updating the optimal class in case any peak is removed:
//---------------------
void OptClass_Update_Rcpp(Rcpp::NumericVector &OptClass,Rcpp::NumericVector const &OldPeakNames,int const &N){
    // First take the Peaks:
    int G=OldPeakNames.size();
    // Then loop and check if any OptClass has to be set to zero:
    for(int i=0;i<N;i++){
        bool Keep=FALSE;//if to keep or sat to zero
        for(int g=0;g<G;g++){
            if(OptClass[i]==OldPeakNames[g]){
                Keep=TRUE;//then it exists
                break;
            }
        }
        if(!Keep) OptClass[i]=0;//if not in then set to zero.
    }
}

// done

//####################################################
//--------------------- EM algorithm and NR algorithm:
//####################################################
//---------------------
//Main SLaplace EM algorithm Function:
//---------------------
Rcpp::List MainEMLoop_Rcpp(Rcpp::NumericMatrix const &DFit,Rcpp::List const &InParam,
                           int const &N,double const &NoisePDF){
    //------------
    // declare functions for EM steps:
    //------------
    void EMstep2D_Rcpp(Rcpp::NumericVector &sdx_g,Rcpp::NumericVector &lambdax_g,
                       Rcpp::NumericVector &mx_g, Rcpp::NumericVector &sdy_g,
                       Rcpp::NumericVector &lambday_g,Rcpp::NumericVector &my_g,
                       Rcpp::NumericVector &p_g,Rcpp::NumericMatrix const &DFit,
                       int const &G,int const &N,bool &OnlyNoise,double &CritVal,
                       double const &NoisePDF,Rcpp::NumericVector &Classification,
                       double &MixLogLik,Rcpp::NumericVector &ClassTot, int &Gupdate);
    //------------
    // Initialize:
    //------------
    // Break input parameters, you have to clone them:
    Rcpp::NumericVector sdx_g=InParam["sdx_g"];
    Rcpp::NumericVector lambdax_g=InParam["lambdax_g"];
    Rcpp::NumericVector mx_g=InParam["mx_g"];
    Rcpp::NumericVector sdy_g=InParam["sdy_g"];
    Rcpp::NumericVector lambday_g=InParam["lambday_g"];
    Rcpp::NumericVector my_g=InParam["my_g"];
    Rcpp::NumericVector p_g=InParam["p_g"];
    // Find total clusters based on if having noise or not:
    // G=total clusters/peaks, Gmax=Total clusters with noise G+1.
    int G=sdx_g.size();//total clusters.
    int Gupdate;//this saves the total number of non NAN peaks
    bool OnlyNoise;// OnlyNoise for counting if the classification returns only noise
    double MixLogLik,CritVal=INFINITY;//mixture log lik, and critical EM value.
    Rcpp::NumericVector Classification(N),ClassTot(G);
    // the classification vector, the n_g(sum of resp),
    // ClassTot=the total obervations in each cluster only for non-noise
    //------------
    // First EMstep:
    //------------
    EMstep2D_Rcpp(sdx_g,lambdax_g,mx_g,sdy_g,lambday_g,my_g,p_g,DFit,
                  G,N,OnlyNoise,CritVal,NoisePDF,Classification,
                  MixLogLik,ClassTot,Gupdate);
    //------------
    //Start EM loop:
    //------------
    int Emit=1;
    while(Emit<=500&&(CritVal>1e-6)){
        //------------
        // new EMstep:
        //------------
        EMstep2D_Rcpp(sdx_g,lambdax_g,mx_g,sdy_g,lambday_g,my_g,p_g,DFit,
                      G,N,OnlyNoise,CritVal,NoisePDF,Classification,
                      MixLogLik,ClassTot,Gupdate);
        // increase step
        Emit+=1;
    }
    //------------
    // check if Gupdate==0:
    //------------
    if(Gupdate==0) OnlyNoise=TRUE; //also only noise just in case
    //------------
    // create output
    //------------
    Rcpp::List Param=Rcpp::List::create(Rcpp::Named("p_g")=p_g,
                                        Rcpp::Named("sdx_g")=sdx_g,
                                        Rcpp::Named("lambdax_g")=lambdax_g,
                                        Rcpp::Named("mx_g")=mx_g,
                                        Rcpp::Named("sdy_g")=sdy_g,
                                        Rcpp::Named("lambday_g")=lambday_g,
                                        Rcpp::Named("my_g")=my_g);
    Rcpp::List Res=Rcpp::List::create(Rcpp::Named("Classification")=Classification,
                                      Rcpp::Named("Param")=Param,
                                      Rcpp::Named("OnlyNoise")=OnlyNoise,
                                      Rcpp::Named("ClassTot")=ClassTot);
    return Res;
}
//Done.

//---------------------
//Estep and Mstep Function for main model fit:
//---------------------
void EMstep2D_Rcpp(Rcpp::NumericVector &sdx_g,Rcpp::NumericVector &lambdax_g,
                   Rcpp::NumericVector &mx_g, Rcpp::NumericVector &sdy_g,
                   Rcpp::NumericVector &lambday_g,Rcpp::NumericVector &my_g,
                   Rcpp::NumericVector &p_g,Rcpp::NumericMatrix const &DFit,
                   int const &G,int const &N,bool &OnlyNoise,double &CritVal,
                   double const &NoisePDF,Rcpp::NumericVector &Classification,
                   double &MixLogLik,Rcpp::NumericVector &ClassTot, int &Gupdate){


    // ------------
    // Declare function:
    // ------------
    // Updating parameters:
    void Mstep2D_Rcpp(Rcpp::NumericVector &sdx_g,Rcpp::NumericVector &lambdax_g,
                      Rcpp::NumericVector &mx_g, Rcpp::NumericVector &sdy_g,
                      Rcpp::NumericVector &lambday_g,Rcpp::NumericVector &my_g,
                      Rcpp::NumericVector &p_g,int const &G,
                      int const &N,Rcpp::NumericVector const &n_g,
                      Rcpp::NumericVector const &A_my_g,
                      Rcpp::NumericVector const &B_my_g,
                      Rcpp::NumericVector const &C_my_g,
                      Rcpp::NumericVector const &RespFi_k_g,
                      Rcpp::NumericVector const &RespFiXi_k_g,
                      Rcpp::NumericMatrix const &Resp,Rcpp::NumericMatrix const &DFit);
    // Finding densities and initiating class totals, also update densities sum:
    void DensRespTot_g_Rcpp(int const &g,int const &i,
                            Rcpp::NumericVector const &sdx_g,Rcpp::NumericVector const &lambdax_g,
                            Rcpp::NumericVector const &mx_g, Rcpp::NumericVector const &sdy_g,
                            Rcpp::NumericVector const &lambday_g,Rcpp::NumericVector const &my_g,
                            Rcpp::NumericVector const &p_g,double &DensSum_i,
                            Rcpp::NumericMatrix const &DFit,
                            Rcpp::NumericMatrix &Resp,int &Gupdate,
                            Rcpp::NumericVector const &LAMBDA_g);
    // function for finding the needed sums for the parameters to update
    void FindParamSums_my_k_g_Rcpp(int const &i,int const &g,
                                   Rcpp::NumericMatrix const &DFit,Rcpp::NumericMatrix const &Resp,
                                   Rcpp::NumericVector const &sdx_g,Rcpp::NumericVector const &lambdax_g,
                                   Rcpp::NumericVector const &mx_g, Rcpp::NumericVector const &sdy_g,
                                   Rcpp::NumericVector const &lambday_g,Rcpp::NumericVector const &my_g,
                                   Rcpp::NumericVector &A_my_g,Rcpp::NumericVector &B_my_g,
                                   Rcpp::NumericVector &C_my_g,
                                   Rcpp::NumericVector &RespFi_k_g,Rcpp::NumericVector &RespFiXi_k_g);
    // function for getting the LAMBDA normalization factor for cluster g:
    double Get_LAMBDA_g_Rcpp(int const &g,Rcpp::NumericVector const &sdx_g,
                             Rcpp::NumericVector const &lambdax_g,
                             Rcpp::NumericVector const &mx_g,
                             Rcpp::NumericVector const &sdy_g,
                             Rcpp::NumericVector const &lambday_g,
                             Rcpp::NumericVector const &my_g);
    //------------
    // Initiate:
    //------------
    // Density of U/Dtags and UDTags.
    double DensSum_i;
    // DensSum_i=sum of densities for i observ over clusters
    double MixLogLikOld=MixLogLik;//for critical value(ONLY IF THIS KIND OF CONVERGENCE USED)
    MixLogLik=0;//sums the mixture log likelihood for BIC, have to start at 0.
    // Note is a cluster is with error (any param==NAN), then it is not included
    OnlyNoise=TRUE;//Counts the total non-noise obs, has to start at 0.
    // NOTE: Resp: The Resp matrix is first filled with the fUDij*p_g values
    // Then it is devided with their sum to finalize them.
    Rcpp::NumericMatrix Resp(N,G+1);
    Rcpp::NumericVector n_g(G+1),A_my_g(G),B_my_g(G),C_my_g(G),
    RespFi_k_g(G),RespFiXi_k_g(G);
    Gupdate=G;//initiate Gupdate to the total number of clusters
    Rcpp::NumericVector LAMBDA_g(G);//for the normalization factor of the densities
    //------------
    // loop and fill the rows of the matrices
    //------------
    for(int i=0;i<N;i++){
        // denine RespSum_i sum of densities:
        DensSum_i=NoisePDF*p_g[0];
        // Resp of noise before deviding with DensSum_i:
        Resp(i,0)=DensSum_i;//0 is for the noise cluster
        //------------
        // loop though parameters:
        // NOTE, having noise or not, those here will start at g=0
        // Because you have the other parameter vectors, which do not include any noise in them.
        //------------
        for(int g=0;g<G;g++){//for each peak cluster only
            // the two following are only needed once per cluster:
            if(i==0){
                // find the normalization factor for cluster g once.
                LAMBDA_g[g]=Get_LAMBDA_g_Rcpp(g,sdx_g,lambdax_g,mx_g,sdy_g,
                                              lambday_g,my_g);
                // seet the ClassTot to zero if i==0, to initiate:
                ClassTot[g]=0;//initiate only once
            }
            // Find densities and update responsibilites:
            DensRespTot_g_Rcpp(g,i,sdx_g,lambdax_g,mx_g,sdy_g,lambday_g,my_g,p_g,
                               DensSum_i,DFit,Resp,Gupdate,LAMBDA_g);

        }//run all clusters for the current i obs.
        //------------
        //now finalize responisbilites by deviding with DensSum_i for obs i:
        //------------
        double maxResp_i=-INFINITY;//max value of resp, for classification
        for(int g=0;g<(G+1);g++){//use Gmax for all the clusters here
            // devide the resp to get their final values:
            Resp(i,g)=Resp(i,g)/DensSum_i;//dont fix small numbers
            // classification:
            if(Resp(i,g)>maxResp_i){
                Classification[i]=g;//save class, given noise, so g>=0
                maxResp_i=Resp(i,g);//save current maximum
            }
            // Compute the rest of the matrices here:
            if(i==0) n_g[g]=0;//initiate only once
            n_g[g]+=Resp(i,g);
            // rest of the matrices:
            // Note this gloc is for the other matrices NOT Resp.
            // the rest are only evaluated for non-noise
            if(g!=0){
                // Find and update the sums that are needed for the parameter update
                FindParamSums_my_k_g_Rcpp(i,g,DFit,Resp,sdx_g,lambdax_g,
                                          mx_g,sdy_g,lambday_g,my_g,
                                          A_my_g,B_my_g,C_my_g,RespFi_k_g,
                                          RespFiXi_k_g);
            }
        }
        // check the classification for observation i:
        if(Classification[i]!=0){
            // save in class totals:
            ClassTot[Classification[i]-1]+=1;//the index of ClassTot is always the Classification[i]-1
            // Check OnlyNoise:
            if(ClassTot[Classification[i]-1]>=2) OnlyNoise=FALSE;//at least two obs on at least one peak
        }
        // add the MixLogLik for obervation i:
        MixLogLik+=std::log(DensSum_i);
    }
    //----------
    // find critical value(IF this is NOT used, then the converge is with the estimates):
    //----------
    CritVal=std::abs(MixLogLik-MixLogLikOld)/std::abs(MixLogLik);
    //----------
    // update parameters
    //----------
    Mstep2D_Rcpp(sdx_g,lambdax_g,mx_g,sdy_g,lambday_g,my_g,p_g,G,
                 N,n_g,A_my_g,B_my_g,C_my_g,RespFi_k_g,RespFiXi_k_g,Resp,DFit);

}
//Done

//---------------------
//Function for getting the LAMBDA normalization factor from cluster g:
//---------------------
double Get_LAMBDA_g_Rcpp(int const &g,Rcpp::NumericVector const &sdx_g,
                         Rcpp::NumericVector const &lambdax_g,
                         Rcpp::NumericVector const &mx_g,
                         Rcpp::NumericVector const &sdy_g,
                         Rcpp::NumericVector const &lambday_g,
                         Rcpp::NumericVector const &my_g){

    //----------
    // compute the 6 parts of the integral F11,F12,F21,F22,F31,F32
    //----------
    // --------- F11:
    double F11_1=(1-lambdax_g[g])*(1-lambday_g[g])/4.0;
    double F11_2=mx_g[g]-my_g[g];
    double F11_3=std::pow(F11_2,2.0);
    double F11_4=std::pow((1-lambday_g[g])*sdy_g[g],2.0);
    double F11_5=std::pow(F11_3+F11_4,0.5);
    double F11=F11_1*(F11_2/F11_5+1.0);
    // --------- F12:
    double F12_1=-2.0*(1-lambdax_g[g])/(3.0*sdy_g[g]);
    double F12_2=std::pow((1.0-lambdax_g[g])*sdx_g[g],2.0);
    double F12_3=std::pow(1.0+std::pow(F11_2-1.0,2.0)/F11_4,-1.5);
    double F12=F12_1*F12_3/std::pow(1.0+F12_2,0.5);
    // --------- F21:
    double F21=-F11_1*F11_2/std::pow(F11_3+F11_4,0.5);
    // --------- F22:
    double F22_1=(1.0+lambdax_g[g])*F11_3/(24.0*sdy_g[g]);
    double F22_2=std::pow((1.0+lambdax_g[g])*sdx_g[g],2.0);
    double F22_3=1.0/std::pow(F11_3+F22_2,0.5);
    double F22_4=4.0*std::pow(1.0+F11_3/(4.0*F11_4),-1.5);
    double F22_5=std::pow(F11_3+4.0*F22_2,0.5);
    double F22=F22_1*(F22_3+F22_4/F22_5);
    // --------- F31:
    double F31=(1-lambdax_g[g])*(1+lambday_g[g])/4.0;
    // --------- F32:
    double F32_1=(1.0+lambdax_g[g])/(4.0*sdy_g[g]);
    double F32_2=-F11_2/std::pow(F11_3+F22_2,0.5);
    double F32_3=std::pow((1+lambday_g[g])*sdy_g[g],2.0);
    double F32_4=16.0*std::pow(1.0+1.0/F32_3,-1.5)*(1-F11_2);
    double F32_5=std::pow(std::pow(1-F11_2,2.0)+F22_2,0.5);
    double F32=F32_1*(F32_2+F32_4/F32_5);
    //----------
    // Sum the parts and invert:
    //----------
    double LAMBDA=std::max(1.0/(F11+F12+F21+F22+F31+F32),0.0);

    return LAMBDA;
}

//---------------------
//Function for the densities, class totals and responsibilites:
//---------------------
void DensRespTot_g_Rcpp(int const &g,int const &i,
                        Rcpp::NumericVector const &sdx_g,Rcpp::NumericVector const &lambdax_g,
                        Rcpp::NumericVector const &mx_g, Rcpp::NumericVector const &sdy_g,
                        Rcpp::NumericVector const &lambday_g,Rcpp::NumericVector const &my_g,
                        Rcpp::NumericVector const &p_g,double &DensSum_i,
                        Rcpp::NumericMatrix const &DFit,
                        Rcpp::NumericMatrix &Resp, int &Gupdate,
                        Rcpp::NumericVector const &LAMBDA_g){

    //------------
    // find Utag_mx_gi, Dtag_my_gi and their sgns:
    //------------
    // x
    double Utag_mx_gi=DFit(i,0)-mx_g[g];
    double sgnx_g=0.0;//zero if x=mode
    if(Utag_mx_gi<0) sgnx_g=-1.0;
    if(Utag_mx_gi>0) sgnx_g=1.0;
    double DenomfUig=std::pow((1.0+sgnx_g*lambdax_g[g])*sdx_g[g],2.0);
    // y
    double Dtag_my_gi=DFit(i,1)-my_g[g];
    double sgny_g=0.0;//zero if y=mode
    if(Dtag_my_gi<0) sgny_g=-1.0;
    if(Dtag_my_gi>0) sgny_g=1.0;
    double DenomfDig=std::pow((1.0+sgny_g*lambday_g[g])*sdy_g[g],2.0);
    //------------
    // Find current densities:
    //------------
    double fUig=std::pow(1.0+std::pow(Utag_mx_gi,2.0)/DenomfUig,-1.5)/
        (2.0*sdx_g[g]);//UTag
    double fDig=std::pow(1.0+std::pow(Dtag_my_gi,2.0)/DenomfDig,-1.5)/
        (2.0*sdy_g[g]);//DTag
    double fUDig=LAMBDA_g[g]*fUig*fDig*p_g[g+1];//multiply with p_g, NB index!
    if(!(std::isnan(fUDig)||std::isinf(fUDig))){
        // save density to vector for the sum:
        DensSum_i+=fUDig;//sum of densities
        Resp(i,g+1)=fUDig;//save for responsinilites
    }else{
        Resp(i,g+1)=0;//the cluster is out
        if(i==0) Gupdate-=1;//remove it from the Gupdate too
    }//else leave it zero, as the cluster is out
}

//Done

//---------------------
//Function for the Sums that are needed for the parameters update for my and k:
//---------------------
void FindParamSums_my_k_g_Rcpp(int const &i,int const &g,
                               Rcpp::NumericMatrix const &DFit,Rcpp::NumericMatrix const &Resp,
                               Rcpp::NumericVector const &sdx_g,Rcpp::NumericVector const &lambdax_g,
                               Rcpp::NumericVector const &mx_g, Rcpp::NumericVector const &sdy_g,
                               Rcpp::NumericVector const &lambday_g,Rcpp::NumericVector const &my_g,
                               Rcpp::NumericVector &A_my_g,Rcpp::NumericVector &B_my_g,
                               Rcpp::NumericVector &C_my_g,
                               Rcpp::NumericVector &RespFi_k_g,
                               Rcpp::NumericVector &RespFiXi_k_g){

    //------------
    // find uxi_g and uyi_g:
    //------------
    // x
    double Utag_mx_gi=DFit(i,0)-mx_g[g-1];
    double Utag_mx_gi2=std::pow(Utag_mx_gi,2.0);
    double sgnx_g=0.0;//zero if x=mode
    if(Utag_mx_gi<0) sgnx_g=-1.0;
    if(Utag_mx_gi>0) sgnx_g=1.0;
    double lsigsgnx_g=std::pow((1.0+sgnx_g*lambdax_g[g-1])*sdx_g[g-1],2.0);
    // y
    double Dtag_my_gi=DFit(i,1)-my_g[g-1];
    double Dtag_my_gi2=std::pow(Dtag_my_gi,2.0);
    double sgny_g=0.0;//zero if y=mode
    if(Dtag_my_gi<0) sgny_g=-1.0;
    if(Dtag_my_gi>0) sgny_g=1.0;
    double lsigsgny_g=std::pow((1.0+sgny_g*lambday_g[g-1])*sdy_g[g-1],2.0);
    //------------
    // find weights for my:
    //------------
    double wxi_g=1.5/(lsigsgnx_g+Utag_mx_gi2);
    double wyi_g=1.5/(lsigsgny_g+Dtag_my_gi2);
    A_my_g[g-1]+=Resp(i,g)*(wxi_g*DFit(i,0)+wyi_g*DFit(i,1));
    B_my_g[g-1]+=Resp(i,g)*wxi_g;
    C_my_g[g-1]+=Resp(i,g)*(wxi_g+wyi_g);
    //------------
    // find weights for k_g:
    //------------
    double fixi_g=2.0*wxi_g;
    RespFi_k_g[g-1]+=Resp(i,g)*fixi_g;
    RespFiXi_k_g[g-1]+=Resp(i,g)*fixi_g*DFit(i,0);
}

//Done

//---------------------
//Function for the Sums that are needed for the parameters update for lambda and sigma:
//---------------------
void FindParamSums_lambda_sigma_g_Rcpp(int const &i,int const &g,
                                       Rcpp::NumericMatrix const &DFit,Rcpp::NumericMatrix const &Resp,
                                       Rcpp::NumericVector const &sdx_g,Rcpp::NumericVector const &lambdax_g,
                                       Rcpp::NumericVector const &mx_g, Rcpp::NumericVector const &sdy_g,
                                       Rcpp::NumericVector const &lambday_g,Rcpp::NumericVector const &my_g,
                                       double &Enum_sdx_g,double &Enum_sdy_g,
                                       double &Pos_lambdax_g,double &Neg_lambdax_g,
                                       double &Pos_lambday_g,double &Neg_lambday_g){

    //------------
    // find uxi_g and uyi_g:
    //------------
    // x
    double Utag_mx_gi=DFit(i,0)-mx_g[g-1];
    double Utag_mx_gi2=std::pow(Utag_mx_gi,2.0);
    double sgnx_g=0.0;//zero if x=mode
    if(Utag_mx_gi<0) sgnx_g=-1.0;
    if(Utag_mx_gi>0) sgnx_g=1.0;
    double lsigsgnx_g=std::pow((1.0+sgnx_g*lambdax_g[g-1])*sdx_g[g-1],2.0);
    double wxi_g=1.5/(lsigsgnx_g+Utag_mx_gi2);
    // y
    double Dtag_my_gi=DFit(i,1)-my_g[g-1];
    double Dtag_my_gi2=std::pow(Dtag_my_gi,2.0);
    double sgny_g=0.0;//zero if y=mode
    if(Dtag_my_gi<0) sgny_g=-1.0;
    if(Dtag_my_gi>0) sgny_g=1.0;
    double lsigsgny_g=std::pow((1.0+sgny_g*lambday_g[g-1])*sdy_g[g-1],2.0);
    double wyi_g=1.5/(lsigsgny_g+Dtag_my_gi2);
    //------------
    // find weights for sdx_g and sdy_g:
    //------------
    double zxi_g=2.0*wxi_g*Utag_mx_gi2*std::pow(sdx_g[g-1],2.0);
    Enum_sdx_g+=Resp(i,g)*zxi_g;
    double zyi_g=2.0*wyi_g*Dtag_my_gi2*std::pow(sdy_g[g-1],2.0);
    Enum_sdy_g+=Resp(i,g)*zyi_g;
    //------------
    // find weights for lambdax_g and lambday_g:
    //------------
    // x
    double Alx_g=1.5/(std::pow(Utag_mx_gi2,-1.0)+std::pow(lsigsgnx_g,-1.0));
    if(sgnx_g>0) Pos_lambdax_g+=Resp(i,g)*Alx_g;
    if(sgnx_g<0) Neg_lambdax_g+=Resp(i,g)*Alx_g;
    // y
    double Aly_g=1.5/(std::pow(Dtag_my_gi2,-1.0)+std::pow(lsigsgny_g,-1.0));
    if(sgny_g>0) Pos_lambday_g+=Resp(i,g)*Aly_g;
    if(sgny_g<0) Neg_lambday_g+=Resp(i,g)*Aly_g;
}

//Done

//---------------------
//Mstep Function for main model fit:
//---------------------
void Mstep2D_Rcpp(Rcpp::NumericVector &sdx_g,Rcpp::NumericVector &lambdax_g,
                  Rcpp::NumericVector &mx_g, Rcpp::NumericVector &sdy_g,
                  Rcpp::NumericVector &lambday_g,Rcpp::NumericVector &my_g,
                  Rcpp::NumericVector &p_g,int const &G,
                  int const &N,Rcpp::NumericVector const &n_g,
                  Rcpp::NumericVector const &A_my_g,
                  Rcpp::NumericVector const &B_my_g,
                  Rcpp::NumericVector const &C_my_g,
                  Rcpp::NumericVector const &RespFi_k_g,
                  Rcpp::NumericVector const &RespFiXi_k_g,
                  Rcpp::NumericMatrix const &Resp,Rcpp::NumericMatrix const &DFit){

    // ------------
    // Declare function:
    // ------------
    // function for finding the needed sums for the parameters to update
    void FindParamSums_lambda_sigma_g_Rcpp(int const &i,int const &g,
                                           Rcpp::NumericMatrix const &DFit,Rcpp::NumericMatrix const &Resp,
                                           Rcpp::NumericVector const &sdx_g,Rcpp::NumericVector const &lambdax_g,
                                           Rcpp::NumericVector const &mx_g, Rcpp::NumericVector const &sdy_g,
                                           Rcpp::NumericVector const &lambday_g,Rcpp::NumericVector const &my_g,
                                           double &Enum_sdx_g,double &Enum_sdy_g,
                                           double &Pos_lambdax_g,double &Neg_lambdax_g,
                                           double &Pos_lambday_g,double &Neg_lambday_g);
    // function for NR steps for lambda x and y:
    void NR_lambdax_g_Rcpp(double &lambdax,double const &Pos_lambdax_g,
                           double const &Neg_lambdax_g,double const &sdx);
    void NR_lambday_g_Rcpp(double &lambday,double const &Pos_lambday_g,
                           double const &Neg_lambday_g,double const &sdy);


    for(int g=0;g<(G+1);g++){
        //----
        // p_g:
        //----
        p_g[g]=n_g[g]/double(N);
        //----
        // rest:
        //----
        if(g!=0){
            //----
            // k_g update:
            //----
            double Ak_g=1.0/std::max(1.0-B_my_g[g-1]/C_my_g[g-1],std::sqrt(DOUBLE_EPS));
            double Bk_g=A_my_g[g-1]/C_my_g[g-1];
            double Ck_g=RespFiXi_k_g[g-1]/RespFi_k_g[g-1];
            double k_g=std::max(std::sqrt(DOUBLE_EPS),Ak_g*(Bk_g-Ck_g));
            //----
            // my_g:
            //----
            my_g[g-1]=Bk_g+std::min(B_my_g[g-1]/C_my_g[g-1],1.0-std::sqrt(DOUBLE_EPS))*k_g;
            // update mx_g:
            mx_g[g-1]=my_g[g-1]-k_g;
            // ----
            // Find constant sums for the other four parameters by updating the weights:
            // ----
            double Enum_sdx_g=0,Enum_sdy_g=0,Pos_lambdax_g=0,Neg_lambdax_g=0,
                Pos_lambday_g=0,Neg_lambday_g=0;
            for(int i=0;i<N;i++){
                FindParamSums_lambda_sigma_g_Rcpp(i,g,DFit,Resp,sdx_g,lambdax_g,mx_g,sdy_g,
                                                  lambday_g,my_g,Enum_sdx_g,Enum_sdy_g,Pos_lambdax_g,
                                                  Neg_lambdax_g,Pos_lambday_g,Neg_lambday_g);
            }
            // ----
            // lambdax_g, lambday_g:
            // ----
            // x:
            double lambdax=lambdax_g[g-1];
            double sdx=sdx_g[g-1];
            NR_lambdax_g_Rcpp(lambdax,Pos_lambdax_g,Neg_lambdax_g,sdx);
            lambdax_g[g-1]=lambdax;
            // y:
            double lambday=lambday_g[g-1];
            double sdy=sdy_g[g-1];
            NR_lambday_g_Rcpp(lambday,Pos_lambday_g,Neg_lambday_g,sdy);
            lambday_g[g-1]=lambday;
            //----
            //sdx_g, sdy_g:
            //----
            sdx_g[g-1]=std::sqrt(Enum_sdx_g/n_g[g]);
            sdy_g[g-1]=std::sqrt(Enum_sdy_g/n_g[g]);
        }
    }
}

//Done

//---------------------
//Newtown rapshon for lambdax:
//---------------------
// main function:
void NR_lambdax_g_Rcpp(double &lambdax,double const &Pos_lambdax_g, double const &Neg_lambdax_g,double const &sdx){

    // define the hyperparameters:
    double hypa=2.0, hypb=40.0;
    int Nit=1;//newton iteractions:
    double lambdax_old=lambdax, Ncrit=INFINITY;//old parameters and criterio:
    // loop
    while(Nit<50&&Ncrit>1e-6){
        // find f:
        double f1=lambdax_old*std::pow(1.0-lambdax_old,3.0)*Pos_lambdax_g;
        double f2=lambdax_old*std::pow(1.0+lambdax_old,3.0)*Neg_lambdax_g;
        double f3=std::pow(sdx,2.0)*std::pow(1.0-lambdax_old,3.0)*std::pow(1.0+lambdax_old,2.0)*
            (lambdax_old*(hypa-1.0)+(1.0+lambdax_old)*(hypb-1.0))/2.0;
        double f=f1-f2+f3;
        // find fder:
        double fder1=Pos_lambdax_g*std::pow(1.0-lambdax_old,2.0)*(1.0-4.0*lambdax_old);
        double fder2=Neg_lambdax_g*std::pow(1.0+lambdax_old,2.0)*(1.0+4.0*lambdax_old);
        double fder3=(hypa-1.0)*std::pow(1.0-lambdax_old,2.0)*(1.0+lambdax_old)*
            (1.0-6.0*std::pow(lambdax_old,2.0)-lambdax_old);
        double fder4=-6.0*(hypb-1.0)*std::pow(1.0+lambdax_old,2.0)*std::pow(1.0-lambdax_old,2.0)*lambdax_old;
        double fder=fder1-fder2+std::pow(sdx,2.0)*(fder3+fder4)/2.0;
        // update:
        lambdax=lambdax_old-f/fder;
        if(lambdax>=0.0) lambdax=-std::sqrt(DOUBLE_EPS);
        if(lambdax<=-1.0) lambdax=-1.0+std::sqrt(DOUBLE_EPS);
        // criterio:
        Ncrit=std::abs(lambdax-lambdax_old);
        // replace lambdax_old:
        lambdax_old=lambdax;
        // iterations:
        Nit+=1;
    }

}
// done

//---------------------
//Newtown rapshon for lambday:
//---------------------
// main function:
void NR_lambday_g_Rcpp(double &lambday,double const &Pos_lambday_g, double const &Neg_lambday_g,double const &sdy){

    // define the hyperparameters:
    double hypa=40.0, hypb=2.0;
    int Nit=1;//newton iteractions:
    double lambday_old=lambday, Ncrit=INFINITY;//old parameters and criterio:
    // loop
    while(Nit<50&&Ncrit>1e-6){
        // find f:
        double f1=lambday_old*std::pow(1.0-lambday_old,3.0)*Pos_lambday_g;
        double f2=lambday_old*std::pow(1.0+lambday_old,3.0)*Neg_lambday_g;
        double f3=std::pow(sdy,2.0)*std::pow(1.0-lambday_old,2.0)*std::pow(1.0+lambday_old,3.0)*
            ((1.0-lambday_old)*(hypa-1.0)-lambday_old*(hypb-1.0))/2.0;
        double f=f1-f2+f3;
        // find fder:
        double fder1=Pos_lambday_g*std::pow(1.0-lambday_old,2.0)*(1.0-4.0*lambday_old);
        double fder2=Neg_lambday_g*std::pow(1.0+lambday_old,2.0)*(1.0+4.0*lambday_old);
        double fder3=-6.0*(hypa-1.0)*std::pow(1.0+lambday_old,2.0)*std::pow(1.0-lambday_old,2.0)*lambday_old;
        double fder4=-(hypb-1.0)*std::pow(1.0+lambday_old,2.0)*(1.0-lambday_old)*
            (1.0-6.0*std::pow(lambday_old,2.0)+lambday_old);
        double fder=fder1-fder2+std::pow(sdy,2.0)*(fder3+fder4)/2.0;
        // update:
        lambday=lambday_old-f/fder;
        if(lambday>=1.0) lambday=1.0-std::sqrt(DOUBLE_EPS);
        if(lambday<=0.0) lambday=std::sqrt(DOUBLE_EPS);
        // criterio:
        Ncrit=std::abs(lambday-lambday_old);
        // replace lambday_old:
        lambday_old=lambday;
        // iterations:
        Nit+=1;
    }

}

// done

//####################################################
//--------------------- Merging functions:
//####################################################
//---------------------
//Main Function for Merging overlapping peaks
//---------------------
void MergeOvPeak_Rcpp(Rcpp::NumericMatrix const &DFit,
                      int const &N,double const &NoisePDF,
                      Rcpp::List &OptParam,Rcpp::NumericVector &OptClass,
                      Rcpp::NumericVector &OptClassTot,
                      Rcpp::NumericVector const &KernSeq){

    //------------
    // declare functions:
    //------------
    // Function for getting the Merging info:
    Rcpp::List Get_MergeClusterInf_Rcpp(Rcpp::List const &OptParam,
                                        double const &window,Rcpp::NumericVector const &OptClassTot);
    // function for initialing when merging:
    Rcpp::List Get_MergeInitials_Rcpp(Rcpp::List const &ClstList,
                                      Rcpp::List const &OptParam,
                                      int const &TotAfterMerged,
                                      Rcpp::NumericVector const &OptClass,
                                      int const &N,Rcpp::NumericMatrix const &DFit,
                                      Rcpp::NumericVector const &KernSeq);
    // Function for running EM for SGT Model
    Rcpp::List MainEMLoop_Rcpp(Rcpp::NumericMatrix const &DFit,
                               Rcpp::List const &InParam_m,int const &N,
                               double const &NoisePDF);
    //------------
    // Initiate common variables:
    //------------
    double window=300;//the merging window
    //------------
    // Find first Merging info and loop:
    //------------
    Rcpp::List MergeInf;
    MergeInf=Get_MergeClusterInf_Rcpp(OptParam,window,OptClassTot);
    // loop and merge iff overlaps:
    int Merge_safe_it=1;//for not getting in infinete loop
    while(Rcpp::as<int>(MergeInf["TotOv"])!=0&&Merge_safe_it<=1000){
        //------------
        // Initiate EM:
        //------------
        // finding new initials, unmerged clusters are given the estimated values
        Rcpp::List InitialsMerge=Get_MergeInitials_Rcpp(MergeInf["ClstList"],
                                                        OptParam,MergeInf["TotAfterMerged"],
                                                                         OptClass,N,DFit,KernSeq);
        //------------
        // Fit model:
        //------------
        Rcpp::List LocRes=MainEMLoop_Rcpp(DFit,InitialsMerge,N,NoisePDF);
        //------------
        // save only if allowed:
        //------------
        if(!Rcpp::as<bool>(LocRes["OnlyNoise"])){
            // Replace the old model with the new one:
            OptParam=LocRes["Param"];
            OptClass=LocRes["Classification"];
            OptClassTot=LocRes["ClassTot"];
            //------------
            // Again find mergeinf test further merging:
            //------------
            MergeInf=Get_MergeClusterInf_Rcpp(OptParam,window,OptClassTot);
        }else{
            break;//break the while and dont merge further, without saving.
        }
        Merge_safe_it+=1;
    }
}
//Done
//---------------------
//Function for getting Merging Cluster Info, clusters to merge etc
//---------------------
Rcpp::List Get_MergeClusterInf_Rcpp(Rcpp::List const &OptParam,
                                    double const &window,Rcpp::NumericVector const &OptClassTot){

    //---------
    //Declare functions:
    //---------
    //for finding pairs of overlaps
    void Get_PairsMat_Rcpp(Rcpp::List const &OptParam,
                           Rcpp::NumericMatrix &PairsMat,
                           int const &G,double const &window,
                           int &TotOv,Rcpp::NumericVector const &OptClassTot);
    //for finding the cluster list.
    void Get_ClstList_Rcpp(Rcpp::List &ClstList,int const &G,
                           Rcpp::NumericMatrix const &PairsMat);
    //---------
    //Initialize output:
    //---------
    Rcpp::List Res=Rcpp::List::create(Rcpp::Named("TotOv")=0,//Total overlaps found
                                      Rcpp::Named("ClstList")=R_NilValue,//List with the cluster IDS to be merged
                                      Rcpp::Named("TotAfterMerged")=R_NilValue);//Total clusters left after merging.
    //---------
    //Enter iff G>1:, alse no cluster to merge
    //---------
    int G=Rcpp::as<Rcpp::NumericVector>(OptParam["sdx_g"]).size();//peaks tot
    if(G>1){
        //---------
        //Find Pairs matrix
        //---------
        // Initiate for the PairsMatrix:
        Rcpp::NumericMatrix PairsMat(G,G);//indicator matrix for overlapping peaks
        int TotOv=0;//number of overlaps found
        //---------
        // call function for getting the PairsMat matrix:
        Get_PairsMat_Rcpp(OptParam,PairsMat,G,window,TotOv,OptClassTot);
        // If no overlaps or all NA or 1 ok rest NA, totov=0.
        //---------
        //Find cluster list:
        //---------
        if(TotOv!=0){
            //---------
            // initiate for the cluster list:
            Rcpp::List ClstList;//The cluster list
            //---------
            // Call function for cluster list:
            Get_ClstList_Rcpp(ClstList,G,PairsMat);
            //---------
            //Update output:
            //---------
            Res["TotOv"]=TotOv;
            Res["ClstList"]=ClstList;
            Res["TotAfterMerged"]=ClstList.size();
        }
    }
    return Res;
}
//Done

//---------------------
//Function for getting the PairsMatrix for merging opperations:
//---------------------
void Get_PairsMat_Rcpp(Rcpp::List const &OptParam,
                       Rcpp::NumericMatrix &PairsMat,
                       int const &G,double const &window,
                       int &TotOv,Rcpp::NumericVector const &OptClassTot){
    // --------------
    // loop and fill in the matrix:
    // --------------
    for(int gr=0;gr<G;gr++){
        // take current mx_g, my_g::
        double mx_gr=Rcpp::as<Rcpp::NumericVector>(OptParam["mx_g"])[gr];
        double my_gr=Rcpp::as<Rcpp::NumericVector>(OptParam["my_g"])[gr];
        double peak_gr=std::round((mx_gr+my_gr)/2.0);
        if(OptClassTot[gr]<2||std::isnan(peak_gr)) continue;//go to next, leave it zero, to be removed from list after
        // else test if overlap:
        for(int gc=gr;gc<G;gc++){
            // take for gc
            double mx_gc=Rcpp::as<Rcpp::NumericVector>(OptParam["mx_g"])[gc];
            double my_gc=Rcpp::as<Rcpp::NumericVector>(OptParam["my_g"])[gc];
            double peak_gc=std::round((mx_gc+my_gc)/2.0);
            if(OptClassTot[gc]<2||std::isnan(peak_gc)) continue;//skip it.
            // else test overlap:
            bool Overlap=std::abs(peak_gc-peak_gr)<=window;
            if(Overlap){
                PairsMat(gr,gc)=1;
                if(gr!=gc) TotOv+=1;//count the overlap, not self.
            }
        }
    }
}
//Done

//---------------------
//Function for getting the ClustList from the PairsMatrix for merging opperations:
//---------------------

void Get_ClstList_Rcpp(Rcpp::List &ClstList,int const &G,
                       Rcpp::NumericMatrix const &PairsMat){

    // ClstList is empty, push elements in it.
    // scan row:
    for(int gr=0;gr<G;gr++){
        // create row non zero elements:
        Rcpp::NumericVector gr_elmt;
        // test if the diagonal==1:
        if(PairsMat(gr,gr)==0) continue;
        // find non zero elements of that row:
        for(int gc=gr;gc<G;gc++){
            if(PairsMat(gr,gc)==1) gr_elmt.push_back(gc);
        }
        // now you have the non zero elements of that row.
        // scan the list to see if the elements belong to any already
        int ClstList_size=ClstList.size();
        if(ClstList_size==0){
            // then add the elements, the list is empty:
            ClstList.push_back(gr_elmt);
        }else{
            // scan the list elements to see if it exists somewhere:
            bool Saved=FALSE;//it any element in gr_elmt is saved in the list
            int SavedPos;//where this element is saved in the list
            for(int lst=0;lst<ClstList_size;lst++){
                Rcpp::NumericVector ClstList_lst=ClstList[lst];
                for(int i=0;i<gr_elmt.size();i++){//elements in vector
                    for(int j=0;j<ClstList_lst.size();j++){//elements in vector of list
                        if(gr_elmt[i]==ClstList_lst[j]){
                            // then inside there
                            Saved=TRUE;
                            SavedPos=lst;//the list element
                        }
                        if(Saved) break;
                    }
                    if(Saved) break;
                }
                if(Saved) break;
            }

            // now check if it is saved or not:
            if(!Saved){
                // then all elements new so save
                ClstList.push_back(gr_elmt);
            }else{
                // else any element in gr_elmt is saved in ClstList[SavedPos]
                Rcpp::NumericVector ClstList_SavedPos=ClstList[SavedPos];
                // save the elements in they are not there:
                for(int i=0;i<gr_elmt.size();i++){//elements in vector
                    bool Element_Saved=FALSE;
                    for(int j=0;j<ClstList_SavedPos.size();j++){//elements in vector of list
                        if(gr_elmt[i]==ClstList_SavedPos[j]){
                            // then inside there
                            Element_Saved=TRUE;
                            break;
                        }
                    }
                    // if element not in then save it:
                    if(!Element_Saved) ClstList_SavedPos.push_back(gr_elmt[i]);
                }
                // replace list element:
                ClstList[SavedPos]=ClstList_SavedPos;
            }
        }
    }
}
//Done
//---------------------
//Function for initializing the Merging parameters:
//---------------------
Rcpp::List Get_MergeInitials_Rcpp(Rcpp::List const &ClstList,
                                  Rcpp::List const &OptParam,
                                  int const &TotAfterMerged,
                                  Rcpp::NumericVector const &OptClass,
                                  int const &N,Rcpp::NumericMatrix const &DFit,
                                  Rcpp::NumericVector const &KernSeq){
    //----------
    // declare functions:
    //----------
    //function for kernels and local maxima:
    Rcpp::NumericMatrix MatchKernelsPairs_Rcpp(Rcpp::NumericVector const &KernSeq,
                                               Rcpp::NumericMatrix const &DFit,
                                               Rcpp::NumericVector const &FeatureID,
                                               int const &Nsub,
                                               int const &KernSeqLength,
                                               double const &BW,
                                               double const &pi);
    //declare function for initiating:
    Rcpp::List InitializeKernels_Rcpp(Rcpp::NumericMatrix const &KernelPairs,
                                      Rcpp::NumericMatrix const &DFit,int const &N);

    //----------
    // Initialize:
    //----------
    //total clusters before noise. And after noise(to be updated), and KernSeqLength.
    int G=TotAfterMerged,KernSeqLength=KernSeq.size();
    const double pi=3.14159265358979323846;//for kernels
    //break the OptParam List:
    Rcpp::NumericVector Merge_mx_g=OptParam["mx_g"];
    Rcpp::NumericVector Merge_my_g=OptParam["my_g"];
    Rcpp::NumericVector Merge_sdx_g=OptParam["sdx_g"];
    Rcpp::NumericVector Merge_sdy_g=OptParam["sdy_g"];
    Rcpp::NumericVector Merge_lambdax_g=OptParam["lambdax_g"];
    Rcpp::NumericVector Merge_lambday_g=OptParam["lambday_g"];
    Rcpp::NumericVector Merge_p_g=OptParam["p_g"];
    // Iniate parameters for the merging:
    Rcpp::NumericVector p_g(G+1),sdx_g(G),lambdax_g(G),mx_g(G),sdy_g(G),
    lambday_g(G),my_g(G);
    // Fill in the noise p_g if any:
    p_g[0]=Merge_p_g[0];
    //----------
    // Loop and fill in the values:
    //----------
    for(int g=0;g<G;g++){//runs only of the clusters, not noise. clusters of the list
        //----------
        // Take list element:
        //----------
        Rcpp::NumericVector ClstList_g=ClstList[g];//list element==clusters to merge
        // Check if single element=> replace the estimates
        // Or if merging clusters=> find Kernels again
        int ClstList_g_size=ClstList_g.size();
        if(ClstList_g_size==1){
            //----------
            //then one cluster, not merging, replace the parameters
            // ClstList_g refers to the index in the old parameters, with noise or not
            // NOTE, I have to have [0], because the ClstList_g is a vector
            //----------
            mx_g[g]=Merge_mx_g[ClstList_g[0]];
            my_g[g]=Merge_my_g[ClstList_g[0]];
            sdx_g[g]=Merge_sdx_g[ClstList_g[0]];
            sdy_g[g]=Merge_sdy_g[ClstList_g[0]];
            lambdax_g[g]=Merge_lambdax_g[ClstList_g[0]];
            lambday_g[g]=Merge_lambday_g[ClstList_g[0]];
            p_g[g+1]=Merge_p_g[ClstList_g[0]];

        }else{
            //----------
            // We have a merging cluster
            // Find the FeatureID_g and Nsub_g:
            //----------
            int Nsub_g=0;
            Rcpp::NumericVector FeatureID_g;
            for(int i=0;i<N;i++){//for each observation.
                bool Found_i=FALSE;
                for(int gi=0;gi<ClstList_g_size;gi++){//for each cluster to be merged
                    // check if observation is in the class
                    if(OptClass[i]==(ClstList_g[gi]+1)){
                        // append the observation:
                        FeatureID_g.push_back(i);
                        // cound:
                        Nsub_g+=1;
                        Found_i=TRUE;
                    }
                    if(Found_i) break;//go to next i
                }
            }
            //----------
            // Now you have everything you need for running kernels:
            // call function for kernel local maxima:
            //----------
            Rcpp::NumericMatrix KernelPairs=MatchKernelsPairs_Rcpp(KernSeq,DFit,
                                                                   FeatureID_g,Nsub_g,
                                                                   KernSeqLength,20.0,pi);
            //----------
            // keep only the Kernel with the max density:
            //----------
            Rcpp::NumericMatrix KernelPairs_g(1,3);
            KernelPairs_g(0,2)=-INFINITY;//set density to -INFINITY
            for(int j=0;j<KernelPairs.nrow();j++){
                if(KernelPairs(j,2)>KernelPairs_g(0,2)){
                    KernelPairs_g(0,2)=KernelPairs(j,2);//replace density
                    KernelPairs_g(0,0)=KernelPairs(j,0);//UTag
                    KernelPairs_g(0,1)=KernelPairs(j,1);//DTag
                }
            }
            //----------
            // Now Initiate a single cluster:
            //----------
            // NOTE: the following will always return noise in p_g, skip the p_g estimates anyway
            Rcpp::List Initials_g=InitializeKernels_Rcpp(KernelPairs_g,DFit,N);
            //----------
            // Replace the new estimates:
            //----------
            // always take the 0 place here
            mx_g[g]=Rcpp::as<Rcpp::NumericVector>(Initials_g["mx_g"])[0];
            my_g[g]=Rcpp::as<Rcpp::NumericVector>(Initials_g["my_g"])[0];
            sdx_g[g]=Rcpp::as<Rcpp::NumericVector>(Initials_g["sdx_g"])[0];
            sdy_g[g]=Rcpp::as<Rcpp::NumericVector>(Initials_g["sdy_g"])[0];
            lambdax_g[g]=Rcpp::as<Rcpp::NumericVector>(Initials_g["lambdax_g"])[0];
            lambday_g[g]=Rcpp::as<Rcpp::NumericVector>(Initials_g["lambday_g"])[0];
            p_g[g+1]=double(Nsub_g)/double(N);
        }

    }
    // -------------
    // create output:
    // -------------
    Rcpp::List Res=Rcpp::List::create(Rcpp::Named("p_g")=p_g,
                                      Rcpp::Named("sdx_g")=sdx_g,
                                      Rcpp::Named("lambdax_g")=lambdax_g,
                                      Rcpp::Named("mx_g")=mx_g,
                                      Rcpp::Named("sdy_g")=sdy_g,
                                      Rcpp::Named("lambday_g")=lambday_g,
                                      Rcpp::Named("my_g")=my_g);
    return Res;
}

//Done
//####################################################
//--------------------- Peak Info Functions:
//####################################################
//---------------------
//Function for getting the peak info of the best model choosen by BIC:
//---------------------
Rcpp::DataFrame GetPeakInf_Rcpp(Rcpp::List const &OptParam,
                                Rcpp::NumericVector const &OptClassTot,
                                int const &ChromSize, int const &Region,
                                std::string const &Chrom,Rcpp::NumericVector &OldPeakNames){

    //---------
    // declare function:
    //---------
    void GetQuantilesCI_Rcpp(double const &sdx, double const &lambdax, double const &mx,
                             double const &sdy,double const &lambday, double const &my,
                             Rcpp::NumericMatrix &GlobalPeaksInfo, int const &PeakID,
                             int const &ChromSize);
    //---------
    //intialize :
    //---------
    // Break input:
    Rcpp::NumericVector sdx_g=OptParam["sdx_g"],lambdax_g=OptParam["lambdax_g"],mx_g=OptParam["mx_g"],
                                       sdy_g=OptParam["sdy_g"],lambday_g=OptParam["lambday_g"],my_g=OptParam["my_g"];
    int G=sdx_g.size(),Gused=0;//total peaks and total used peaks
    //---------
    // Find used clusters with non-zero totals:
    //---------
    for(int g=0;g<G;g++){
        if(OptClassTot[g]>=2&&!std::isnan(sdx_g[g])&&
           !std::isnan(lambdax_g[g])&&!std::isnan(mx_g[g])&&
           !std::isnan(sdy_g[g])&&!std::isnan(lambday_g[g])&&
           !std::isnan(my_g[g])) Gused+=1;
    }
    // break if no used peak:
    if(Gused==0) return R_NilValue;//(it will return an empty 0x0 dataframe)
    //---------
    // initiate the output matrix:
    //---------
    Rcpp::NumericMatrix GlobalPeaksInfo(Gused,17);
    int PeakID=0;//PeakID(for peak name and matrix index)
    for(int g=0;g<G;g++){//loop to all and skip the empty
        if(OptClassTot[g]>=2&&!std::isnan(sdx_g[g])&&
           !std::isnan(lambdax_g[g])&&!std::isnan(mx_g[g])&&
           !std::isnan(sdy_g[g])&&!std::isnan(lambday_g[g])&&
           !std::isnan(my_g[g])){//non-empty:NOTE: This removes both empty clusters and NAN clusters.

            GlobalPeaksInfo(PeakID,0)=Region;
            GlobalPeaksInfo(PeakID,1)=PeakID+1;
            GlobalPeaksInfo(PeakID,2)=OptClassTot[g];
            GlobalPeaksInfo(PeakID,3)=std::round((mx_g[g]+my_g[g])/2.0);//Peak.Summit
            GlobalPeaksInfo(PeakID,4)=mx_g[g];//Up.Summit
            GlobalPeaksInfo(PeakID,5)=my_g[g];//Down.Summit
            // Find Quantiles Info:
            GetQuantilesCI_Rcpp(sdx_g[g],lambdax_g[g],mx_g[g],sdy_g[g],
                                lambday_g[g],my_g[g],GlobalPeaksInfo,PeakID,
                                ChromSize);
            // Add the 4 other parameters, note that mx and my are named as up/down summits:
            GlobalPeaksInfo(PeakID,13)=sdx_g[g];
            GlobalPeaksInfo(PeakID,14)=lambdax_g[g];
            GlobalPeaksInfo(PeakID,15)=sdy_g[g];
            GlobalPeaksInfo(PeakID,16)=lambday_g[g];

            // save old peak name to fix the classes:
            OldPeakNames.push_back(g+1);
            // increase PeakID and index:
            PeakID+=1;
        }
    }
    // Create the data frame out put:
    Rcpp::StringVector ChromVect(Gused,Chrom);
    Rcpp::DataFrame Res=Rcpp::DataFrame::create(Rcpp::Named("Chrom")=ChromVect,
                                                Rcpp::Named("Region")=GlobalPeaksInfo(_,0),
                                                Rcpp::Named("Peak")=GlobalPeaksInfo(_,1),
                                                Rcpp::Named("Pets")=GlobalPeaksInfo(_,2),
                                                Rcpp::Named("Peak.Summit")=GlobalPeaksInfo(_,3),
                                                Rcpp::Named("Up.Summit")=GlobalPeaksInfo(_,4),
                                                Rcpp::Named("Down.Summit")=GlobalPeaksInfo(_,5),
                                                Rcpp::Named("CIQ.Up.start")=GlobalPeaksInfo(_,6),
                                                Rcpp::Named("CIQ.Up.end")=GlobalPeaksInfo(_,7),
                                                Rcpp::Named("CIQ.Up.size")=GlobalPeaksInfo(_,8),
                                                Rcpp::Named("CIQ.Down.start")=GlobalPeaksInfo(_,9),
                                                Rcpp::Named("CIQ.Down.end")=GlobalPeaksInfo(_,10),
                                                Rcpp::Named("CIQ.Down.size")=GlobalPeaksInfo(_,11),
                                                Rcpp::Named("CIQ.Peak.size")=GlobalPeaksInfo(_,12),
                                                Rcpp::Named("sdx")=GlobalPeaksInfo(_,13),
                                                Rcpp::Named("lambdax")=GlobalPeaksInfo(_,14),
                                                Rcpp::Named("sdy")=GlobalPeaksInfo(_,15),
                                                Rcpp::Named("lambday")=GlobalPeaksInfo(_,16));

    return Res;

}

//Done

//---------------------
//Function for getting the Quantiles Info for each cluster:
//---------------------
void GetQuantilesCI_Rcpp(double const &sdx, double const &lambdax, double const &mx,
                         double const &sdy,double const &lambday, double const &my,
                         Rcpp::NumericMatrix &GlobalPeaksInfo, int const &PeakID,
                         int const &ChromSize){

    double SQRTx,SQRTy;
    //--------------
    //Upper Tags quantile:
    //--------------
    double UpPcheck=(1.0-lambdax)/2.0;
    // compute 0.05 quantile:
    if(0.05<=UpPcheck){
        SQRTx=std::pow(std::pow(2.0*0.05/(1.0-lambdax)-1.0,-2.0)-1.0,-0.5);
        GlobalPeaksInfo(PeakID,6)=mx-SQRTx*(1.0-lambdax)*sdx;
    }else{
        SQRTx=std::pow(std::pow(0.05-(1.0-lambdax)/2.0,-2.0)-4.0/std::pow(1.0+lambdax,2.0),-0.5);
        GlobalPeaksInfo(PeakID,6)=mx+SQRTx*2.0*sdx;
    }
    //fix border.
    if(GlobalPeaksInfo(PeakID,6)<1) GlobalPeaksInfo(PeakID,6)=1;
    // compute 0.95 quantile:
    if(0.95<=UpPcheck){
        SQRTx=std::pow(std::pow(2.0*0.95/(1.0-lambdax)-1.0,-2.0)-1.0,-0.5);
        GlobalPeaksInfo(PeakID,7)=mx-SQRTx*(1.0-lambdax)*sdx;
    }else{
        SQRTx=std::pow(std::pow(0.95-(1.0-lambdax)/2.0,-2.0)-4.0/std::pow(1.0+lambdax,2.0),-0.5);
        GlobalPeaksInfo(PeakID,7)=mx+SQRTx*2.0*sdx;
    }
    //fix border:
    if(GlobalPeaksInfo(PeakID,7)>ChromSize) GlobalPeaksInfo(PeakID,7)=ChromSize;
    // size:
    GlobalPeaksInfo(PeakID,8)=std::round(GlobalPeaksInfo(PeakID,7))-
        std::round(GlobalPeaksInfo(PeakID,6))+1;
    //--------------
    //Down Tags quantile:
    //--------------
    double DownPcheck=(1.0-lambday)/2.0;
    // compute 0.05 quantile:
    if(0.05<=DownPcheck){
        SQRTy=std::pow(std::pow(2.0*0.05/(1.0-lambday)-1.0,-2.0)-1.0,-0.5);
        GlobalPeaksInfo(PeakID,9)=my-SQRTy*(1.0-lambday)*sdy;
    }else{
        SQRTy=std::pow(std::pow(0.05-(1.0-lambday)/2.0,-2.0)-4.0/std::pow(1.0+lambday,2.0),-0.5);
        GlobalPeaksInfo(PeakID,9)=my+SQRTy*2.0*sdy;
    }
    //fix border.
    if(GlobalPeaksInfo(PeakID,9)<1) GlobalPeaksInfo(PeakID,9)=1;
    // compute 0.95 quantile:
    if(0.95<=DownPcheck){
        SQRTy=std::pow(std::pow(2.0*0.95/(1.0-lambday)-1.0,-2.0)-1.0,-0.5);
        GlobalPeaksInfo(PeakID,10)=my-SQRTy*(1.0-lambday)*sdy;
    }else{
        SQRTy=std::pow(std::pow(0.95-(1.0-lambday)/2.0,-2.0)-4.0/std::pow(1.0+lambday,2.0),-0.5);
        GlobalPeaksInfo(PeakID,10)=my+SQRTy*2.0*sdy;
    }
    //fix border.
    if(GlobalPeaksInfo(PeakID,10)>ChromSize) GlobalPeaksInfo(PeakID,10)=ChromSize;
    // size:
    GlobalPeaksInfo(PeakID,11)=std::round(GlobalPeaksInfo(PeakID,10))
        -std::round(GlobalPeaksInfo(PeakID,9))+1;
    // Find the peak size:
    GlobalPeaksInfo(PeakID,12)=std::round(GlobalPeaksInfo(PeakID,10))-
    std::round(GlobalPeaksInfo(PeakID,6))+1;
}

//Done
