#include <Rcpp.h>
using namespace Rcpp;

//####################################################
//Functions for filtering the fastq files
//####################################################
//---------------------
//Main function for filtering: (exported in R)
//---------------------
// [[Rcpp::export]]
SEXP FilterFastqYield_fun_Rcpp(int &Curfastqyieldsize,
                               std::vector<std::string> &SreadFastq1,
                               Rcpp::NumericVector &WidthSreadFastq1,
                               std::vector<std::string> &SreadFastq2,
                               Rcpp::NumericVector &WidthSreadFastq2,
                               std::string &S0_LinkerA,
                               std::string &S0_LinkerB,
                               int &S0_LinkerOccurence,
                               int &S0_MinReadLength,
                               int &S0_MaxReadLength){

    //---------------------
    //Function declaration:
    //---------------------
    void Get_linker_type_fun_Rcpp(int const &LinkAPos,int const &LinkBPos,
                                  int &LinkType,int &LinkPos);
    // ---------------------------------------
    // --------Prepare the data for filtering:
    // ---------------------------------------
    // Create output vector: 0 means ambiguous, 1 usable and 2 chimeric, -1 NNs
    // Reads too sort/long are ambiguous at once.
    // Reads with different linkers are chimeric
    // Reads with same linker and correct size are usable,
    // If S0_LinkerOccurence=0, reads with any adaptor missing become
    // ambiguous.
    // If S0_LinkerOccurence=1, reads with read 1 linker missing but read 2
    // linker in place move to usable
    // If S0_LinkerOccurence=2, the other way around
    // If S0_LinkerOccurence=3, reads with no linker in any read become usable.
    // If S0_LinkerOccurence=4, no linker in BOTH reads.
    Rcpp::NumericVector FilterClasses(Curfastqyieldsize);
    // vector for narrowing positions:
    Rcpp::NumericMatrix NarrowingPos(Curfastqyieldsize,2);
    // counters:
    int TotusableYield=0, TotchimericYield=0,TotambiYield=0,TotNsYield=0;
    std::string NNs ("N");
    // ---------------------------------------
    // --------loop to identify adaptors:
    // ---------------------------------------
    // For everything, remember that c++ starts at 0
    for(int i=0;i<Curfastqyieldsize;i++){
        // ---------------------------------------
        // First check if any read includes NNs
        // ---------------------------------------
        if(SreadFastq1[i].find(NNs)!=-1||SreadFastq2[i].find(NNs)!=-1){
            TotNsYield+=1;
            FilterClasses[i]=-1;
            NarrowingPos(i,0)=WidthSreadFastq1[i];
            NarrowingPos(i,1)=WidthSreadFastq2[i];
            continue;
        }
        // ---------------------------------------
        // Check the linker positions. If double linkers are
        // occured keep the leftmost one.
        // Its position specifies the trimming too.
        // ---------------------------------------
        int LinkAPos_1=SreadFastq1[i].find(S0_LinkerA);
        int LinkAPos_2=SreadFastq2[i].find(S0_LinkerA);
        int LinkBPos_1=SreadFastq1[i].find(S0_LinkerB);
        int LinkBPos_2=SreadFastq2[i].find(S0_LinkerB);
        // define linker types for the reads:
        // 0, no linker, 1 linker A, 2 linker B
        // And the final cut positions.
        int LinkType_1, LinkType_2, LinkPos_1, LinkPos_2;
        // ---------------------------------------
        // Find the type of the linker and save its position:
        // ---------------------------------------
        // First read:
        Get_linker_type_fun_Rcpp(LinkAPos_1,LinkBPos_1,LinkType_1,LinkPos_1);
        // second read:
        Get_linker_type_fun_Rcpp(LinkAPos_2,LinkBPos_2,LinkType_2,LinkPos_2);
        // ---------------------------------------
        // If the linker exists in the reads, check their lengths after trimming.
        // Those who do not pass are ambiguous anyway:
        // ---------------------------------------
        bool BorderCheck_1, BorderCheck_2;
        if(LinkPos_1!=-1){
            // linker there:
            BorderCheck_1=(LinkPos_1<S0_MinReadLength||
                LinkPos_1>S0_MaxReadLength);
        }else{
            // linker missing:
            BorderCheck_1=(WidthSreadFastq1[i]<S0_MinReadLength||
                WidthSreadFastq1[i]>S0_MaxReadLength);
        }
        if(LinkPos_2!=-1){
            // linker there:
            BorderCheck_2=(LinkPos_2<S0_MinReadLength||
                LinkPos_2>S0_MaxReadLength);
        }else{
            // linker missing:
            BorderCheck_2=(WidthSreadFastq2[i]<S0_MinReadLength||
                WidthSreadFastq2[i]>S0_MaxReadLength);
        }

        if(BorderCheck_1||BorderCheck_2){
            // then the reads have the linkers, but the borders feil so ambiguous:
            // go to next, they are already to 0
            TotambiYield+=1;
            NarrowingPos(i,0)=WidthSreadFastq1[i];
            NarrowingPos(i,1)=WidthSreadFastq2[i];
            continue;
        }
        // ---------------------------------------
        // Compare the linker occurences:
        // ---------------------------------------
        if((LinkType_1==1&&LinkType_2==1)||(LinkType_1==2&&LinkType_2==2)){
            // then same linkers, so usable:
            TotusableYield+=1;
            FilterClasses[i]=1;
            NarrowingPos(i,0)=LinkPos_1;
            NarrowingPos(i,1)=LinkPos_2;
        }else if((LinkType_1==1&&LinkType_2==2)||(LinkType_1==2&&LinkType_2==1)){
            // then definetely chimeric:
            FilterClasses[i]=2;
            TotchimericYield+=1;
            NarrowingPos(i,0)=WidthSreadFastq1[i];
            NarrowingPos(i,1)=WidthSreadFastq2[i];
        }else if(LinkType_1==0&&LinkType_2!=0&&
            (S0_LinkerOccurence==1||S0_LinkerOccurence==3)){
            // then linker 1 missing, but 2 there, and keep them as usable:
            TotusableYield+=1;
            FilterClasses[i]=1;
            NarrowingPos(i,0)=WidthSreadFastq1[i];
            NarrowingPos(i,1)=LinkPos_2;

        }else if(LinkType_1!=0&&LinkType_2==0&&
            (S0_LinkerOccurence==2||S0_LinkerOccurence==3)){
            // then linker 2 missing, but 1 there, and keep them as usable:
            TotusableYield+=1;
            FilterClasses[i]=1;
            NarrowingPos(i,0)=LinkPos_1;
            NarrowingPos(i,1)=WidthSreadFastq2[i];

        }else if(LinkType_1==0&&LinkType_2==0&&
            (S0_LinkerOccurence==3||S0_LinkerOccurence==4)){
            // both missing, but keep them as usable:
            TotusableYield+=1;
            FilterClasses[i]=1;
            NarrowingPos(i,0)=WidthSreadFastq1[i];
            NarrowingPos(i,1)=WidthSreadFastq2[i];
        }else{
            // All the other categories are ambiguous:
            TotambiYield+=1;
            NarrowingPos(i,0)=WidthSreadFastq1[i];
            NarrowingPos(i,1)=WidthSreadFastq2[i];
        }
    }

    Rcpp::List Res=Rcpp::List::create(Rcpp::Named("NarrowingPos")=NarrowingPos,
                                      Rcpp::Named("FilterClasses")=FilterClasses,
                                      Rcpp::Named("TotusableYield")=TotusableYield,
                                      Rcpp::Named("TotchimericYield")=TotchimericYield,
                                      Rcpp::Named("TotambiYield")=TotambiYield,
                                      Rcpp::Named("TotNsYield")=TotNsYield);

    return Rcpp::wrap(Res);
}
// done

//---------------------
//Function for choosing the linker type and its position:
//---------------------
void Get_linker_type_fun_Rcpp(int const &LinkAPos,int const &LinkBPos,
                              int &LinkType,int &LinkPos){

    if(LinkAPos==-1&&LinkBPos==-1){
        // No linkers in read:
        LinkType=0;
        LinkPos=-1;
    }else if(LinkAPos!=-1&&LinkBPos==-1){
        // then only linker A:
        LinkType=1;
        LinkPos=LinkAPos;
    }else if(LinkAPos==-1&&LinkBPos!=-1){
        // then only linker B:
        LinkType=2;
        LinkPos=LinkBPos;
    }else if(LinkAPos<=LinkBPos){
        // then both there, choose A:
        LinkType=1;
        LinkPos=LinkAPos;
    }else{
        // then both there, choose B:
        LinkType=2;
        LinkPos=LinkBPos;
    }
}
// done
