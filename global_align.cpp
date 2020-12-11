/* @file global_align.cpp
*  @brief This script is an example of global sequence alignment with dynamic programming.*  
*  It was developed with C++ and depended on Eigen library. And It refer the implement in course_bioinfo_training of
*  GuangChuang Yu (https://github.com/GuangchuangYu/course_bioinfo_training).
*  
*  So user must install Eigen library firstly with the following code on Ubuntu.*  
*  ```
*  sudo apt-get install libeigen3-dev
*  ```  
*  Second, it also need be compiled before run it.
*  ```
*  g++ global_align.cpp -o global_align
*  ```
*  Finally, you can view the help information with the following code. 
*  ```
*  ./global_align -h
*  ```
*
*  @author Shuangbin Xu
*
*  @date 12/06/2020
*/
#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include <eigen3/Eigen/Eigen>
using namespace std;

class GlobalAlign
{
public:
    string seq1;
    string seq2;
    array<int, 3> score;
    Eigen::MatrixXd matrix;
    string aln1;
    string aln2;
    vector<int> gap_position;
    vector<int> mis_position;
    void show();
};

void GlobalAlign::show(){
    int m = seq1.length();
    int n = seq2.length();
    int l = aln1.length();
    printf("Sequence X:   %s\n", seq1.c_str());
    printf("Sequence Y:   %s\n", seq2.c_str());   
    printf("Scoring system:  %d for match; %d for mismatch; %d for gap.\n", score[0], score[1], score[2]);
    printf("\n");
    printf("Dynamic programming matrix:\n");
    printf("\n");
    cout << matrix << endl;
    string aln = ""; 
    for (int i=0; i<l; i ++){
        if (count(gap_position.begin(), gap_position.end(), i)){
            aln += " ";
        }else if (count(mis_position.begin(), mis_position.end(), i)){
            aln += "?";
        }else{
            aln += "|";
        }
    }
    printf("\n");
    printf("Alignment:\n");
    printf("X: %s\n", aln1.c_str());
    printf("   %s\n", aln.c_str());
    printf("Y: %s\n", aln2.c_str());
    printf("\n");
    printf("Optimum alignment score: %f\n", matrix(m, n));

}

vector<int> findgapindex(const string& sample1, const string& sample2){
    vector<int> charindex;
    int l = sample1.size();
    for(int i = 0; i < l; i++){
        if (sample1[i] == '-' || sample2[i]=='-')
            charindex.push_back(i);
    }
    return (charindex);
}

vector<int> finddiffindex(const string& sample1, const string& sample2){
    vector<int> diffindex;
    int l = sample1.size();
    for (int i = 0; i < l; i++){
        if ((sample1[i] != '-' && sample2[i] != '-') && sample1[i] != sample2[i])
            diffindex.push_back(i);
    }
    return (diffindex);
}

GlobalAlign Global_Align(string seq1, string seq2, double match, double mismatch, double gap){
    GlobalAlign res;
    string x = "0" + seq1;
    string y = "0" + seq2;
    int q = x.length();
    int r = y.length();
    // initial score matrix
    Eigen::MatrixXd mat(q,r);
    for (int i = 0; i < q; i++){
        mat(i, 0) = gap * i;
    }
    for (int j = 0; j < r; j++){
        mat(0, j) = gap * j;
    }
    // the dynamic programming matrix
    for (int i = 1; i < q; i++){
        for (int j = 1; j < r; j++){
            if (x[i]==y[j]){
                mat(i, j) = mat(i-1, j-1) + match;
            }else{
                double score1 = mat(i-1, j) + gap;
                double score2 = mat(i, j-1) + gap;
                double score3 = mat(i-1, j-1) + mismatch;
                vector<double> scores{score1, score2, score3};
                mat(i, j) = *max_element(scores.begin(), scores.end());
            }
        }
    }
    // Traceback
    string aln1="";
    string aln2="";
    int i = q - 1;
    int j = r - 1;
    bool xz = false;
    bool yz = false;
    
    while(i >= 0 && j >= 0){
        if (i==0 && j == 0){
            break;
        }
        if (i == 0 && j > 0){
            xz = true;
            i = i + 1;
        }
        if (i > 0 &&  j == 0){
            yz = true;
            j = j + 1;
        }
        // case1
        double sc = mat(i-1, j-1);
        if (x[i] == y[j]){
            sc = sc + match;
        }else{
            sc = sc + mismatch;
        }
        if (sc == mat(i,j) && (i == 1 && j > 0 && xz)){
            aln1 = "-" + aln1;
            aln2 = y[j] + aln2;
            i = i - 1;
            j = j - 1;
            continue;
        }
        if (sc == mat(i,j) && (i > 0 &&  j == 1 && yz)){    
            aln1 = x[i] + aln1;
            aln2 = "-" + aln2;
            i = i - 1;
            j = j - 1;
            continue;
        }
        if (sc == mat(i,j) && (!(i == 1 && j > 0 && xz) || !(i > 0 &&  j == 1 && yz))){
            aln1 = x[i] + aln1;
            aln2 = y[j] + aln2;
            i = i - 1;
            j = j - 1;
            continue; 
        }
        // case2
        if ((mat(i-1,j) + gap) == mat(i,j)){
            aln1 = x[i] + aln1;
            aln2 = "-" + aln2;
            i = i - 1;
            continue;
        }
        // case3
        if ((mat(i,j-1) + gap) == mat(i,j)){
            aln1 = "-" + aln1;
            aln2 = y[j] + aln2;
            j = j -1;
            continue;
        }
    }
    res.seq1 = seq1;
    res.seq2 = seq2;
    res.aln1 = aln1;
    res.aln2 = aln2;
    array <int, 3> scs ={0};
    scs[0] = match;
    scs[1] = mismatch;
    scs[2] = gap;
    res.score = scs;
    res.matrix = mat;
    res.gap_position = findgapindex(aln1, aln2);
    res.mis_position = finddiffindex(aln1, aln2);
    return(res);
}

void usage(){
    printf("This script was designed to implement NW algorithm with dynamic programming.\n");
    printf("Author: Shuangbin Xu\n");
    printf("Email: xshuangbin@163.com\n");
    printf("Usage: global_align -q <seq1> -r <seq2> -m [5] -n [-2] -g [-6]\n");
    printf("        -q                      the first nucleotide sequence.\n");
    printf("        -r                      the second nucleotide sequence.\n");
    printf("        -m                      the score of match, default is 5.\n");
    printf("        -n                      the score of mismatch, default is -2.\n");
    printf("        -g                      the score of gap, default is -6.\n");
    printf("        -h                      print the help information.\n");
}

int main(int argc, char *argv[]){
    if (argc < 2){
        usage();
        return (0);
    };
    int opt=0;
    string seq1;
    string seq2;
    double match=5;
    double mismatch=-2;
    double gap=-6;
    while((opt = getopt(argc, argv, "hq:r:m:n:g:")) != -1) {
        switch(opt){
            case 'h':
                usage();
                return(0);
            case 'q':
                seq1 = string(optarg);
                break;
            case 'r':
                seq2 = string(optarg);
                break;
            case 'm':
                match = atoi(optarg);
                break;
            case 'n':
                mismatch = atoi(optarg);
                break;
            case 'g':
                gap = atoi(optarg);
                break;
            case '?':
                printf("unrecognized option\n");
                break;
        }
    }
    GlobalAlign Aln = Global_Align(seq1, seq2, match, mismatch, gap);
    //cout << Aln.matrix << endl;
    Aln.show();
    return(0);
}
