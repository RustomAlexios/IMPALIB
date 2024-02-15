#include <limits>
#include <string>
#include <bits/stdc++.h>
#include <typeinfo>
#include <chrono>
#include <iostream>
#include <vector>
#include<fstream>
#include <algorithm>
#include <numeric>
#include <functional>
#include "../../../build/cnpylib/include/cnpy.h"
using namespace std;

#define IMPALIB_TYPE double

#ifdef IMPALIB_TYPE
    typedef IMPALIB_TYPE impalib_type;
#else
    typedef float impalib_type;
#endif

const int N_BINS = 100;
const int N_ITEMS = 10000;

impalib_type zero_value = 0;

//export LD_LIBRARY_PATH=../../../build/cnpylib
//g++ -o test_EOP test_EOP.cpp -L../../../build/cnpylib -lcnpy -lz --std=c++11

void extrinsic_output_package_update(vector<vector<impalib_type>> , vector<impalib_type> );

void extrinsic_output_package_update( vector<vector<impalib_type>> rExtrinsicOutputBin,
                                    vector<impalib_type> rOric2PackageM){
        
        vector<impalib_type> ExtrinsicOutputPackage(N_ITEMS, zero_value);

        copy(rOric2PackageM.begin(), rOric2PackageM.end(), ExtrinsicOutputPackage.begin());
        
        
        
        for (int bin_index = 0; bin_index < rExtrinsicOutputBin.size(); bin_index++){
            //vector<impalib_type> intermediate_vector(N_ITEMS, zero_value);
            
            //transform(rExtrinsicOutputBin[bin_index].begin(), rExtrinsicOutputBin[bin_index].end(), 
            //                        ExtrinsicOutputPackage.begin(), intermediate_vector.begin(), std::plus<impalib_type>());
            
            //copy(intermediate_vector.begin(), intermediate_vector.end(), ExtrinsicOutputPackage.begin());
            //vector<impalib_type> intermediateVector(rExtrinsicOutputBin[bin_index]);
            transform(rExtrinsicOutputBin[bin_index].begin(), rExtrinsicOutputBin[bin_index].end(), 
                                    ExtrinsicOutputPackage.begin(), ExtrinsicOutputPackage.begin(), std::plus<impalib_type>());
            }

        vector<impalib_type> ExtrinsicOutputPackageDummy(N_ITEMS, 0.0);

        for (int bin_index=0; bin_index<N_BINS; bin_index++){
            for (int item_index=0; item_index<N_ITEMS; item_index++){
                ExtrinsicOutputPackageDummy[item_index] = ExtrinsicOutputPackageDummy[item_index] + rExtrinsicOutputBin[bin_index][item_index];
            }
        }
        
        for (int item_index=0; item_index<N_ITEMS; item_index++){
                ExtrinsicOutputPackageDummy[item_index] = ExtrinsicOutputPackageDummy[item_index] + rOric2PackageM[item_index];
                if (ExtrinsicOutputPackageDummy[item_index] != ExtrinsicOutputPackage[item_index]){
                    cout<<"Failed"<<endl;
                    cout<<ExtrinsicOutputPackageDummy[item_index]<<endl;
                    cout<<ExtrinsicOutputPackage[item_index]<<endl;
                    exit(7);
                }
            }
        
        cout<<"Succeeded"<<endl;
        
        /*fstream file_output("extrinsic_output_package_wrapper", ios::out | ios::binary | ios:: trunc);
                if (file_output.is_open()) {
                    for (int i=0; i<N_ITEMS; i++){
                        file_output.write((char*)(&ExtrinsicOutputPackage[i]), sizeof(ExtrinsicOutputPackage[i]));}
                        file_output.close();}
                else {cout << "Error! File cannot be opened!" << endl;} */
}

int main(){
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> doubleDist(10, 300);

        vector<impalib_type> oric_to_package_m(N_ITEMS,zero_value);
        vector<vector<impalib_type>> extrinsic_output_bin(N_BINS, vector<impalib_type>(N_ITEMS,zero_value));

        for (int i=0; i<N_ITEMS; i++){
            oric_to_package_m[i] = doubleDist(gen);
            for (int j=0; j<N_BINS; j++){
                extrinsic_output_bin[j][i] = doubleDist(gen);
            }   
        }

        /*cnpy::NpyArray input1 = cnpy::npy_load("oric_to_package_m.npy");
        impalib_type* oric_to_package_m_pure = input1.data<impalib_type>();
        vector<impalib_type> oric_to_package_m(N_ITEMS,zero_value);
        copy(oric_to_package_m_pure, oric_to_package_m_pure + N_ITEMS, oric_to_package_m.begin());

        cnpy::NpyArray input2 = cnpy::npy_load("extrinsic_output_bin.npy");
        impalib_type* extrinsic_output_bin_pure = input2.data<impalib_type>();
        

        for (int bin_index=0; bin_index<N_BINS; bin_index++){
        copy ( extrinsic_output_bin_pure + N_ITEMS*bin_index, extrinsic_output_bin_pure+N_ITEMS*(bin_index+1), extrinsic_output_bin[bin_index].begin() );
        }*/

        extrinsic_output_package_update(extrinsic_output_bin, oric_to_package_m);
}