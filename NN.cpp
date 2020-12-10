//
//  NN.cpp
//  Proj_print_survey
//
//  Created by Raffaele Marino on 26.03.20.
//  Copyright © 2020 Raffaele Marino. All rights reserved.
//

#include "NN.hpp"

void NN::open_files(){
    mat dataset_training_surveys;
    mat dataset_training_solutions;
    mat dataset_test_surveys;
    mat dataset_test_solutions;
    mat dataset_last_test_surveys;
    mat dataset_last_test_solutions;
    data::Load("./training/SurveysK=3N=10000alpha=4.2_training.txt", dataset_training_surveys, true);
    data::Load("./training/SurveysK=3N=10000alpha=4.23_test.txt", dataset_test_surveys, true);
    data::Load("./training/SurveysK=3N=10000alpha=4.24_test.txt", dataset_last_test_surveys, true);
    data::Load("./training/Sol_SP_for_nnK=3N=10000alpha=4.23_test.txt", dataset_test_solutions, true);
    data::Load("./training/Sol_SP_for_nnK=3N=10000alpha=4.2_training.txt", dataset_training_solutions, true);
    data::Load("./training/Sol_SP_for_nnK=3N=10000alpha=4.24_test.txt", dataset_last_test_solutions, true);
    trainData = dataset_training_surveys.submat(0, 0, dataset_training_surveys.n_rows - 3, dataset_training_surveys.n_cols - 1);
    trainLabels = dataset_training_solutions.submat(0, 0, dataset_training_solutions.n_rows - 3, dataset_training_solutions.n_cols - 1);
    testData = dataset_test_surveys.submat(0, 0, dataset_test_surveys.n_rows - 3, dataset_test_surveys.n_cols - 1);
    testLabels = dataset_test_solutions.submat(0, 0, dataset_test_solutions.n_rows - 3, dataset_test_solutions.n_cols - 1);
    test_last_Data = dataset_last_test_surveys.submat(0, 0, dataset_last_test_surveys.n_rows - 3, dataset_last_test_surveys.n_cols - 1);
    test_last_Labels = dataset_last_test_solutions.submat(0, 0, dataset_last_test_solutions.n_rows - 3, dataset_last_test_solutions.n_cols - 1);

    /*for (uword i=0; i<trainData.n_rows; i=i+5) {
        for (uword j=0; j<trainData.n_cols; j++) {
            trainData(0+i, j)=0.;
            trainData(1+i, j)=0.;
            trainData(2+i, j)=0.;
            
        }
    }

    for (uword i=0; i<testData.n_rows; i=i+5) {
        for (uword j=0; j<testData.n_cols; j++) {
            testData(0+i, j)=0.;
            testData(1+i, j)=0.;
            testData(2+i, j)=0.;
        }
    }
    
    for (uword i=0; i<test_last_Data.n_rows; i=i+5) {
        for (uword j=0; j<test_last_Data.n_cols; j++) {
            test_last_Data(0+i, j)=0.;
            test_last_Data(1+i, j)=0.;
            test_last_Data(2+i, j)=0.;
        }
    }
    
    */
    for (uword i=0; i<trainLabels.n_rows; ++i) {
        for (uword j=0; j<trainLabels.n_cols; ++j) {
            if(trainLabels(i,j)>0)trainLabels(i,j)=1;
            else trainLabels(i,j)=0;
        }
    }
    for (uword i=0; i<testLabels.n_rows; ++i) {
        for (uword j=0; j<testLabels.n_cols; ++j) {
            if(testLabels(i,j)>0)testLabels(i,j)=1;
            else testLabels(i,j)=0;
        }
    }
    for (uword i=0; i<test_last_Labels.n_rows; ++i) {
        for (uword j=0; j<test_last_Labels.n_cols; ++j) {
            if(test_last_Labels(i,j)>0)test_last_Labels(i,j)=1;
            else test_last_Labels(i,j)=0;
        }
    }
}

/*public member for training parameters of neural network*/
void NN::Training_parameters(){
    mat X_batches; /*single data*/
    mat Y_batches; /*single label*/
    AdaDelta opt(1.0, 1, 0.99, 1e-8, 1000, 1e-9, true);
#ifdef TRAIN_ANAL
    mat X_test_err; /*single label*/
    mat Y_test_err; /*single label*/
    mat Z_test_err; /*single label*/
    vector<double> accurancy;
    vector<double> sdvaccurancy;
    accurancy.reserve(testLabels.n_rows);
    sdvaccurancy.reserve(testLabels.n_rows);
#endif
    for (int epoch = 0; epoch < N_Epoch; ++epoch)
    {
        cout<<"Epoch Training: "<<epoch<<endl;
        for (int b=0; b<trainLabels.n_rows; ++b) {
            X_batches = trainData.submat(b*4, 0, (b*4)+3, trainLabels.n_cols-1);
            Y_batches = trainLabels.submat(b, 0, b, trainLabels.n_cols-1);
            model.Train(X_batches, Y_batches, opt);
#ifdef TRAIN_ANAL
            double sum=0.;
            vector<double> stdv;
            stdv.resize(testLabels.n_cols);
            for (int d=0; d<testLabels.n_cols; ++d) {
                unsigned counter=0;
                  for(int c=0; c<testLabels.n_rows; ++c){
                    X_test_err =testData.submat(c*4, d, (c*4)+3, d);
                    Z_test_err=testLabels.submat(c, d, c, d);
                    model.Predict(X_test_err, Y_test_err);
                    if(Y_test_err(0,0)>0.5)Y_test_err(0,0)=1.;
                    else Y_test_err(0,0)=0.;
                    if(Y_test_err(0,0)==Z_test_err(0,0))counter++;
                }
                sum+=(double)counter/(double)(testLabels.n_rows);
                stdv[d]=(double)counter/(double)(testLabels.n_rows);
                
            }
            accurancy.push_back(sum/(double)(testLabels.n_cols));
            double sum_sdv=0.;
            for (unsigned i=0; i<stdv.size(); ++i) {
                sum_sdv+=(stdv[i]-(sum/(double)(testLabels.n_cols)))*(stdv[i]-(sum/(double)(testLabels.n_cols)));
            }
            sdvaccurancy.push_back(sqrt(sum_sdv/((double)stdv.size()-1.)));
            cout<<b<<" "<<sum/(double)(testLabels.n_cols)<<" "<<sqrt(sum_sdv/((double)stdv.size()-1.))<<endl;
#endif
        }
    }
    cout<<"I finished the trainging, I test my learning values"<<endl;
    mat assignments;
    mat true_val, pred_val, sigmoidal_val;
    sigmoidal_val.zeros(1, testLabels.n_rows);
    pred_val.zeros(1, testLabels.n_rows);
    true_val.zeros(1, testLabels.n_rows);
    for (int epoch = 0; epoch < testLabels.n_cols; ++epoch)
    {
        cout<<"Epoch Testing: "<<epoch<<endl;
        
        for (int b=0; b<testLabels.n_rows; ++b) {
            X_batches = testData.submat(b*4, epoch, (b*4)+3, epoch);
            Y_batches = testLabels.submat(b, epoch, b, epoch);
            
            model.Predict(X_batches, assignments);
            sigmoidal_val(0,b)=assignments(0,0);
            if(assignments(0,0)>0.5)assignments(0,0)=1.;
            else assignments(0,0)=0.;
            pred_val(0,b)=assignments(0,0);
            true_val(0,b)=Y_batches(0,0);
        }
        unsigned counter=0;
        ofstream out_distrib("distribution.txt");
        for (unsigned int i=0; i<testLabels.n_rows; ++i) {
            if(pred_val(0,i)==true_val(0,i)){
            counter++;
                out_distrib<<"V "<<pred_val(0,i)<<" "<<true_val(0,i)<<" "<<sigmoidal_val(0,i)<<endl;
            }else{
                out_distrib<<"F "<<pred_val(0,i)<<" "<<true_val(0,i)<<" "<<sigmoidal_val(0,i)<<endl;
            }
        }
        cout<<"Fraction of true labels obtained by NN: "<<(double)counter/(double)testLabels.n_rows<<endl;
        
        
    }
    
    for (int epoch = 0; epoch < test_last_Labels.n_cols; ++epoch)
    {
        cout<<"Epoch LAST Testing: "<<epoch<<endl;
        
        for (int b=0; b<test_last_Labels.n_rows; ++b) {
            X_batches = test_last_Data.submat(b*4, epoch, (b*4)+3, epoch);
            Y_batches = test_last_Labels.submat(b, epoch, b, epoch);
            
            model.Predict(X_batches, assignments);
            if(assignments(0,0)>0.5)assignments(0,0)=1.;
            else assignments(0,0)=0.;
            pred_val(0,b)=assignments(0,0);
            true_val(0,b)=Y_batches(0,0);
        }
        unsigned counter=0;
        for (unsigned int i=0; i<testLabels.n_rows; ++i) {
            if(pred_val(0,i)==true_val(0,i))counter++;
        }
        cout<<"Fraction of true labels obtained by NN: "<<(double)counter/(double)testLabels.n_rows<<endl;
        
    }
    data::Save("model.xml", "model", model, false);
#ifdef TRAIN_ANAL
    ofstream out("Test_error.txt");
    for (unsigned i=0; i<accurancy.size(); ++i) {
        out<<i<<" "<<accurancy[i]<<" "<<sdvaccurancy[i]<<endl;
    }
#endif
}
