//
//  NN.hpp
//  Proj_print_survey
//
//  Created by Raffaele Marino on 26.03.20.
//  Copyright © 2020 Raffaele Marino. All rights reserved.
//



#ifndef NN_h
#define NN_h
#include <mlpack/core.hpp>
#include <mlpack/methods/ann/layer/layer.hpp>
#include <mlpack/methods/ann/loss_functions/sigmoid_cross_entropy_error.hpp>
#include <mlpack/methods/ann/layer/flexible_relu.hpp>
#include <mlpack/methods/ann/ffn.hpp>
#include <ensmallen.hpp>
using namespace mlpack;
using namespace mlpack::ann;
using namespace arma;
using namespace ens;
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
/***************************** START NN CLASS **********************************/
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/



/*class Neural Network for predicting the values of variables*/
/*The neural network is composed by 4 layers, the INPUT linear layer with sigmoid activation function, 2 hidden linear layers with sigmoid activation function, and one output linear layer with sigmoid activation function. The goal of the neural network is approximate a function of parameters that outputs the probability that a variable node is true */
class NN{
public:
    /*constructor*/
    NN(){
        open_files();
        /*we are building the neural network*/
        model.Add<Linear<> >(4,40);
        model.Add<SigmoidLayer<> >();
        model.Add<Linear<> >(40,40);
        model.Add<SigmoidLayer<> >();
        model.Add<Linear<> >(40,40);
        model.Add<SigmoidLayer<> >();
        model.Add<Linear<> >(40,40);
        model.Add<SigmoidLayer<> >();
        model.Add<Linear<> >(40,1);
        model.Add<SigmoidLayer<> >();
        /*end of the neural network*/
        X_input.zeros(4,1);
    }
    
    void Training_parameters();
    void load();
    bool Predict(double pp, double pm, unsigned n_n, unsigned n_nn);
    
    
private:
    
    /*open files and store data*/
    void open_files();
    const int N_Epoch=1; //epoch for NN
    FFN<SigmoidCrossEntropyError<>, RandomInitialization> model;
    /*definition storing data matrices*/
    mat trainData; /*training data*/
    mat trainLabels; /*training labels*/
    mat testData; /*test data*/
    mat testLabels; /*test labels*/
    mat test_last_Data; /*last test data*/
    mat test_last_Labels; /*last test labels*/
    mat Y_output; /*return the predicted value of the NN*/
    mat X_input; /*INPUT vector for prediction*/
    
    
};
/*public inline member for outputting the value of the neural network*/
inline bool NN::Predict(double pp, double pm, unsigned n_n, unsigned n_nn){
    X_input(0,0)=pp;
    X_input(1,0)=pm;
    X_input(2,0)=(double) n_n;
    X_input(3,0)=(double) n_nn;
    model.Predict(X_input, Y_output);
    if(Y_output(0,0)>0.5) return true;
        else return false;
}

/*public inline member for load parameters of the neural network*/
inline void NN::load(){
    data::Load("/scratch/rmarino/slurm/model.xml", "model", model);
}

/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
/******************************  END NN CLASS  **********************************/
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/




#endif /* NN_h */
