//
//  main.cpp
//  MAX-SurveyPropagation
/*
 Copyright 2018 Raffaele Marino
 This file is part of BSP.
 
 MAX-SP is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 MAX-SP is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with BSP; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

//
//  Created by Raffaele Marino on 2018-08-25.
//  Copyright © 2018 Raffaele Marino. All rights reserved.
//


                                /*******************************************/
                                /*******************************************/
                                /*IL PERDER TEMPO A CHI PIU' SA PIU' SPIACE*/
                                /*******************************************/
                                /*******************************************/

#define VERSION "MAX_E_K_SAT_VERSION_2020"
#include "Header.h"
#include "Vertex.hpp"
#include "Graph.hpp"






/*This code describes message passing algorithms that are useful for solving MAX-E-K-SAT problems.
 This code is an open source code, it can be modified and improved with new strategies.
 For any question about this code (and if you will find bugs), please, contact me at:
 raffaele.marino@epfl.ch or marino.raffaele@mail.huji.ac.il or marinoraffaele.nunziatella@gmail.com */




void help();/*help function*/

int main(int argc,  char * const argv[]) {
    
    ofstream outfile("output.txt", ios_base::app);
    Graph G(argc,argv);/*declaration object Graph*/
    vector<Vertex> V;/*declaration vertex vector*/
    cout<<"I START "<<endl;
    
    /***********************************************************************************/
    /***********************************************************************************/
    /***********************************************************************************/
    /***************************** START WRITE OR LOAD *********************************/
    /***********************************************************************************/
    /***********************************************************************************/
    /***********************************************************************************/
    
    /* At this point the code choses between loading a SAT CNF instance or building one. */
    int c;
    while ((c = getopt (argc, argv, "l:w:h?")) != -1)/* write or load an instance K-SAT*/
        switch (c) {
            case 'l':/* load an instance K-SAT*/
                switch (argc) {
                    case 3:
                        G.read_from_file_graph(); /*read CNF instance from file*/
                        break;
                    default:
                        help(); /*describe how to give INPUT from terminal*/
                        exit(-1);
                        break;
                }
                break;
            case 'w':/* write an instance for SAT problem*/
                G.write_on_file_graph(); /*build and write a CNF instance on file*/
                break;
                
            case 'h':/* help function is called*/
                
            default:
                
                fprintf(stderr, "%s [options]\n"
                        "\n"
                        "-l= reading input file"
                        "\n"
                        "-w= writing input file"
                        "\n"
                        "+++++++++++++++++++++++"
                        "\n"
                        ,argv[0]);
                help();/*help function*/
                exit(-1);
                break;
        }
    /***********************************************************************************/
    /***********************************************************************************/
    /***********************************************************************************/
    /****************************** END WRITE OR LOAD **********************************/
    /***********************************************************************************/
    /***********************************************************************************/
    /***********************************************************************************/
#ifdef NeurNET
    NN NeurNet; /*definition Neural Network object*/
#ifdef NeurNETraining
    cout<<"I AM LEARNING PARAMETERS AT THIS STEP"<<endl;
    cout<<"\t \t \t start training \t \t \t"<<endl;
    NeurNet.Training_parameters();/*Training parameters neural network*/
#endif

#ifdef NeurNETLoad
    cout<<"I AM LOADING PARAMETERS FROM model.xml FILE"<<endl;
    NeurNet.load();/*load parametes neural networks before trained*/
#endif
    
#endif
    /*At this point the code initializes the vertex Vector. */
    V.resize(G.N());/*Initialization Vector V*/
    for (unsigned int i=0; i<G.N(); ++i) {
        V[i]._vertex=i+1; /*label each variable node with a number form 1 to N*/
        V[i]._vertex_lli=(long int)(i+1);
#ifdef NeurNET
        V[i].Get_NN_ptr(NeurNet);
#endif
    }
    
    /*In graph G the vector V is stored in a vector of pointers*/
    G.get_V(V);/*pass vector V to object G*/
    
    /*An instance of MAX-K-SAT, or NAESAT, is given as a conjunction of one or more clauses,
     where a clause is a disjunction of literals.
     Each clause corresponds to a function node, each variable to a variable node,
     and an edge connects a function node and a variable node if and only if the clause contains the variable.
     All these information are collected in split_and_collect_information(). This public member of class Graph
     helps to split a CNF instance into Vertex objects and Clause objects.
     These objects store all types of information about the factor graph.
     The reader can see the specific for each class member in Graph.hpp, Vertex.hpp and Graph.cpp,
     Vertex.cpp files.*/
    
    G.split_and_collect_information();/*split and collect information into graph G*/

    
    cout<<"START SP:"<<endl;
    
    goto SP;
    /***********************************************************************************/
    /***********************************************************************************/
    /***********************************************************************************/
    /********************************* START SP ***************************************/
    /***********************************************************************************/
    /***********************************************************************************/
    /***********************************************************************************/
    
    
    /*The heart of message passing procedures are iterative equations.
     Survey propagation (SP) algorithm is based on equations  derived by cavity method.
     
     For cavity method and analytic formulation of SP equations we refer to :
     
     [1] Marino Raffaele, Giorgio Parisi, and Federico Ricci-Tersenghi. "The backtracking survey propagation algorithm for solving random K-SAT problems." Nature communications 7 (2016): 12996.
     [2] Mézard Marc, Giorgio Parisi, and Riccardo Zecchina. "Analytic and algorithmic solution of random satisfiability problems." Science 297.5582 (2002): 812-815.
     [3] Braunstein, Alfredo, and Riccardo Zecchina. "Survey propagation as local equilibrium equations." Journal of Statistical Mechanics: Theory and Experiment 2004.06 (2004): P06007.
     [4] Mézard Marc, and Giorgio Parisi. "The cavity method at zero temperature." Journal of Statistical Physics 111.1-2 (2003): 1-34.
     [5] Mézard Marc, and Riccardo Zecchina. "Random k-satisfiability problem: From an analytic solution to an efficient algorithm." Physical Review E 66.5 (2002): 056126.
     [6] Parisi Giorgio. "On the survey-propagation equations for the random K-satisfiability problem." arXiv preprint cs/0212009 (2002).
     [7] Parisi Giorgio. "A backtracking survey propagation algorithm for K-satisfiability." arXiv preprint cond-mat/0308510 (2003).
     [8] Aurell Erik, Uri Gordon, and Scott Kirkpatrick. "Comparing beliefs, surveys, and random walks." Advances in Neural Information Processing Systems (2005).
     
     
     SP equations return messages that go from each caluse to each variable node,
     and compute for each variable node surveys. A message sent from a clause c to a variable i
     informs the variable what is the best choice to do for satisfying clause c. This message is computed
     from the messages recived by remaining variables j into clause c, but distinct from i.
     All messages are computed iteratively till a fixed point shows up.
     In this code this procedure is obtained by the public member of class Graph convergence_messages().
     For details, we refer to technical specifications into class Graph.
     Once a convergece is found, surveys for each variable node are computed. They are three,
     and each of them describes the probability that a variable node can be assigned to
     a true value, s_T, to a false value, s_F, or can be indeterminate s_I.
     This calculation, in this code, is performed by the public member surveys().
     For details, we refer to technical specifications in class Graph and class Vertex.*/

SP:
    
    G.convergence_messages();/*find messages convergence*/
    G.surveys();/*compute surveys for variable nodes*/
    G.fix_all_var(); /*fix all variables*/
    cout<<G.N()<<" "<<G.M()<<" "<<G.alpha()<<" "<<G._comp_init<<" "<<G.time_conv()<<" "<<G.check_instance()<<" "<<G.seed()<<" "<<(double)G.check_instance()/(double)G.M()<<" "<<G.unit_prop_print()<<" "<<G._counter_conv<<" "<<G._n_cl_non_conv<<" "<<G.err_medio<<endl;/*check the number of calsuses satisfied*/
    outfile<<G.N()<<" "<<G.M()<<" "<<G.alpha()<<" "<<G._comp_init<<" "<<G.time_conv()<<" "<<G.check_instance()<<" "<<G.seed()<<" "<<(double)G.check_instance()/(double)G.M()<<" "<<G.unit_prop_print()<<" "<<G._counter_conv<<" "<<G._n_cl_non_conv<<" "<<G.err_medio<<endl;/*check the number of calsuses satisfied*/
    return 0;
}


void help(){    /*help function*/
    
    cout << "************** Help Function **************"<<endl;
    
    cout << "\n"<<endl;
    
    cout << "************************************"<<endl;
    
    cout << "\n"<<endl;
    
    cout << "If you use -w: [K] [alpha] [N variables]"<<endl;
    
    cout << "\n"<<endl;
    
    cout << "If you use -l: [formula_INPUT_File_Name.cnf]"<<endl;
    
}
