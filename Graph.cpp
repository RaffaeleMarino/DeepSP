//
//  Graph.cpp
//  TheBacktrackingSurveyPropagation
/*
 Copyright 2018 Raffaele Marino
 This file is part of BSP.
 
 BSP is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 BSP is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with BSP; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

//
//  Created by Raffaele Marino on 26/11/2018.
//  Copyright © 2018 Raffaele Marino. All rights reserved.
//

#include "Graph.hpp"

/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
/***************************** START CLAUSE CLASS **********************************/
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/


/*public member class Clause. This member helps to clean the graph from satisfied clauses and
 from unsitisfied literal in unsitisified clauses. It starts to store the value of the message from the clause j
 to the variable i, the one which has been fixed previously. The member store the survey in survey_cl_to_i_f
 and set in variable i the survey _surveys_cl_to_i[j] equl to 1. Then check if the clause is satisfied or not*/

void Clause::clean(Vertex * &p_V, unsigned int &j){/*clean the clasue*/
    bool flag_deg=true;
    survey_cl_to_i_f[*imag(p_V->_I_am_in_cl_at_init[j])]=p_V->_surveys_cl_to_i[j];/*store clause to variable survey in the clause*/
    if(cp_lit[*imag(p_V->_I_am_in_cl_at_init[j])]){
        (p_V->_who_I_am)?p_V->_surveys_cl_to_i[j]=0.:p_V->_surveys_cl_to_i[j]=1.;/*set survey clause to variable  into a Vertex object to one*/
        _vb[*imag(p_V->_I_am_in_cl_at_init[j])]=p_V->_who_I_am;/*update bool vector  _vb*/
    }else{
        (!p_V->_who_I_am)?p_V->_surveys_cl_to_i[j]=0.:p_V->_surveys_cl_to_i[j]=1.;/*set survey clause to variable  into a Vertex object to one*/
        _vb[*imag(p_V->_I_am_in_cl_at_init[j])]=(!(p_V->_who_I_am));/*update bool vector  _vb*/
    }
    _v_ref[*imag(p_V->_I_am_in_cl_at_init[j])]=_cast_b_to_i(_vb[*imag(p_V->_I_am_in_cl_at_init[j])]);
    _go_forward[*imag(p_V->_I_am_in_cl_at_init[j])]=0;
    _lit.erase(p_V->_lit_list_i[j]);/*erase literal from the clause*/
    V.erase(p_V->where_I_am[j]);/*erase variable form the clause*/
    _it_list_V_begin=V.begin();
    _size_cl--;
    survey_cl_to_i.erase(p_V->_where_surveys_are_in_cl[j]);/*erase the survey clause to i from the clause*/
    survey_begin=survey_cl_to_i.begin();
    if(I_am_a_cl_true){
        flag_deg=false;
    }
    I_am_a_cl_true=_logic_operator();/*check if the clause has been satisfied or not*/
    if(flag_deg and I_am_a_cl_true){
        /*if the clause has been satisfied then save surveys, set the other variables surveys to 0  and */
        for (unsigned int i=0; i<survey_cl_to_i_f.size(); ++i) {
            if(survey_cl_to_i_f[i]==-1.){
                survey_cl_to_i_f[i]=*survey_cl_to_i_f_ptr[i];
                *survey_cl_to_i_f_ptr[i]=0.;
            }
        }
        /*update the variable nodes degree of the other variables that have not been fixed yet.*/
        for (list<Vertex *>::iterator __it=V.begin(); __it!=V.end(); ++__it) {
            (*__it)->_degree_i--;
            if((*__it)->_degree_i<0){
                cout<<"Error degree node j less than 0"<<endl;
                cout<<"Il grado e':"<<(*__it)->_degree_i<<" "<<(*__it)->_vertex<<endl;
                exit(-1);
            }
        }
        
    }
}

/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
/******************************  END CLAUSE CLASS **********************************/
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/



/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
/*****************************  START GRAPH CLASS **********************************/
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/


/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
/********************************** START SP ***************************************/
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/




void Graph::clean(Vertex * &V_to_clean){
    unsigned int k=0;
    unsigned int temp;
    for (unsigned int j=0; j<V_to_clean->_surveys_cl_to_i.size(); ++j) {
        k=*real(V_to_clean->_I_am_in_cl_at_init[j]);
        _cl[k].clean(V_to_clean, j);/*clean the graph*/
        if(_cl[k].I_am_a_cl_true and _cl[k]._I_am_in_list_unsat){
            /*At this point the algorithm  picks a clause sat and set it on the top of the vector vec_list*/
            /*The vector vec_list_cl contains pointers to clauses sat and unsat. The value _m divides the sat clause, from 0 to _m-1, to unsat clause,
             from _m to  _M*/
            
            temp=_cl[k]._l;
            swap(vec_list_cl[_cl[k]._l], vec_list_cl[_m]);
            (vec_list_cl[_cl[k]._l])->_l=temp;
            (vec_list_cl[_m])->_l=_m;
            _cl[k]._I_am_in_list_unsat=false;
            _cl_list.erase(_cl[k]._it_list);
            _m++;
        }
    }
}



/*public memeber class Graph. This memeber computes the surveys for each variable node.
 It is called when a convergence of all messages from a clause to variable is found.
 Moreover this member computes the variable complexity associated to each variable node*/
void Graph::surveys(){/*compute surveys for each variable node*/
    complexity_variables=0.;/*set complexity variable to 0*/
    unsigned int _size_init=static_cast<unsigned int>(_list_fixed_element.size());
    for (unsigned int i=_size_init; i<_N; ++i) {
        ptrV[i]->compute_s();/*compute surveys for each variable node*/
        complexity_variables+=(static_cast<double>(ptrV[i]->_degree_i)-1.)*log(ptrV[i]->complexity_variable);/*update graph variable complexity*/
    }
    complexity=complexity_clauses-complexity_variables;/*compute graph total complexity*/
    if(_numb_of_dec_moves==1 and _numb_of_back_moves==1) _comp_init=complexity;
}

/*public member class Graph which describe unit propagation algorithm. This member is called when a clause is
 composed only by one literal, and therefore has to be satisfied, otherwise we obtain a contraddiction*/
void Graph::unit_propagation(unsigned int  &c, unsigned long &i){ /*fix the variable into clause c to true and clean the graph*/
    //cout<<"unit propagation"<<endl;
    _unit_prop++;/*update counter unit propagation*/
     Vertex * cp=(_cl[c].v_V[i]); /*variable node copy*/
    if (_cl[c].v_lit[i]){
        cp->_who_I_am=true; /*set the variable to the value that satisfied the literal into the clause*/
        cp->_sT=1.;/*set the associated survey to one*/
        cp->_sF=0.;
        cp->_sI=0.;
    }else {
        cp->_who_I_am=false;/*set the variable to the value that satisfied the literal into the clause*/
        cp->_sT=0.;
        cp->_sF=1.;/*set the associated survey to one*/
        cp->_sI=0.;
    }
    
   
    cp->_sC=1.;/*set the certitude to 1*/
    cp->_degree_i=0;/*set variable node degree to 0*/
    _list_fixed_element.push_back(cp);/*store the Vertex in the listt of fixed element*/
    cp->_it_list_fixed_elem=--_list_fixed_element.end();/*save its address into the varibale node*/
    cp->_I_am_a_fixed_variable=true;/*fix the variable*/
    /*clean the graph*/
    /*
     Clean the graph means:
     
     1) erasing from the factor graph all satisfied clauses;
     2) erasing literals associated to variable i, which are present in a clause that is not satisfied by the variable node assignement.
     */
    clean(cp);
    _N_t=_N-static_cast<unsigned int>(_list_fixed_element.size());/*update the number of un-fixed variable nodes*/
}

/*public member class Graph. This member update all messages from clauses to variables and stops if:
 a) a convergence is found : SUCCESS;
 b) a contradiction is found : exit FAILURE;
 c) no convergence is found after t_max iteration: exit FAILURE
 */
void Graph::convergence_messages(){/*compute convergence messages for message passing algorithm*/
    bool conv_f=false;
    unsigned long i,l, __k;
    unsigned int C;
    update_products();
    /*convergence*/
    for (unsigned int t=0; t<t_max; ++t) {
        conv_f=false;
        _counter_conv=0;/*set counter convergence to zero*/
        complexity_clauses=0.;/*set clauses complexity to zero*/
        i=0;
        l=0;
        for(__k=_m; __k<_M; ++__k){/*for on clause objects*/
            C=(vec_list_cl[__k])->_c;
            s=_cl[C]._size_cl_init;
            i=0;
            while (1) {
                if (_cl[C]._go_forward[i]) {
                    _new=1.;
                    norm=1.;
                    l=0;
                    while (1) {
                        if(i!=l && _cl[C]._go_forward[l]){
                            _prod_S=_Pr_S(_cl[C].v_V[l],_cl[C].v_lit[l], _cl[C].div_s[l]);
                            _prod_U=_Pr_U(_cl[C].v_V[l],_cl[C].v_lit[l]);
                            _new*=__pu();/*new message from cl to variable is computed*/
                            norm*=__norm();
                        }
                        ++l;
                        if(l==s)break;
                    }
                    _new=compute_message(_new,norm);/*set to 0 a message iff the message is smaller than 1e-16*/
                    _cl[C].old_s[i]=_cl[C].update[i];
                    if(_new!=1.){
                    _cl[C].update[i]=_new;/*update new message value*/
                    }else{
                    _cl[C].update[i]=1.;/*update new message value*/
                    unit_propagation(C, i); /*fix the variable into clause c to true and clean the graph*/
                        t=0;
                        break;
                    }
                    if(_counter_conv==0){
                        if(_conv(_cl[C].update[i], _cl[C].old_s[i])){
                            ++_counter_conv;
                        }
                    }
                    ++i;
                }else{
                    ++i;
                }
                if(i==s)break;
            }
            i=0;
            while (1) {
                *_cl[C].v_survey_cl_to_i[i]=_cl[C].update[i];
                if(_cl[C].update[i]!=1.)_cl[C].div_s[i]=Div_s(_cl[C].update[i]);
                ++i;
                if(i==s)break;
            }
        }
        update_products();/*update products into vertex node for speeding up the algorithm*/
        if(_counter_conv==0){/*if counter convergence is zero, a convergence is found*/
            _M_t=static_cast<unsigned int>(_cl_list.size());
           // cout<<"I found a convergence at "<<t<<" "<<_m<<endl;
            _time_conv_print=t;
            update_complexity_clauses();/*update complexity_clause*/
            conv_f=true;
            break;
        }
    }
    if(!conv_f){/*if after t_max iterations no convergence is found, the algorithm return exit failure */
        // cout<<"I found a convergence at "<<t<<" "<<_m<<endl;
        
        _time_conv_print=t_max;
        number_conv_surv();
        update_complexity_clauses();/*update complexity_clause*/
        
    }
}

/*public member class Graph. This member  updates only clause complexity.*/
void Graph::update_complexity_clauses(){
    double _ps=1.,_pu=1.;/*products for clause complexity*/
    unsigned long j;
    unsigned int C;
    //cout<<"This is the value where I start: "<<_M-_m<<endl;
    for(unsigned int __k=_m; __k<_M; ++__k) {
        C=(vec_list_cl[__k])->_c;
        _ps=1.;/*initialize to one product __ps*/
        _pu=1.;/*initialize to one product __pu*/
        /*updating complexity clauses*/
        for (j=0; j<_cl[C]._size_cl_init; ++j) { /*for each literal in a clause clauses complexity is updated */
            if(_cl[C]._go_forward[j]==1){

                _prod_S=_Pr_S(_cl[C].v_V[j],_cl[C].v_lit[j], _cl[C].div_s[j]);
                _prod_U=_Pr_U(_cl[C].v_V[j],_cl[C].v_lit[j]);
                _ps*=__norm();/*update product for clause complexity*/
                _pu*=__pu();/*update product for clause complexity*/
            }
        }
        complexity_clauses+=log(_ps-_pu);/*update clause complexity*/
    }
    cout<<complexity_clauses<<endl;
}

/*public member class Graph which splits the global factor graph information in different objects and vectors.*/
void Graph::split_and_collect_information(){/*split the graph in different vectors and lists*/
    /*create clauses*/
    vec_list_cl.resize(_M);
    list<Vertex *>::iterator _it; /*iterator list of Vertex pointers*/
    double _rn=0.; /*random number*/
    unsigned int pos=0; /*position into a clause*/
    _cl[0].I_am_a_cl_true=false;/*set to false bool variable I_am_a_cl_true*/
    _cl_list.push_front(&_cl[0]);/*insert pointer of clause into the unsitisfied clauses list*/
    _cl[0]._it_list=_cl_list.begin();/*store the iterator of the unsitisfied clauses list*/
    _cl[0]._I_am_in_list_unsat=true;/*set to true the boolean varibale I_am_in_list_unsat*/
    vec_list_cl[0]=&_cl[0]; /*pointer to an element of a vector of Clause objects*/
    _m=0;
    for (unsigned int i=0, l=0; i<_ivec.size(); ++i) {
        if(_ivec[i]==0){
            pos=0; /*set position to 0*/
            ++l; /*clause label increment*/
            if(l==_M)break;
            _cl[l].I_am_a_cl_true=false;/*set to false bool variable I_am_a_cl_true*/
            _cl_list.push_front(&_cl[l]);/*insert pointer of clause into the unsitisfied clauses list*/
            _cl[l]._it_list=_cl_list.begin();/*store the iterator of the unsitisfied clauses list*/
            _cl[l]._I_am_in_list_unsat=true;/*set to true the boolean varibale I_am_in_list_unsat*/
            vec_list_cl[l]=&_cl[l];
        }else{
            _cl[l]._c=l;/*clause labeling*/
            _cl[l]._l=l;/*position in vec_listy*/
            _cl[l].survey_cl_to_i_f_ptr.push_back(NULL); /*initialization vector of frozen survey pointers*/
            _cl[l].survey_cl_to_i_f.push_back(-1.);/*initialization vector of frozen surveys*/
            _cl[l]._vb.push_back(false); /*initialization Boolean vector*/
            _cl[l]._v_ref.push_back(-1);
            _cl[l].V.push_back(ptrV[abs(_ivec[i])-1]);/*Vertex pointer stored in list V*/
            _cl[l]._it_list_V_begin=_cl[l].V.begin();/*iterator stored into the clause*/
            _cl[l]._size_cl++; /*size of the clause*/
            _cl[l]._size_cl_init++; /*size of the clause at the beginning*/
            _cl[l].cpV.push_back(ptrV[abs(_ivec[i])-1]);/*copy of V in a vertex vector*/
            _cl[l]._vecpos.push_back(pos);/*vector of literal position in the clause*/
            _it=--_cl[l].V.end();/*pick iterator of V in cl*/
            ptrV[abs(_ivec[i])-1]->where_I_am.push_back(_it);/*store iterator*/
            ptrV[abs(_ivec[i])-1]->_degree_i++;/*update degree vertex i*/
            /*initialization surveys, real t, imag t-1*/
            _rn=random()/(RAND_MAX+1.0);/*choose a random number between 0-1*/
            ptrV[abs(_ivec[i])-1]->_surveys_cl_to_i.push_back(_rn);/*survey random initialization*/
            _cl[l].old_s.push_back(_rn);
            _cl[l].update.push_back(_rn);
            _cl[l].div_s.push_back(Div_s(_rn));
            pos++;/*update position*/
            if(_ivec[i]>0){/*check if literal is negated or not*/
                _cl[l]._lit.push_back(true);/*store literal in Boolean list _lit*/
                _cl[l].cp_lit.push_back(true);/*store literal in a boolean vector*/
                _cl[l].cp_lit_int.push_back(1);/*store literal as integer*/
            }else{
                _cl[l]._lit.push_back(false);/*store literal in Boolean list _lit*/
                _cl[l].cp_lit.push_back(false);/*store literal in a boolean vector*/
                _cl[l].cp_lit_int.push_back(0);/*store literal as integer*/
                
            }
            ptrV[abs(_ivec[i])-1]->_lit_list_i.push_back(--_cl[l]._lit.end());/*store iterator in Vertex i*/
        }
    }
    /*To readers: please do not chenge this part, because we are using dynamic vectors and when one needs to store pointers of these vectors one needs to make safty choices. More precisely, each time that vector<T> constructor is called for allocating new memory, pointers can change and problems appear*/
    bool flag_s=true;
    _cl[0].v_survey_cl_to_i.resize(_cl[0].size());
    _cl[0]._go_forward.resize(_cl[0].size());
    _cl[0].v_lit.resize(_cl[0].size());
    _cl[0].v_V.resize(_cl[0].size());
    _cl[0]._var.resize(_cl[0].size());
    for (unsigned int i=0,k=0, l=0; i<_ivec.size(); ++i) {
        if(_ivec[i]==0){
            k=0;
            ++l; /*increment of one clause label*/
            flag_s=true;
            if(l==_M)break;
            _cl[l]._go_forward.resize(_cl[l].size());
            _cl[l].v_survey_cl_to_i.resize(_cl[l].size());
            _cl[l].v_lit.resize(_cl[l].size());
            _cl[l].v_V.resize(_cl[l].size());
            _cl[l]._var.resize(_cl[l].size());
        }else{
        
            _cl[l].v_V[k]=ptrV[abs(_ivec[i])-1];/*Vertex pointer stored in list V*/
            _cl[l]._go_forward[k]=1;
            ptrV[abs(_ivec[i])-1]->_I_am_in_cl_at_init.push_back(complex<unsigned int *>(_cl[l]._ptr_c(),_cl[l]._ptr_pos(k)));/*store initial position of Vertex i into clause l*/
            ptrV[abs(_ivec[i])-1]->_bvec_lit.push_back(_cl[l]._ptr_lit_int(k));/*pointers vector literal variablee node*/
            _cl[l].survey_cl_to_i.push_back(ptrV[abs(_ivec[i])-1]->ptr_survey());/*store pointer of survey in cl l*/
            _cl[l].v_survey_cl_to_i[k]=ptrV[abs(_ivec[i])-1]->ptr_survey();
            _cl[l]._var[k]=_ivec[i];
            if(flag_s)_cl[l].survey_begin=_cl[l].survey_cl_to_i.begin();
            flag_s=false;
            if(_ivec[i]>0){
                _cl[l].v_lit[k]=true;
                ptrV[abs(_ivec[i])-1]->_surveys_cl_to_i_plus.push_back(ptrV[abs(_ivec[i])-1]->ptr_survey());
                
            }else{
                _cl[l].v_lit[k]=false;
                ptrV[abs(_ivec[i])-1]->_surveys_cl_to_i_minus.push_back(ptrV[abs(_ivec[i])-1]->ptr_survey());
            }
            _cl[l].survey_cl_to_i_f_ptr[k++]=ptrV[abs(_ivec[i])-1]->ptr_survey();
            ptrV[abs(_ivec[i])-1]->_where_surveys_are_in_cl.push_back(--_cl[l].survey_cl_to_i.end());
            ptrV[abs(_ivec[i])-1]->_i++;
        }
    }
    
    update_products();/*update products*/
}

void Graph::update_products(){/*update products*/
    unsigned int _size_init=static_cast<unsigned int>(_list_fixed_element.size());
    for (unsigned int i=_size_init; i<_N; ++i) {
        //        if(!ptrV[i]->_I_am_a_fixed_variable){
        ptrV[i]->make_products();/*update products*/
        ptrV[i]->_i=0/*set variable _i in Vertex i to 0*/;
        //        }
    }
}

/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
/************************************ END SP ***************************************/
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/

/*MAX SAT fixing*/
#ifdef MAXSAT
/*public member that fixes all variable only once*/
void Graph::fix_all_var(){
    for (unsigned int i=0; i<ptrV.size(); ++i) {
        if(!ptrV[i]->_I_am_a_fixed_variable){
        _list_fixed_element.push_back(ptrV[i]);/*store the variable node into fixed element list*/
        ptrV[i]->_it_list_fixed_elem=--_list_fixed_element.end();/*save its address into the varibale node*/
        ptrV[i]->_I_am_a_fixed_variable=true;/*fix the variable*/
        ptrV[i]->fix_var_i();/*set the variable to the value predicted by the rule described in main.cpp*/
        if(_last_certitude>ptrV[i]->_sC)_last_certitude=ptrV[i]->_sC;/*store last certitude*/
        
        /*clean the graph*/
        /*
         Clean the graph means:
         
         1) erasing from the factor graph all satisfied clauses;
         2) erasing literals associated to variable i, which are present in clauses not satisfied by the variable node assignement.
         */
        clean(ptrV[i]);
        ptrV[i]->_degree_i=0;/*set to 0 the degree of the variable node*/
        }
    }
    
    _N_t=_N-static_cast<unsigned int>(_list_fixed_element.size());/*update the number of un-fixed variable nodes*/
   // cout<<"How many variables do I fixed? "<< _N-_N_t<<endl;
}
    


/*public member that checks how many sat variables satisfied exist */
double Graph::check_instance(){
    double counter_unsat=0;
    /*check if the global solution is correct*/
    for(unsigned int c=0; c<_M;++c){    /*problem check*/
        if (!_cl[c]._check()) {
            counter_unsat++;
        }
    }
    return counter_unsat;
}

void Graph::number_conv_surv(){
    _counter_conv=0;
    _n_cl_non_conv=0;
    err_medio=0.;
    unsigned int C=0;
    for(unsigned int __k=_m; __k<_M; ++__k){/*for on clause objects*/
        C=(vec_list_cl[__k])->_c;
        bool flag=true;
        for (unsigned int i=0; i<_cl[C]._size_cl_init; ++i) {
            if(_conv(_cl[C].update[i], _cl[C].old_s[i])){
                ++_counter_conv;
                err_medio+=abs(_cl[C].update[i]-_cl[C].old_s[i]);
                flag=false;
            }
        }
        if (!flag)_n_cl_non_conv++;
    }
    err_medio=err_medio/(double)_counter_conv;
    
}

#endif


/*public member class Graph. This member prints the variable surveys at the first iteration*/
void Graph::print_only_one_sol(){
    /*output file*/
    const string title="Sol_SP_for_nn";
    const string txt=".dat";
    string str;
    ostringstream _seeds, _Ks, _sN, _salpha;
    _seeds<<_seed;
    _Ks<<_K;
    _salpha<<_alpha;
    _sN<<_N;
    str=directory;
    str=title;
    str+="K=";
    str+=_Ks.str();
    str+="N=";
    str+=_sN.str();
    str+="alpha=";
    str+=_salpha.str();
    str+=txt;
    ofstream outfilesol_(str.c_str(), ios_base::app);
    vector<long int> mysol;
    mysol.resize(_N+1);
    for (list<Vertex*>::iterator __it=_list_fixed_element.begin(); __it!=_list_fixed_element.end(); ++__it) {
        unsigned int var=(*__it)->_vertex;
        if((*__it)->_who_I_am)mysol[(*__it)->_vertex]=static_cast<long int>(var);
        else mysol[(*__it)->_vertex]=-(static_cast<long int>(var));
    }
    for (unsigned i=1; i<=_N; ++i) {
        if(mysol[i]!=0)outfilesol_<<mysol[i]<<" ";
    }
    outfilesol_<<_seed<<" "<<_comp_init/(double)_N<<" "<<endl;
    
    
}



/*public member class Graph. This member prints the variable surveys at the first iteration*/
void Graph::print_surveys(){
    /*output file*/
    const string title="Surveys";
    const string txt=".dat";
    string str;
    ostringstream _seeds, _Ks, _sN, _salpha;
    _seeds<<_seed;
    _Ks<<_K;
    _salpha<<_alpha;
    _sN<<_N;
    str=directory;
    str+=title;
    str+="K=";
    str+=_Ks.str();
    str+="N=";
    str+=_sN.str();
    str+="alpha=";
    str+=_salpha.str();
    //str+=sat;
    //str+=_seeds.str();
    //str+=_Ks.str();
    //str+=sat;
    //str+=_seeds.str();
    str+=txt;
    ofstream outfilewff(str.c_str(), ios_base::app);
    vector<long int> mysol;
    mysol.resize(_N+1);
    for (list<Vertex*>::iterator __it=_list_fixed_element.begin(); __it!=_list_fixed_element.end(); ++__it) {
        unsigned int var=(*__it)->_vertex;
        if((*__it)->_who_I_am)mysol[(*__it)->_vertex]=static_cast<long int>(var);
        else mysol[(*__it)->_vertex]=-(static_cast<long int>(var));
    }
    for (unsigned i=0; i<_N; ++i) {
        if(mysol[i+1]!=0)outfilewff<<_v_sT[i]<<" "<<_v_sI[i]<<" "<<_v_sF[i]<<" ";
    }
    outfilewff<<_seed<<" "<<_comp_init/(double)_N<<endl;
    str.clear();
    str=directory;
    str+="Complexity_var";
    str+="K=";
    str+=_Ks.str();
    str+="N=";
    str+=_sN.str();
    str+="alpha=";
    str+=_salpha.str();
    //str+=sat;
    //str+=_seeds.str();
    //str+=_Ks.str();
    //str+=sat;
    //str+=_seeds.str();
    str+=txt;
    ofstream outfilewff_comp(str.c_str(), ios_base::app);
    for (unsigned i=0; i<_N; ++i) {
        if(mysol[i+1]!=0)outfilewff_comp<<_v_comp_var[i]<<" ";
    }
    outfilewff_comp<<_seed<<" "<<_comp_init/(double)_N<<endl;
}
/*public member class Graph. This member prints the values on file*/
void Graph::print(){
    const string title="Value_Anal_Compl";
    const string txt=".txt";
    string str;
    ostringstream _seeds, _Ks, _sN, _salpha;
    _seeds<<_seed;
    _Ks<<_K;
    _salpha<<_alpha;
    _sN<<_N;
    str=directory;
    str=title;
    str+="_seed_";
    str+=_seeds.str();
    str+="K=";
    str+=_Ks.str();
    str+="N=";
    str+=_sN.str();
    str+="alpha=";
    str+=_salpha.str();
    //str+=sat;
    //str+=_seeds.str();
    //str+=_Ks.str();
    //str+=sat;
    //str+=_seeds.str();
    str+=txt;
    ofstream outfilecomp_(str.c_str(), ios_base::app);
    for (unsigned i=0; i<_v_Nt.size(); ++i) {
        outfilecomp_<<_v_Nt[i]<<" "<<_v_M_t[i]<<" "<<_v_c[i]<<" "<<_v_time[i]<<" "<<_N<<endl;
    }
    
}

/*public member class Graph. In this member an instance of the problem is built and is printed on file*/
void Graph::write_on_file_graph(){/*build a graph and write a CNF formula*/
    /* get the values from INPUT and store them in the right place*/
    _N=static_cast<unsigned int>(stoul(_argv[_argc-1],nullptr,0));
    _alpha=stod(_argv[_argc-2],nullptr);
    double m=((double)_N*_alpha);
    _M=static_cast<unsigned int>(m);
    while(static_cast<double>(_M)<m){
        if(((double)_M)==m)break;
        _M++;
    }
    _alpha=(double)_M/(double)_N;
    _K=static_cast<unsigned int>(stoul(_argv[_argc-3],nullptr,0));
    cout<<setprecision(9);
    cout<<"I am going to build an instance with N: "<<_N<<" variables and "<<_M<<" clauses with a clause density equal to "<<_alpha<<endl;
    _ivec.reserve(_K*_M);
    _N_t=_N;
    _M_t=0;
    /*string for output file*/
    const string title="Formulas_training_CNF";
    const string sat="-SAT_seed=";
    const string txt=".cnf";
    string str;
    ostringstream _seeds, _Ks, _sN, _salpha;
    _seeds<<_seed;
    _Ks<<_K;
    _salpha<<_alpha;
    _sN<<_N;
    str=directory;
    str+=title;
    str+="K=";
    str+=_Ks.str();
    str+="N=";
    str+=_sN.str();
    str+="alpha=";
    str+=_salpha.str();
    str+=sat;
    str+=_seeds.str();
    str+=txt;
    
    /*I am modifing how to write the output file*/
#ifdef PRINT_FORMULA_CNF
    /*output file*/
    ofstream outfilewff(str.c_str(), ios_base::app);
    outfilewff << "c seed="<<_seed<<endl;
    outfilewff << "p cnf";
    outfilewff << ' ' << _N << ' ' << _M << endl;
#endif
    /******************************************************/
    /******************************************************/
    /********************** START *************************/
    /******** Henry Kautz & Bart Selman makewff.c *********/
    /******************************************************/
    /******************************************************/
    /******************************************************/
    
    /*minimal modifiications for compatibility in c++*/
    
    int i, j, k;
    int lit;
    int cl[MAX_CLEN];
    bool dup;
    
    /*build the CNF instance*/
    for (i=0; i<_M; i++){
        
        for (j=0; j<_K; j++){
            
            do {
                
                
                lit =  static_cast<int>(random() % _N) + 1;
                
                
                dup = false;
                
                for (k=0; k<j; k++)
                    
                    if (lit == cl[k]) dup = true;
                
            } while(dup);
            
            cl[j] = lit;
            
        }
        /* flip the literal*/
        for (j=0; j<_K; j++){
            
            if (_flip()) cl[j] *= -1;
            
            _ivec.push_back(cl[j]);
            
        }
        _ivec.push_back(0);
        /******************************************************/
        /******************************************************/
        /******************************************************/
        /******** Henry Kautz & Bart Selman makewff.c *********/
        /*********************** END **************************/
        /******************************************************/
        /******************************************************/
        
        /* print on file*/
#ifdef PRINT_FORMULA_CNF
        for (j=0; j<_K; j++)outfilewff<<cl[j]<<" ";
        
        outfilewff<<"0"<<" "<<endl;
#endif
    }
   
}

/*public member class Graph. This member reads an instance of a problem from a file*/
void Graph::read_from_file_graph(){/*read an instance given as INPUT*/
    /*open file to read*/
    ifstream infile(_argv[_argc-1].c_str());
    if (!infile) {
        cout <<"File not found"<<endl;
        cout<<"Please check if the name is correct or if the file exists."<<endl;
        exit(-1);
    }else{
        string ogg1;
        int ogg2, ogg3;
        bool flag=false;
        int var1;
        /*read file and store variables in _ivec*/
        /* A CNF formula is composed by by +- numbers. 0 identifies the end of a clause*/
        while (!infile.eof()) {
            if(infile.eof())break;
            if (!flag) {
                infile>>ogg1;
                if (ogg1=="c"){
                    infile>>ogg1;
                    string str;
                    str=ogg1.substr(5);
                    _seed_out=atoi(str.c_str());
                }
                if (ogg1== "cnf") {
                    infile>>ogg2>>ogg3;
                    _N=ogg2;
                    _M=ogg3;
                    _N_t=_N;
                    _M_t=_M;
                    _ivec.reserve(_M*10);
                    flag=true;
                }
            }
            if (flag) {
                infile>>var1;
                if(infile.eof())break;
                _ivec.push_back(var1);
            }
        }
        vector<long int>(_ivec).swap(_ivec);
    }
}

/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
/******************************  END GRAPH CLASS  **********************************/
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/

