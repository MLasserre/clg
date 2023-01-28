#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <set>
#include <algorithm>

#include <agrum/agrum.h>
#include <agrum/BN/BayesNet.h>
#include <agrum/tools/multidim/potential.h>
#include <agrum/tools/variables/continuousVariable.h>
#include <agrum/tools/variables/rangeVariable.h>
#include <agrum/tools/variables/labelizedVariable.h>
#include <agrum/tools/core/set.h>

#include "gaussian.h"
#include "linear_gaussian.h"
#include "canonical_form.h"
#include "continuous_variable.h"
#include "scope.h"
#include "variable_elimination.h"

using namespace std;
using namespace Eigen;

int main(void){

    // Variables discrètes du problème
    auto F = gum::LabelizedVariable("F", "Filter State", {"intact","defect"});
    auto W = gum::LabelizedVariable("W", "Waste Type", {"industrial","household"});
    auto B = gum::LabelizedVariable("B", "Burning Regimen", {"stable","unstable"});

    // Variables continues du problème
    ContinuousVariable M_in("M_in", "Metals in Waste");
    ContinuousVariable M_out("M_out", "Metals Emission");
    ContinuousVariable E("E", "Filter Efficiency");
    ContinuousVariable D("D", "Dust Emission");
    ContinuousVariable C("C", "CO2 Concentration in Emission");
    ContinuousVariable L("L", "Light Penetrability");

    // Tables de formes canoniques associées à chaque variable
    gum::Potential<CanonicalForm> F_potential;
    gum::Potential<CanonicalForm> W_potential;
    gum::Potential<CanonicalForm> B_potential;

    gum::Potential<CanonicalForm> M_in_potential;
    gum::Potential<CanonicalForm> M_out_potential;
    gum::Potential<CanonicalForm> E_potential;
    gum::Potential<CanonicalForm> D_potential;
    gum::Potential<CanonicalForm> C_potential;
    gum::Potential<CanonicalForm> L_potential;


    // Ajout des variables discrètes aux potentiels
    F_potential.add(F);
    W_potential.add(W);
    B_potential.add(B);

    M_in_potential.add(W);

    E_potential.add(F);
    E_potential.add(W);

    D_potential.add(B);
    D_potential.add(W);

    C_potential.add(B);


    // REMPLISSAGE DES TABLES

    // Potentiel associé à la variable F
    gum::Instantiation I_F(F_potential);
    F_potential.set(I_F.chgVal("F", "intact"), CanonicalForm(Gaussian(0.95)));
    F_potential.set(I_F.chgVal("F", "defect"), CanonicalForm(Gaussian(0.05)));
    cout << "Potentiel de F : " << F_potential << endl;

    // Potentiel associé à la variable W
    gum::Instantiation I_W(W_potential);
    W_potential.set(I_W.chgVal("W", "industrial"), CanonicalForm(Gaussian(2./7.)));
    W_potential.set(I_W.chgVal("W", "household"), CanonicalForm(Gaussian(5./7.)));
    cout << "Potentiel de W : " << W_potential << endl;

    // Potentiel associé à la variable B
    gum::Instantiation I_B(B_potential);
    B_potential.set(I_B.chgVal("B", "stable"), CanonicalForm(Gaussian(0.85)));
    B_potential.set(I_B.chgVal("B", "unstable"), CanonicalForm(Gaussian(0.15)));
    cout << "Potentiel de B : " << B_potential << endl;

    // Potentiel associé à la variable E
    gum::Instantiation I_E(E_potential);
    E_potential.set(I_E.chgVal("F", "intact").chgVal("W", "household"),
                    CanonicalForm(Gaussian(Scope({E}),
                                           Matrix<double,1,1>(0.00002),
                                           Matrix<double,1,1>(-3.2))));
    E_potential.set(I_E.chgVal("F", "defect").chgVal("W", "household"),
                    CanonicalForm(Gaussian(Scope({E}),
                                           Matrix<double,1,1>(0.0001),
                                           Matrix<double,1,1>(-0.5))));
    E_potential.set(I_E.chgVal("F", "intact").chgVal("W", "industrial"),
                    CanonicalForm(Gaussian(Scope({E}),
                                           Matrix<double,1,1>(0.00002),
                                           Matrix<double,1,1>(-3.9))));
    E_potential.set(I_E.chgVal("F", "defect").chgVal("W", "industrial"),
                    CanonicalForm(Gaussian(Scope({E}),
                                           Matrix<double,1,1>(0.0001),
                                           Matrix<double,1,1>(-0.4))));
    cout << "Potentiel de E : " << E_potential << endl;

    // Potentiel associé à la variable D
    gum::Instantiation I_D(D_potential);
    D_potential.set(I_D.chgVal("B", "stable").chgVal("W", "industrial"),
                    CanonicalForm(LinearGaussian(Scope({D}),
                                                 Scope({E}),
                                                 Matrix<double,1,1>(0.03),
                                                 Matrix<double,1,1>(1),
                                                 Matrix<double,1,1>(6.5))));
    D_potential.set(I_D.chgVal("B", "stable").chgVal("W", "household"),
                    CanonicalForm(LinearGaussian(Scope({D}),
                                                 Scope({E}),
                                                 Matrix<double,1,1>(0.04),
                                                 Matrix<double,1,1>(1),
                                                 Matrix<double,1,1>(6.0))));
    D_potential.set(I_D.chgVal("B", "unstable").chgVal("W", "industrial"),
                    CanonicalForm(LinearGaussian(Scope({D}),
                                                 Scope({E}),
                                                 Matrix<double,1,1>(0.1),
                                                 Matrix<double,1,1>(1),
                                                 Matrix<double,1,1>(7.5))));
    D_potential.set(I_D.chgVal("B", "unstable").chgVal("W", "household"),
                    CanonicalForm(LinearGaussian(Scope({D}),
                                                 Scope({E}),
                                                 Matrix<double,1,1>(0.1),
                                                 Matrix<double,1,1>(1),
                                                 Matrix<double,1,1>(7.0))));
    cout << "Potentiel de D : " << D_potential << endl;

    // Potentiel associé à la variable C
    gum::Instantiation I_C(C_potential);
    C_potential.set(I_C.chgVal("B", "stable"),
                    CanonicalForm(Gaussian(Scope({C}),
                                           Matrix<double,1,1>(0.1),
                                           Matrix<double,1,1>(-2))));
    C_potential.set(I_C.chgVal("B", "unstable"),
                    CanonicalForm(Gaussian(Scope({C}),
                                  Matrix<double,1,1>(0.3),
                                  Matrix<double,1,1>(-1))));
    cout << "Potentiel de C : " << C_potential << endl;

    // Potentiel associé à la variable M_in
    gum::Instantiation I_M_in(M_in_potential);
    M_in_potential.set(I_M_in.chgVal("W", "industrial"),
                       CanonicalForm(Gaussian(Scope({M_in}),
                                              Matrix<double,1,1>(0.1),
                                              Matrix<double,1,1>(0.5))));
    M_in_potential.set(I_M_in.chgVal("W", "household"),
                       CanonicalForm(Gaussian(Scope({M_in}),
                                              Matrix<double,1,1>(0.005),
                                              Matrix<double,1,1>(-0.5))));
    cout << "Potentiel de M_in : " << endl << M_in_potential << endl;

    // Potentiel associé à la variable L
    gum::Instantiation I_L(L_potential);
    L_potential.set(I_L, CanonicalForm(LinearGaussian(Scope({L}),
                                                      Scope({D}),
                                                      Matrix<double,1,1>(0.25),
                                                      Matrix<double,1,1>(-1),
                                                      Matrix<double,1,1>(3))));
    cout << "Potentiel de L : " << L_potential << endl;

    // Potentiel associé à la variable L
    gum::Instantiation I_M_out(M_out_potential);
    M_out_potential.set(I_M_out, CanonicalForm(
                                     LinearGaussian(Scope({M_out}),
                                                    Scope({D, M_in}),
                                                    Matrix<double,1,1>(0.002),
                                                    Matrix<double,2,1>(1,1),
                                                    Matrix<double,1,1>(0))));
    cout << "Potentiel de M_out : " << M_out_potential << endl;

    cout << "Somme potentiel W : " << W_potential.sum() << endl;

    // Creation d'un vecteur contenant tous les potentiels.
    vector< gum::Potential< CanonicalForm > > potential_set
        = {F_potential, E_potential, W_potential, M_in_potential, M_out_potential,
           D_potential, B_potential, C_potential, L_potential};

    // Création d'un BN fantome pour trouver l'ordre d'élimination
    auto BN = gum::BayesNet<double>::fastPrototype(
                  "F->E->D->L;W->D->M_out;W->M_in->M_out;W->E;B->D;B->C");
    cout << BN << endl;

    auto ce = ContinuousEvidence({{L, 0.1}});
    auto de = DiscreteEvidence({{&W,"industrial"}});
    auto cOrder = vector<ContinuousVariable>({});
    auto dOrder = vector<gum::LabelizedVariable*>({&F});

    auto res = sum_product_ve(cOrder, dOrder, ce, de, potential_set);
    //cout << res << endl;

    //gum::Set<const gum::LabelizedVariable*> s;
    //s.insert(&F);
    //cout << "Test : " << F_potential.margSumOut(s) << endl;

 /*

    ContinuousVariable X("X");
    ContinuousVariable Y("Y");
    ContinuousVariable Z("Z");

    auto D = gum::LabelizedVariable("D", "D", {"True","False"});

    gum::Potential<CanonicalForm> X_potential;
    gum::Potential<CanonicalForm> Y_potential;
    gum::Potential<CanonicalForm> Z_potential;
    gum::Potential<CanonicalForm> D_potential;

    X_potential.add(D);
    D_potential.add(D);

    gum::Instantiation I_X(X_potential);
    gum::Instantiation I_Y(Y_potential);
    gum::Instantiation I_Z(Z_potential);
    gum::Instantiation I_D(D_potential);

    X_potential.set(I_X.chgVal("D", "True"),
                       CanonicalForm(Gaussian(Scope({X}),
                                              Matrix<double,1,1>(1),
                                              Matrix<double,1,1>(-1))));
    X_potential.set(I_X.chgVal("D", "False"),
                       CanonicalForm(Gaussian(Scope({X}),
                                              Matrix<double,1,1>(1),
                                              Matrix<double,1,1>(1))));
    
    Y_potential.set(I_Y, CanonicalForm(Gaussian(Scope({Y}),
                                                Matrix<double,1,1>(9),
                                                Matrix<double,1,1>(20))));
    
    Z_potential.set(I_Z, CanonicalForm(LinearGaussian(Scope({Z}),
                                                      Scope({Y}),
                                                      Matrix<double,1,1>(10),
                                                      Matrix<double,1,1>(-100),
                                                      Matrix<double,1,1>(50))));

    D_potential.set(I_D.chgVal("D", "True"), CanonicalForm(Gaussian(0.5)));
    D_potential.set(I_D.chgVal("D", "False"), CanonicalForm(Gaussian(0.5)));


    vector<gum::Potential< CanonicalForm > > potential_set = {X_potential,
                                                              D_potential};

    vector<gum::Potential< CanonicalForm > > pot_set = {Y_potential,
                                                              Z_potential};
    //auto cev = ContinuousEvidence({{Z, 18.}});
    //auto dev = DiscreteEvidence();
    //auto coOrder = vector<ContinuousVariable>();
    //auto diOrder = vector<gum::LabelizedVariable>();

    auto ce = ContinuousEvidence({{L, 0.1}});
    auto de = DiscreteEvidence({{W,"industrial"}});
    auto cOrder = vector<ContinuousVariable>({});
    auto dOrder = vector<gum::LabelizedVariable>();

    auto res = sum_product_ve(cOrder, dOrder, ce, de, potential_set);
    cout << res << endl;
    */

    //auto resu = sum_product_ve(coOrder, diOrder, cev, dev, pot_set);
    //gum::Instantiation I_resu(resu);
    //cout << resu.get(I_resu) << endl;

    //auto a=gum::LabelizedVariable("a","a");
    //auto b=gum::RangeVariable("b","b",0,5);
    //P.add(a);
    //P.add(b);
    //gum::Instantiation I(P);
    //cout<<I<<endl;
    //I.inc();
    //cout<<I<<endl;
    //I.inc();
    //cout<<I<<endl;

    //P.set(I,p);
    //I.inc();
    //P.set(I,q);
    //cout<<p<<std::endl;
   

    //gum::Potential<CanonicalForm> Q;
    //cout << "P : " << P << endl;
    //cout << Q << endl;
    //cout << P*Q << endl;


    //typedef Eigen::Matrix<int, Eigen::Dynamic, 1> IndexType;
    //IndexType my_perm(3);
    //my_perm << 1, 2, 0;
    //PermutationMatrix<Dynamic,Dynamic> perm(my_perm);
    //MatrixXd A(3,3);
    //A << 1, 2, 3,
         //4, 5, 6,
         //7, 8, 9;
    //cout << "A matrix " << endl << A << endl;
    //cout << "Permuted A matrix " << endl << perm.transpose()*A*perm << endl;

    return 0;
}
