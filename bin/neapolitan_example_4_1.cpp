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

    // Variables continues du problème
    ContinuousVariable WH("Work hours");   // Work hours
    ContinuousVariable In("Income");   // Income

    // Tables de formes canoniques associées à chaque variable
    gum::Potential<CanonicalForm> WH_potential;
    gum::Potential<CanonicalForm> In_potential;


    // REMPLISSAGE DES TABLES

    gum::Instantiation I_WH(WH_potential);
    WH_potential.set(I_WH, CanonicalForm(
                                     Gaussian(Scope({WH}),
                                                    Matrix<double,1,1>(5),
                                                    Matrix<double,1,1>(40))));
    cout << "Potentiel de WH : " << WH_potential << endl;

    gum::Instantiation I_In(In_potential);
    In_potential.set(I_In, CanonicalForm(
                                     LinearGaussian(Scope({In}),
                                                    Scope({WH}),
                                                    Matrix<double,1,1>(30),
                                                    Matrix<double,1,1>(10),
                                                    Matrix<double,1,1>(0))));
    cout << "Potentiel de WH : " << WH_potential << endl;

    // Creation d'un vecteur contenant tous les potentiels.
    vector< gum::Potential< CanonicalForm > > potential_set
        = {WH_potential, In_potential};

    auto ce = ContinuousEvidence({});
    auto de = DiscreteEvidence({});
    auto cOrder = vector<ContinuousVariable>({WH});
    auto dOrder = vector<gum::LabelizedVariable*>({});

    auto res = sum_product_ve(cOrder, dOrder, ce, de, potential_set);
    gum::Instantiation I_res(res);
    I_res.setFirst();

    //cout << "Resulting potential : " << Gaussian(res) << endl;
    cout << "Resulting potential : " << I_res << endl;

    return 0;
}
