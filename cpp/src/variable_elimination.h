#ifndef VARIABLE_ELIMINATION_H
#define VARIABLE_ELIMINATION_H

#include <vector>
#include <map>
#include <agrum/agrum.h>
#include <agrum/tools/multidim/potential.h>
#include <agrum/tools/variables/labelizedVariable.h>
#include "continuous_variable.h"
#include "canonical_form.h"


void sum_product_eliminate_continuous_var(
        std::vector<gum::Potential< CanonicalForm > > &cf_set,
        ContinuousVariable &variable);

void sum_product_eliminate_discrete_var(
        std::vector<gum::Potential< CanonicalForm > > &cf_set,
        gum::LabelizedVariable &variable);

void reduce_discrete_var(
        std::vector< gum::Potential< CanonicalForm > > &pot_set,
        std::pair< gum::LabelizedVariable*, std::string > anEvidence);


gum::Potential<CanonicalForm> sum_product_ve(
        std::vector<ContinuousVariable> &continuous_elim_order,
        std::vector<gum::LabelizedVariable*> &discrete_elim_order,
        ContinuousEvidence &continuous_evidence,
        DiscreteEvidence &discrete_evidence,
        std::vector< gum::Potential< CanonicalForm> > cf_set);
#endif // VARIABLE_ELIMINATION_H
