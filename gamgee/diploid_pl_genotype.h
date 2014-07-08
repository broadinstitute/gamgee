#ifndef __gamgee__diploid_pl_genotype__
#define __gamgee__diploid_pl_genotype__

namespace gamgee {
/**
 * @brief simple enum to keep the indices of the genotypes in the PL field of diploid individuals
 */
enum class DiploidPLGenotype { HOM_REF = 0, HET = 1, HOM_VAR = 2};

}  // end of namespace

#endif
