#include <cmath>
#include <iostream>
#include <boost/format.hpp>

extern "C" {
  #include <ViennaRNA/part_func.h>
  #include <ViennaRNA/structure_utils.h>
}

using f = boost::format;

int main(int argc, char **argv) {
	string seq = "ACGUGAAAACGU";
	auto *fc = vrna_fold_compound(seq.c_str(), NULL, VRNA_OPTION_PF);
	int N = fc->length;

	double rt_37 = 1.9858775e-3 * 310;  // kcal/mol at 37°C
	double aptamer_kd = 0.32;           // μM
	double std_conc = 1e6;              // 1M in μM
	double aptamer_dg = rt_37 * log(aptamer_kd / std_conc);

	vrna_sc_add_hi_motif(
			fc,
			"GAUACCAGCCGAAAGGCCCUUGGCAGC",
			"(...((.(((....)))....))...)",
			aptamer_dg,
			VRNA_OPTION_DEFAULT);

	vrna_pf(fc, NULL);

	cout << "Sequence:" << endl;
	cout << seq << endl;
	cout << endl;

	cout << "plist (eigen):" << endl;

	auto *bppm = fc->exp_matrices->probs;
	float const cutoff = 0.001;

	for(int i = 1; i < N; i++) {
		for(int j = i+1; j <= N; j++) {
			double bpp = bppm[fc->iindx[i] - j];
			if (bpp > cutoff) {
				cout << f("(%d,%d): %f") % i % j % bpp << endl;
			}
		}
	}

	cout << endl;
	cout << "plist (vrna):" << endl;

	auto *plist = vrna_plist_from_probs(fc, cutoff);

	for(int i = 0; true; i++) {
		auto bp = plist[i];
		if (bp.i == 0) break;
		cout << f("(%d,%d): %f") % bp.i % bp.j % bp.p << endl;
	}





	return 0;
}


