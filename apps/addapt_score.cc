#include <iostream>
#include <string>

#include <docopt/docopt.h>
#include <yaml-cpp/yaml.h>

#include "config.hh"
#include "model.hh"
#include "sampling.hh"
#include "scoring.hh"
#include "utils.hh"

using namespace std;
using namespace addapt;

// This 

static char const USAGE[] = R"""(
Calculate the score of the sequence specified in the config files.  The exact 
meaning of the score can be configured, but in general terms the more likely 
the RNA device is to fold into the desired conformations, the higher the score 
will be.

Usage:
  addapt_score <config>... [options]

Options:
  --version
    Display the version of ``addapt`` being used.
    
  -h, --help
    Display this usage information.
)""";

int main(int argc, char **argv) {
	try {
		map<string, docopt::value> args = docopt::docopt(
				USAGE+1, {argv + 1, argv + argc}, true, VERSION);
		vector<string> config_files = args["<config>"].asStringList();

		// Create the device.
		DevicePtr device = device_from_yaml(config_files);
		
		// Create the score function.
		ScoreFunctionPtr scorefxn = scorefxn_from_yaml(config_files);

		// Score the device.
		EvaluatedScoreFunction score_table;
		double score = scorefxn->evaluate(device, score_table);

		// Print the score and (in --verbose mode) the individual score terms.
		cout << score << endl;
		for(auto score_term: score_table) {
			cout << score_term.name << ":\t" << score_term.term << endl;
		}
		return 0;
	}
	catch(YAML::Exception exc) {
		cerr << "YAML Error: " << exc.msg << endl;
		return 1;
	}
	catch(string error_message) {
		cerr << "Error: " << error_message << endl;
		return 1;
	}
}
