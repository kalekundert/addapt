#include <cmath>
#include <chrono>
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

static char const USAGE[] = R"""(
Insert the theophylline aptamer into the first hairpin of the sgRNA and run a 
Monte Carlo design simulation to optimize variable regions in the nexus, the 
ruler, and the rest of the first hairpin.  The design goal is to only form the 
wildtype nexus and hairpin base pairs when theophylline is bound.

Usage:
  addapt <config>... [options]

Options:
  -n --num-moves <num>           [default: 100]
    The number of moves to attempt in the design simulation.  I haven't yet 
    determined how many moves are required to reach convergence.
    
  -T --temperature <schedule>    [default: auto]
    The temperature to use for the Metropolis criterion, which affects the 
    likelihood of accepting a negative move.  If T=0, only positive moves will 
    be accepted.  In the limit that T=inf, every move will be accepted.  You 
    can specify a fixed temperature (e.g. "5"), a multi-cooled simulated 
    annealing schedule (e.g. "5 10=>0"), or schedule that tries to achieve 
    a certain acceptance rate (e.g. "auto 50%").
    
  -r --random-seed <seed>        [default: 0]
    The seed for the random number generator.  If running in parallel, this 
    should be different for each job.
    
  -o --output <path>             [default: traj.tsv]
    The path where the trajectory of the design simulation will be saved.  This 
    trajectory includes scores and sequences for every step of the simulation.
    
  -i --output-interval <steps>   [default: 1]
    How often a new snapshot in the trajectory should be recorded.
    
  -v --version
    Display the version of ``addapt`` being used.
    
  -h --help
    Display this usage information.
)""";

int main(int argc, char **argv) {
	try {
		map<string, docopt::value> args = docopt::docopt(
				USAGE+1, {argv + 1, argv + argc}, true, "0.0");
		vector<string> config_files = args["<config>"].asStringList();
		
		// Create the device.
		DevicePtr device = device_from_yaml(config_files);

		// Create the score function.
		ScoreFunctionPtr scorefxn = scorefxn_from_yaml(config_files);

		// Create the Monte Carlo sampler.
		MonteCarloPtr sampler = make_shared<MonteCarlo>();
		*sampler += make_shared<UnbiasedMutationMove>();

		ThermostatPtr thermostat = args["--temperature"]? 
			thermostat_from_str(args["--temperature"].asString()) :
			thermostat_from_yaml(config_files);

		ReporterPtr progress_bar = make_shared<ProgressReporter>();
		ReporterPtr traj_reporter = make_shared<TsvTrajectoryReporter>(
				args["--output"].asString(),
				stoi(args["--output-interval"].asString()));

		sampler->num_steps(stoi(args["--num-moves"].asString()));
		sampler->scorefxn(scorefxn);
		sampler->thermostat(thermostat);
		sampler->add_reporter(progress_bar);
		sampler->add_reporter(traj_reporter);

		std::mt19937 rng(stoi(args["--random-seed"].asString()));

		// Run the design simulation.
		sampler->apply(device, rng);
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
