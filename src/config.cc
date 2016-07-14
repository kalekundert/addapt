#include <regex>

#include <yaml-cpp/yaml.h>

#include "config.hh"
#include "model.hh"
#include "sampling.hh"
#include "scoring.hh"
#include "utils.hh"

namespace addapt {

using std::stod;
using std::stoi;

YAML::Node
find_section(vector<string> config_files, string name) {
	YAML::Node section;
	bool section_found = false;

	for(string file: config_files) {
		YAML::Node config = YAML::LoadFile(file);

		if(config[name]) {
			if(section_found) {
				throw (f("found 2 '%s' configurations") % name).str();
			}

			section = config[name];
			section_found = true;
		}
	}

	if(not section_found) {
		throw (f("no '%s' configuration") % name).str();
	}

	return section;
}

DevicePtr
device_from_yaml(vector<string> config_files) {
	DevicePtr device;

	// Load the sequence.
	YAML::Node seq_section = find_section(config_files, "sequence");
	device = make_shared<Device>(seq_section.as<string>());

	// Load any macrostates that are defined.
	YAML::Node macro_section = find_section(config_files, "macrostates");
	for(auto item: macro_section) {
		device->add_macrostate(
				item.first.as<string>(),
				item.second.as<string>());
	}

	return device;
}

ScoreFunctionPtr
scorefxn_from_yaml(vector<string> config_files) {
	ScoreFunctionPtr scorefxn = make_shared<ScoreFunction>();

	// Load score terms for the apo and holo states.
	YAML::Node obj_section = find_section(config_files, "objective");
	*scorefxn += score_term_from_str(
			ConditionEnum::APO, obj_section["apo"].as<string>());
	*scorefxn += score_term_from_str(
			ConditionEnum::HOLO, obj_section["holo"].as<string>());

	// Load the aptamer parameters.
	YAML::Node apt_section = find_section(config_files, "aptamer");
	scorefxn->aptamer(make_shared<Aptamer>(
			apt_section["sequence"].as<string>(),
			apt_section["fold"].as<string>(),
			stod(apt_section["affinity"].as<string>())));

	return scorefxn;
}

ScoreTermPtr
score_term_from_str(ConditionEnum condition, string spec) {
	std::regex pattern("(not )?(\\w+)");
	std::smatch match;

	if(std::regex_match(spec, match, pattern)) {
		return make_shared<MacrostateProbTerm>(
				match[2], condition, match[1].matched?
				FavorableEnum::NO : FavorableEnum::YES);
	}

	throw (f("can't understand objective: '%s'") % spec).str();
}

ThermostatPtr
thermostat_from_yaml(vector<string> config_files) {
	YAML::Node section = find_section(config_files, "thermostat");
	return thermostat_from_str(section.as<string>());
}

ThermostatPtr
thermostat_from_str(string spec) {
	std::regex fixed_pattern(
			"([0-9.e+-]+)" 	    // A floating-point number (the temperature).
	);
	std::regex annealing_pattern(
			"([0-9.e+-]+)"      // A floating point number (the high temperature).
			" to "
			"([0-9.e+-]+)"      // A floating point number (the low temperature).
			" in "
			"([0-9]+)"          // An integer (the number of steps per cycle).
			" steps"
	);
	std::regex auto_scaling_pattern(
			"auto"
			"(?:"						    // Optional argument.
			"\\s+"
			"([0-9.]+)%"        // A percentage (the target acceptance rate).
				"(?:"         
				"\\s+"
				"([0-9]+)"        // An integer (the training period).
					"(?:"
					"\\s+"
					"([0-9.e+-]+)"  // A floating point number (the initial temperature).
					")?"
				")?"
			")?"
	);

	std::smatch match;

	if(std::regex_match(spec, match, fixed_pattern)) {
		double temperature = stod(match[1]);
		return make_shared<FixedThermostat>(temperature);
	}

	if(std::regex_match(spec, match, annealing_pattern)) {
		double high_temp = stod(match[1]);
		double low_temp = stod(match[2]);
		int cycle_len = stoi(match[3]);
		return make_shared<AnnealingThermostat>(cycle_len, high_temp, low_temp);
	}

	if(std::regex_match(spec, match, auto_scaling_pattern)) {
		double accept_rate = stod(match[1].length()? match[1].str() : "50") / 100;
		int training_period = stod(match[2].length()? match[2].str() : "100");
		double initial_temp = stod(match[3].length()? match[3].str() : "1");
		return make_shared<AutoScalingThermostat>(
				accept_rate, training_period, initial_temp);
	}

	throw (f("can't make a thermostat from '%s'") % spec).str();
}


}
