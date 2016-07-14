#include <yaml-cpp/yaml.h>

#include "model.hh"
#include "sampling.hh"
#include "scoring.hh"
#include "utils.hh"

namespace addapt {

YAML::Node
find_section(vector<string> config_files, string name);

DevicePtr
device_from_yaml(vector<string> config_files);

ScoreFunctionPtr
scorefxn_from_yaml(vector<string> config_files);

ScoreTermPtr
score_term_from_str(ConditionEnum condition, string spec);

ThermostatPtr
thermostat_from_yaml(vector<string> config_files);

ThermostatPtr
thermostat_from_str(string spec);


}
