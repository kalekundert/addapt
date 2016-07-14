#include <yaml-cpp/yaml.h>

#include "model.hh"
#include "sampling.hh"
#include "scoring.hh"
#include "utils.hh"

namespace addapt {

DevicePtr
device_from_yaml(vector<string>);

ScoreFunctionPtr
scorefxn_from_yaml(vector<string>);

ScoreTermPtr
score_term_from_str(ConditionEnum, string);

ThermostatPtr
thermostat_from_yaml(vector<string>);

ThermostatPtr
thermostat_from_str(string);


}
