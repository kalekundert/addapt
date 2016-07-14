#include <algorithm>

#include "model.hh"
#include "utils.hh"

namespace addapt {

Device::Device(string seq): my_seq(seq) {}

string
Device::seq() const {
	return my_seq;
}

char
Device::seq(int index) const {
	index = normalize_index(my_seq, index, IndexEnum::ITEM);
	return my_seq[index];
}

string
Device::macrostate(string name) const {
	if(my_macrostates.find(name) == my_macrostates.end()) {
		throw (f("no macrostate '%s'") % name).str();
	}
	return my_macrostates.at(name);
}

unordered_map<string,string>
Device::macrostates() const {
	return my_macrostates;
}

void
Device::add_macrostate(string name, string macrostate) {
	if(macrostate.length() != len()) {
		throw "constraint length doesn't match sequence length";
	}
	my_macrostates[name] = macrostate;
}

int
Device::len() const {
	return my_seq.length();
}

void
Device::mutate(int index, char const mutation) {
	index = normalize_index(my_seq, index, IndexEnum::ITEM);
	my_seq[index] = mutation;
}

DevicePtr
Device::copy() const {
	DevicePtr other = std::make_shared<Device>(my_seq);
	other->my_macrostates = my_macrostates;
	return other;
}

void
Device::assign(DevicePtr other) {
	my_seq = other->my_seq;
	my_macrostates = other->my_macrostates;
}


Aptamer::Aptamer(string seq, string fold, double affinity):
	my_seq(seq), my_fold(fold), my_affinity(affinity) {}

string
Aptamer::seq() const {
	return my_seq;
}

string
Aptamer::fold() const {
	return my_fold;
}

double
Aptamer::affinity() const {
	return my_affinity;
}


}
