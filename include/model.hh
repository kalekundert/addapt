#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "utils.hh"

namespace addapt {

class Device;
using DevicePtr = std::shared_ptr<Device>;
using DeviceConstPtr = std::shared_ptr<Device const>;

class Aptamer;
using AptamerPtr = std::shared_ptr<Aptamer>;
using AptamerConstPtr = std::shared_ptr<Aptamer const>;

struct Context;
using ContextPtr = std::shared_ptr<Context>;
using ContextConstPtr = std::shared_ptr<Context const>;

class Device {

public:

	/// @brief Default deviceor.
	Device(string);

	/// @brief Return the sequence of this device.
	string seq() const;

	/// @brief Return the nucleotide at the given position of this device.
	char seq(int) const;

	/// @brief Return constraints that define the given macrostate.
	string macrostate(string) const;

	/// @brief Return all the macrostates associated with this device.
	unordered_map<string,string> macrostates() const;

	/// @brief Add constraints that define a particular macrostate.
	void add_macrostate(string, string);

	/// @brief Return the length of this device.
	int len() const;

	/// @brief Make a point mutation in this device.
	void mutate(int, char const);

	/// @brief Return a deep-copy of this device.
	DevicePtr copy() const;

	/// @brief Make this device equivalent to the given one.
	void assign(DevicePtr);

private:

	string my_seq;
	unordered_map<string,string> my_macrostates;
	ContextPtr my_context;

};

class Aptamer {

public:

	/// @brief Device with a sequence, fold, and ΔG (kcal/mol).
	Aptamer(string, string, double);

	/// @brief Return the sequence of this aptamer.
	string seq() const;

	/// @brief Return a pseudo-dot-bracket string specifying the base pairs 
	/// formed by this aptamer in the holo condition.
	string fold() const;

	/// @brief Return the affinity of this aptamer for its ligand (i.e. its  
	/// dissociation constant in units of μM).
	double affinity() const;

private:

	string my_seq;
	string my_fold;
	double my_affinity;

};

struct Context {
	string before;
	string after;
};


}
