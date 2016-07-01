#include <iostream>
#include <regex>

using namespace std;

int main() {
	regex auto_scaling_pattern(
			"auto"
			"(?:"						    // Optional argument.
			"\\s+"
			"([0-9.]+)%"        // A percentage (the target acceptance rate).
				"(?:"         
				"\\s+"
				"([0-9.e+-]+)"    // A floating point number (the initial temperature).
					"(?:"
					"\\s+"
					"([0-9]+)"      // An integer (the training period).
					")?"
				")?"
			")?"
	);

	smatch match;
	string haystack("auto 50%");

	regex_match(haystack, match, auto_scaling_pattern);

	cout << match[1] << endl;
	cout << (match[2].length()? match[2].str() : "100") << endl;
	cout << match[3] << endl;

	return 0;
}
