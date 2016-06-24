#include <iostream>
#include <vector>
#include <initializer_list>

using std::string;
using std::vector;

class Foo {

public:

	vector<string> my_a, my_b;

	Foo(vector<string> a, vector<string> b): my_a(a), my_b(b) {}

};

std::ostream & operator<<(std::ostream &out, vector<string> const &vec) {
	for(auto x: vec) {
		out << x << " ";
	}
	return out;
}

int main() {
	Foo f({"hello", "world"}, {"a", "b", "c"});
	std::cout << f.my_a << std::endl;
	std::cout << f.my_b << std::endl;
	return 0;
}




