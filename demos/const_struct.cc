#include <iostream>

using namespace std;

struct Foo {
	const string a;
	const string b;
};

void greet(Foo f) {
	cout << f.a << " " << f.b << endl;
}

int main() {
	Foo f = {"hello", "world"};
	f.a = "goodbye";
	greet(f);
	return 0;
}

