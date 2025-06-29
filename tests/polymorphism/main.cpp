#include <iostream>
#include <memory>
using namespace std;

class Base1
{
public:
	Base1() = default;
	~Base1() = default;
	
	virtual void func1()
	{
		cout << "func1 from base1\n";
	}
};

class Base2 : public Base1
{
public:
	Base2() = default;
	~Base2() = default;
};

class Derived : public Base2
{
	void func1()
	{
		cout << "func1 from derived\n";
	}
};



int main()
{
	shared_ptr<Base1> obj = make_shared<Derived>();
	obj->func1();
	
	return 0;
}
