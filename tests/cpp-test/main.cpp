#include <algorithm>
#include <iostream>
#include <array>
using namespace std;

int main()
{
    array<double, 9> r{0,0,0,1.5, 1.6, 1.7, 0,0,0};
    for_each(&r[3], &r[5], [](double &x){
        x *= 2;
    });

    for (auto x : r)
        cout << x << ',';
    cout << endl;
}
