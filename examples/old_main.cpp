#include <iostream>
#include <string>
#include <filesystem>

#include "transflow.hpp"

using std::cout;
using std::endl;

int main()
{
//    Material mat(1, 2, 3);
    cout << "hei" << endl;
//    Material mat = Material::concrete;
//    cout << mat.density() << endl;
    BurialMedium ans = BurialMedium::soil;
    cout << ans.conductivity() << endl;
}
