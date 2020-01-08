#include "heattransfer/material.hpp"

//const Material Material::concrete = Material( 2.9,    3400,   650);
//const Material Material::steel    = Material(50,      7800,   590);
//const Material Material::coating  = Material( 0.74,   1300,  1900);
//const Material Material::soil     = Material( 2,      2000,  1000); // water saturated sand

//const Material Material::seawater = Material( 0.571,  1020,  4187);
//const Material Material::air      = Material( 0.0257, 1.225, 1012);

//Name                             Density    Thermal Conductivity  Heat Capacity  Description
//                                 (kg/m3)    (W/(m2.K/m))          (J/mol.K)
//-----------------------------------------------------------------------------------------------------------------
//Stee7800-50.00                   7800,000   50,000                500,000        Steel 7800-50.00 FP
//Coat1300-00.74                   1300,000   0,7400                1900,000       Coating 1300-00.74 FP
//Conc3400-02.90                   3400,000   2,900                 650,000        Concrete 3400-02.90 FP
//Soil2000-02.00                   2000,000   2,000                 1000,000       Saturated sand 2000 2 W/mK
//Seaw1010-00.57                   1010,000   0,5710                4187,000       Sea water 1010

//Stee7800-50.00                   7800,000   50,000                500,000        Steel 7800-50.00 FP
//Coat1300-00.74                   1300,000   0,7400                1900,000       Coating 1300-00.74 FP
//Conc3400-02.90                   3400,000   2,900                 1500,000       Concrete 3400-02.90 FP
//Soil2000-02.00                   2000,000   2,200                 1000,000       Saturated sand 2000 2 W/mK
//Seaw1010-00.57                   1010,000   0,5710                4187,000       Sea water 1010
//Coat0930-00.22                   930,000    0,2200                1800,000       Polypropylene
//Air_Ambient                      1,225      0,0257                1012,000       Ambient air
