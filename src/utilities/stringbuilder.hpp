#pragma once

#include <string>
#include <sstream>

namespace utils
{

// stringbuilder from http://stackoverflow.com/a/5686975/1850917
// typical usage:
//     string s(stringbuilder() << "someString" << somethingElsePerhapsAVariable << "anotherString");
struct stringbuilder
{
   std::stringstream ss;
   template<typename T>
   stringbuilder& operator << (const T &data)
   {
        ss << data;
        return *this;
   }
   operator std::string() { return ss.str(); }
};

}
