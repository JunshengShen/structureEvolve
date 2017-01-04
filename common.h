
#pragma once
#ifdef __GNUG__
     // g++, clang++: [[ gnu::always_inline ]] does not also imply C++ inline
     // [[ gnu::always_inline ]] is an additional implementation-defined
     // attribute using the C++11 unified standard attribute syntax. note: gnu::
     #define FORCEINLINE [[ gnu::always_inline ]] inline

#elif defined _MSC_VER
     // msc++: ____forceinline implies C++ inline
     // msc++ does not provide any additional implementation-defined
     // attributes using the C++11 unified standard attribute syntax
     #define FORCEINLINE __forceinline // inline is implied

#else
     #define FORCEINLINE inline
#endif
#include <vector>