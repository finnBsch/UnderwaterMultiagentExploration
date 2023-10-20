//
// Created by finn on 5/4/23.
// Only include in source file if building a library!!
//

#ifndef GMRF_MSG_ASSERT_H
#define GMRF_MSG_ASSERT_H
#include <iostream>

// Assert with Message. Will only be compiled in DEBUG MODE. Can throw message.
#ifndef NDEBUG
#   define M_Assert(Expr, Msg) \
    __M_Assert(#Expr, Expr, __FILE__, __LINE__, Msg)
#else
#   define M_Assert(Expr, Msg) ((void)0);
#endif

void __M_Assert(const char* expr_str, bool expr, const char* file, int line, const char* msg)
{
    if (!expr)
    {
        std::cerr << "Assert failed:\t" << msg << "\n"
                  << "Expected:\t" << expr_str << "\n"
                  << "Source:\t\t" << file << ", line " << line << "\n";
        abort();
    }
}
#endif //GMRF_MSG_ASSERT_H
