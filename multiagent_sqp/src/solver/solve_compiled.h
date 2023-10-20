/* This file was automatically generated by CasADi 3.6.3+.
 *  It consists of: 
 *   1) content generated by CasADi runtime: not copyrighted
 *   2) template code copied from CasADi source: permissively licensed (MIT-0)
 *   3) user code: owned by the user
 *
 */
#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

extern "C" int solve_compiled(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem);
extern "C" int solve_compiled_alloc_mem(void);
extern "C" int solve_compiled_init_mem(int mem);
extern "C" void solve_compiled_free_mem(int mem);
extern "C" int solve_compiled_checkout(void);
extern "C" void solve_compiled_release(int mem);
extern "C" void solve_compiled_incref(void);
extern "C" void solve_compiled_decref(void);
extern "C" casadi_int solve_compiled_n_in(void);
extern "C" casadi_int solve_compiled_n_out(void);
extern "C" casadi_real solve_compiled_default_in(casadi_int i);
extern "C" const char* solve_compiled_name_in(casadi_int i);
extern "C" const char* solve_compiled_name_out(casadi_int i);
extern "C" const casadi_int* solve_compiled_sparsity_in(casadi_int i);
extern "C" const casadi_int* solve_compiled_sparsity_out(casadi_int i);
extern "C" int solve_compiled_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w);
#define solve_compiled_SZ_ARG 684
#define solve_compiled_SZ_RES 696
#define solve_compiled_SZ_IW 645
#define solve_compiled_SZ_W 38468
