#include <vector>
#include <itpp/itbase.h>

#ifndef HAVE_UTILS
#define HAVE_UTILS


extern itpp::cvec raw2iq(const std::vector<int8_t> &raw);

// 定义 sinc 函数，模拟 MATLAB 中的 sinc 函数
extern double sinc(double x);

extern itpp::cvec zadoff_chu(int l, int i);

extern itpp::cvec pss(int n_id_2);

extern void pss_gen(itpp::cmat &td_pss);

extern itpp::cmat pss_fo_set_gen(itpp::cmat &pss, const itpp::vec &fo_search_set);

extern itpp::cvec shift(const itpp::cvec &vec, int len);

extern itpp::mat fft_corr(itpp::cvec s, itpp::cmat &td_pss, const itpp::vec &fo_search_set);

extern itpp::vec absx2(const itpp::cvec &capbuf);

// mat reshape(vec &sp, int m, int n);

extern itpp::vec sum_by_columns(itpp::mat matrix);

extern void tshiftmy(itpp::vec &x, int t);
#endif