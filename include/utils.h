#ifndef HAVE_UTILS
#define HAVE_UTILS

#include <vector>
#include <itpp/itbase.h>

itpp::cvec raw2iq(const std::vector<int8_t> &raw);

// 定义 sinc 函数，模拟 MATLAB 中的 sinc 函数
double sinc(double x);

itpp::cvec zadoff_chu(int l, int i);

itpp::cvec pss(int n_id_2);

void pss_gen(itpp::cmat &td_pss);

itpp::cmat pss_fo_set_gen(itpp::cmat &pss, const itpp::vec &fo_search_set);

itpp::cvec shift(const itpp::cvec &vec, int len);

itpp::mat fft_corr(itpp::cvec s, itpp::cmat &td_pss, const itpp::vec &fo_search_set);

itpp::vec absx2(const itpp::cvec &capbuf);

// mat reshape(vec &sp, int m, int n);

itpp::vec sum_by_columns(itpp::mat matrix);

#endif