#include <itpp/itbase.h>
#include <itpp/base/vec.h>
#include <itpp/signal/transforms.h>
#include <itpp/stat/misc_stat.h>
#include <itpp/signal/freq_filt.h>
#include <math.h>
#include <list>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <boost/math/special_functions/gamma.hpp>
#include <sys/time.h>
#include <curses.h>
#include "common.h"
#include "lte_lib.h"
#include "constants.h"
#include "itpp_ext.h"
#include "dsp.h"

using namespace itpp;
using namespace std;

cvec raw2iq(const vector<int8_t> &raw)
{
    size_t num_samples = raw.size() / 2;
    cvec iq(num_samples);

    for (size_t i = 0; i < num_samples; ++i)
    {
        iq[i] = complex<double>(raw[2 * i] / 128.0, raw[2 * i + 1] / 128.0);
    }

    return iq;
}

// 定义 sinc 函数，模拟 MATLAB 中的 sinc 函数
double sinc(double x)
{
    return (x == 0) ? 1.0 : (std::sin(x) / x);
}

cvec zadoff_chu(int l, int i)
{
    cvec y(l);
    if (l & 1)
    { // 如果 l 是奇数
        for (int n = 0; n < l; ++n)
        {
            y[n] = std::exp(-std::complex<double>(0, 1) * i * M_PI / l * n * (1 + n));
        }
    }
    else
    {
        for (int n = 0; n < l; ++n)
        {
            y[n] = std::exp(-std::complex<double>(0, 1) * i * M_PI / l * n * n);
        }
    }
    return y;
}

cvec pss(int n_id_2)
{
    int lut[] = {25, 29, 34};
    cvec y = zadoff_chu(63, lut[n_id_2 + 1]);
    y.del(32);
    return y;
}

void pss_gen(cmat &td_pss)
{
    cmat fd_pss(128, 3);
    // cmat td_pss(128 + 9, 3);
    for (size_t i = 0; i < 3; i++)
    {
        cvec temp = pss(i); // 假设 pss 函数返回 vec 类型
        cvec t = to_cvec(vec(0), vec(0));
        fd_pss.set_col(i, concat(t, temp.get(31, 128 - 1), zeros_c(65), temp.get(0, 30)));
        cvec temp_td = ifft(fd_pss.get_col(i)) * std::sqrt(128.0 / 62.0);
        td_pss.set_col(i, concat(temp_td.get(128 - 8, 128 - 1), temp_td.get(0, 128 - 9)));
    }
}

cmat pss_fo_set_gen(cmat &pss, const itpp::vec &fo_search_set)
{
    size_t num_pss = pss.cols();
    size_t len_pss = pss.rows();

    double sampling_rate = 1.92e6;
    size_t num_fo = fo_search_set.length();

    cmat pss_fo_set = zeros_c(len_pss, num_fo * num_pss);

    for (size_t i = 0; i < num_pss; i++)
    {
        size_t sp = i * num_fo;
        size_t ep = sp + num_fo - 1;
        for (size_t col = sp; col <= ep; col++)
        {
            cmat tmp = kron(ones_c(1, num_fo), pss(0, pss.rows() - 1, i, i)) * itpp::exp(std::complex<double>(0, 1) * 2 * M_PI * (1.0 / sampling_rate) * itpp::linspace(0, len_pss - 1, len_pss) * fo_search_set);
            pss_fo_set.set_submatrix(0, sp, conj(tmp) / len_pss);
        }
    }
    return pss_fo_set;
}

cvec shift(const cvec &vec, int len)
{
    cvec temp(vec(len, vec.size() - 1));
    return concat(temp, vec(0, len - 1));
}

mat fft_corr(cvec s, cmat &td_pss, const vec &fo_search_set)
{
    double sampling_rate = 1.92e6;
    size_t len = s.size();
    double freq_step = sampling_rate / len;

    size_t len_pss = td_pss.rows();
    size_t num_pss = td_pss.cols();
    size_t num_fo = fo_search_set.size();
    size_t num_fo_pss = num_fo * num_pss;
    mat corr_store(len, num_fo_pss);

    s = fft(s);
    cmat fd_pss(len_pss, num_fo_pss);
    for (int i = 0; i < num_fo_pss; ++i)
    {
        fd_pss.set_col(i, itpp::fft(itpp::conj(itpp::reverse(td_pss.get_row(i)) / (int)len_pss), len_pss));
    }
    for (int i = 0; i < num_fo_pss; ++i)
    {
        int pss_idx = (i / num_fo) + 1;
        int fo_idx = i - (pss_idx - 1) * num_fo;
        double fo = fo_search_set(fo_idx);
        double fd_fo_shift_len = fo / freq_step;

        // Perform circular shift and element-wise multiplication
        cvec tmp = s.transpose() * shift(fd_pss.get_col(pss_idx - 1), fd_fo_shift_len);
        corr_store.set_col(i, sqrt(itpp::abs(ifft(tmp))));
    }
    corr_store = corr_store(len_pss, corr_store.rows() - 1, 0, corr_store.cols() - 1);
    return corr_store;
}

vec absx2(const cvec &capbuf)
{
    return elem_mult(real(capbuf), real(capbuf)) + elem_mult(imag(capbuf), imag(capbuf));
}

// mat reshape(vec &sp, int m, int n)
// {
//     mat result(m, n);
//     for (int i = 0; i < n; i++)
//     {
//         result.set_row(i, sp(i * m, i * m + n));
//     }
//     return result;
// }

vec sum_by_columns(mat matrix)
{
    vec result(matrix.cols());
    for (size_t i = 0; i < matrix.cols(); i++)
    {
        result[i] = sum(matrix.get_col(i));
    }
    return result;
}