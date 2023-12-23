

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
#include "utils.h"

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
    {
        y = itpp::exp(elem_mult(-std::complex<double>(0, 1) * i * M_PI / l * linspace(0, l - 1, l), to_cvec(linspace(1, l, l))));
    }
    else
    {
        cvec t = -std::complex<double>(0, 1) * i * M_PI / l * linspace(0, l - 1, l);
        y = itpp::exp(elem_mult(t, t));
    }
    return y;
}

cvec pss(int n_id_2)
{
    int lut[] = {25, 29, 34};
    cvec y = zadoff_chu(63, lut[n_id_2]);
    y.del(31);
    return y;
}

void pss_gen(cmat &td_pss)
{
    cmat fd_pss(128, 3);
    for (size_t i = 0; i < 3; i++)
    {
        cvec temp = pss(i); // 假设 pss 函数返回 vec 类型
        cvec t = zeros_c(1);
        cvec row = concat(t, temp.get(31, -1), zeros_c(65), temp.get(0, 30));
        fd_pss.set_col(i, row);
        cvec temp_td = idft(fd_pss.get_col(i)) * std::sqrt(128.0 / 62.0);
        td_pss.set_col(i, concat(temp_td.get(128 - 9, 128 - 1), temp_td));
    }
}

cmat pss_fo_set_gen(cmat &pss, const itpp::vec &fo_search_set)
{
    size_t num_pss = pss.cols();
    size_t len_pss = pss.rows();

    double sampling_rate = 1.92e6;
    size_t num_fo = fo_search_set.length();

    cmat pss_fo_set = zeros_c(len_pss, num_fo * num_pss);

    cvec tmp = cvec(len_pss);
    for (size_t i = 0; i < len_pss; i++)
    {
        tmp[i] = complex<double>(0, 2 * M_PI * (1.0 / sampling_rate) * i);
    }

    for (size_t i = 0; i < num_pss; i++)
    {
        size_t sp = i * num_fo;
        size_t ep = sp + num_fo - 1;
        pss_fo_set.set_submatrix(0, sp, conj(elem_mult(kron(ones_c(1, num_fo), cmat(pss.get_col(i))), exp(outer_product(tmp, to_cvec(fo_search_set))))) / len_pss);
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
    cmat fd_pss(len, num_pss);
    for (int i = 0; i < num_pss; ++i)
    {
        fd_pss.set_col(i, itpp::fft(itpp::conj(itpp::reverse(td_pss.get_col(i)) / (int)len_pss), len));
    }
    cout << "fd_pss" << fd_pss(0, 0) << fd_pss(0, 1) << fd_pss(1, 0) << fd_pss(1, 1) << endl;

    for (int i = 0; i < num_fo_pss; ++i)
    {
        int pss_idx = floor(i / num_fo);
        int fo_idx = i - pss_idx * num_fo;
        double fo = fo_search_set(fo_idx);
        double fd_fo_shift_len = fo / freq_step;

        // Perform circular shift and element-wise multiplication
        vec tmp = abs(ifft(elem_mult(s, shift(fd_pss.get_col(pss_idx), fd_fo_shift_len)), len));
        corr_store.set_col(i, elem_mult(tmp, tmp));
    }
    corr_store = corr_store(len_pss, -1, 0, -1);
    cout << "corr_store" << corr_store(0, 0) << corr_store(0, 1) << corr_store(1, 0) << corr_store(1, 1) << endl;
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

void tshiftmy(vec &x, int t)
{
    vec y(x);
    if (t == floor(t))
    {
        t = mod(t, length(x));
        for (int i = 0; i < length(x) - t; ++i)
        {
            y[t + i] = x[i];
        }

        // Second assignment: y(1:t) = x(end-t+1:end);
        for (int i = 0; i < t; ++i)
        {
            y[i] = x[length(x) - t + i];
        }
        return;
    }
}