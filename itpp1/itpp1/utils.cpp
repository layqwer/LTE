#include <itpp/itbase.h>
#include <itpp/base/vec.h>
#include <itpp/signal/transforms.h>
#include <itpp/stat/misc_stat.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <list>
#include <iomanip>
#include <algorithm>
#include <vector>
//// #include <boost/math/special_functions/gamma.hpp>
#include <time.h>
// #include <curses.h>
#include <regex>
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
	const size_t num_samples = raw.size() / 2;
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
    return x == 0 ? 1.0 : std::sin(x) / x;
}

cvec zadoff_chu(int len, int ii)
{
	cvec y(len);
	vec line = linspace(0, len - 1, len);
	const complex<double> val = -complex<double>(0, 1) * ii * M_PI / len;
	if (len & 1) {
		for (int i = 0; i < len; ++i) {
			y[i] = exp(val * line[i] * (line[i] + 1));
		}
	}
	else {
		itpp::cvec t(len);
		for (int i = 0; i < len; ++i) {
			complex<double> tmp = val * line[i];
			t[i] = tmp * tmp;
		}
		y = itpp::exp(t);
	}

	return y;
}

cvec pss(int n_id_2)
{
	const int lut[] = {25, 29, 34};
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
	const size_t num_pss = pss.cols();
	const size_t len_pss = pss.rows();

	const double sampling_rate = 1.92e6;
	const size_t num_fo = fo_search_set.length();

    cmat pss_fo_set = zeros_c(len_pss, num_fo * num_pss);

    cvec tmp = cvec(len_pss);
    for (size_t i = 0; i < len_pss; i++)
    {
        tmp[i] = complex<double>(0, 2 * M_PI * (1.0 / sampling_rate) * i);
    }

    for (size_t i = 0; i < num_pss; i++)
    {
	    const size_t sp = i * num_fo;
        size_t ep = sp + num_fo - 1;
        pss_fo_set.set_submatrix(0, sp, conj(elem_mult(kron(ones_c(1, num_fo), cmat(pss.get_col(i))), exp(outer_product(tmp, to_cvec(fo_search_set))))) / len_pss);
    }
    return pss_fo_set;
}

cvec shift(const cvec &vec, int len)
{
	const int size = vec.size();
    len = (len % size + size) % size; // 处理负数位移

    cvec temp(size);

    for (int i = 0; i < size; ++i)
    {
        temp[(i + len) % size] = vec[i];
    }

    return temp;
    // cvec temp(vec(abs(len), -1));
    // return concat(temp, vec(0, len - 1));
}

mat fft_corr(cvec s, cmat &td_pss, const vec &fo_search_set)
{
	const double sampling_rate = 1.92e6;
	const size_t len = s.size();
	const double freq_step = sampling_rate / len;

	const size_t len_pss = td_pss.rows();
	const size_t num_pss = td_pss.cols();
	const size_t num_fo = fo_search_set.size();
	const size_t num_fo_pss = num_fo * num_pss;
    mat corr_store(len, num_fo_pss);
    cmat corr_store_tmp(len, num_fo_pss);

    s = fft(s);
    cmat fd_pss(len, num_pss);
    for (int i = 0; i < num_pss; ++i)
    {
        fd_pss.set_col(i, itpp::fft(itpp::conj(itpp::reverse(td_pss.get_col(i)) / (int)len_pss), len));
    }
    for (int i = 0; i < num_fo_pss; ++i)
    {
	    const int pss_idx = floor(i / num_fo);
	    const int fo_idx = i - pss_idx * num_fo;
	    const double fo = fo_search_set(fo_idx);
	    const double fd_fo_shift_len = fo / freq_step;
        corr_store_tmp.set_col(i, elem_mult(s, shift(fd_pss.get_col(pss_idx), fd_fo_shift_len)));
    }
    for (int i = 0; i < num_fo_pss; i++)
    {
        vec t = abs(ifft(corr_store_tmp.get_col(i)));
        corr_store.set_col(i, elem_mult(t, t));
    }
    corr_store = corr_store(len_pss - 1, -1, 0, -1);
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

vec tshiftmy(vec &x, int t)
{
    vec y(x);
    if (t == floor(t))
    {
        t = mod(t, length(x));
        const int n = x.size();
        // Convert negative shift to positive equivalent
        t = (t % n + n) % n;
        // Circular shift to the right
        for (int i = 0; i < n - t; ++i)
        {
            y[i + t] = x[i];
        }
        // Circular shift to the left
        for (int i = 0; i < t; ++i)
        {
            y[i] = x[n - t + i];
        }
    }
    return y;
}

cvec getSubVec(cvec s, int start)
{
    return s(itpp_ext::matlab_range(start, 2, s.size()));
}

vec getSubVec(vec s, int start)
{
    return s(itpp_ext::matlab_range(start, 2, s.size()));
}

ivec mod(ivec s, int i)
{
    ivec res(s.size());
    for (size_t j = 0; j < s.size(); j++)
    {
        res[j] = mod(s[j], i);
    }
    return res;
}

ivec sss(int n_id_1, int n_id_2, int slot_num)
{
    if (slot_num != 0 && slot_num != 10)
    {
        // cout << "slot_num must be either 0 or 10" << endl;
    }
    const int qp = floor(n_id_1 / 30);
    const int q = floor((n_id_1 + qp * (qp + 1) / 2) / 30);
    const int mp = n_id_1 + q * (q + 1) / 2;
    const int m0 = mod(mp, 31);
    const int m1 = mod(m0 + floor(mp / 31) + 1, 31);

    ivec s_td = "0 0 0 0 1 0 0 1 0 1 1 0 0 1 1 1 1 1 0 0 0 1 1 0 1 1 1 0 1 0 1";
    s_td = 1 - 2 * s_td;

    ivec c_td = "0 0 0 0 1 0 1 0 1 1 1 0 1 1 0 0 0 1 1 1 1 1 0 0 1 1 0 1 0 0 1";
    c_td = 1 - 2 * c_td;

    ivec z_td = "0 0 0 0 1 1 1 0 0 1 1 0 1 1 1 1 1 0 1 0 0 0 1 0 0 1 0 1 0 1 1";
    z_td = 1 - 2 * z_td;

    const ivec s0_m0 = s_td(mod(itpp_ext::matlab_range(m0, 30 + m0), 31));
    const ivec s1_m1 = s_td(mod(itpp_ext::matlab_range(m1, 30 + m1), 31));

    const ivec c0 = c_td(mod(itpp_ext::matlab_range(n_id_2, 30 + n_id_2), 31));
    const ivec c1 = c_td(mod(itpp_ext::matlab_range(n_id_2 + 3, 30 + n_id_2 + 3), 31));

    const ivec z1_m0 = z_td(mod(itpp_ext::matlab_range(0, 30) + mod(m0, 8), 31));
    const ivec z1_m1 = z_td(mod(itpp_ext::matlab_range(0, 30) + mod(m1, 8), 31));
    ivec s(62);
    if (slot_num == 0)
    {
        ivec t1 = elem_mult(s1_m1, c1, z1_m0);
        ivec t2 = elem_mult(s0_m0, c0);
        for (size_t i = 0; i < 31; i++)
        {
            s[i * 2 + 1] = t1[i];
            s[i * 2] = t2[i];
        }
    }
    else if (slot_num == 10)
    {
        ivec t1 = elem_mult(s0_m0, c1, z1_m1);
        ivec t2 = elem_mult(s1_m1, c0);
        for (size_t i = 0; i < 31; i++)
        {
            s[i * 2 + 1] = t1[i];
            s[i * 2] = t2[i];
        }
    }
    return s;
}

double angle(complex<double> d)
{
    return std::atan2(d.imag(), d.real());
}

ivec wrap(ivec x, double sm, double lg)
{
    ivec res(x.size());
    for (size_t i = 0; i < x.size(); i++)
    {
        res[i] = mod(x[i] - sm, lg - sm) + sm;
    }
    return res;
}

void sss_detect(Cell &peak, cvec capbuf, int thresh2_n_sigma, double fc, int sampling_carrier_twist, double tdd_flag)
{
    double fs_lte = 30720000.0;
    int32 peak_loc = peak.ind;
    double peak_freq = peak.freq;
    int8 n_id_2_est = peak.n_id_2;

    double k_factor;
    if (sampling_carrier_twist == 1)
    {
        k_factor = (fc - peak.freq) / fc;
    }
    else
    {
        k_factor = peak.k_factor;
    }
    int min_idx, sss_ext_offset, sss_nrm_offset;

    if (tdd_flag == 1)
    {
        min_idx = 3 * (128 + 32) + 32;
        sss_ext_offset = 3 * (128 + 32);
        sss_nrm_offset = 412;
    }
    else
    {
        min_idx = 163 - 9;
        sss_ext_offset = 128 + 32;
        sss_nrm_offset = 128 + 9;
    }
    if (peak_loc < min_idx)
    {
        peak_loc = peak_loc + 9600 * k_factor;
    }
    vec pss_loc_set = itpp_ext::matlab_range((double)peak_loc, 9600 * k_factor, (double)(length(capbuf) - 125 - 9));
    int n_pss = pss_loc_set.size();
    vec pss_np(n_pss);
    cmat h_raw(n_pss, 62);
    cmat h_sm(n_pss, 62);
    cmat sss_nrm_raw(n_pss, 62);
    cmat sss_ext_raw(n_pss, 62);

    for (size_t k = 0; k < n_pss; k++)
    {
        double pss_loc = std::round(pss_loc_set(k));
        double pss_dft_location = pss_loc + 9 - 2;
        cvec dft_in = fshift(capbuf(pss_dft_location, pss_dft_location + 127), -peak_freq, fs_lte / 16);
        cvec dft_out = dft(concat(dft_in(2, -1), dft_in(0, 1)));
        h_raw.set_row(k, concat(dft_out(dft_out.length() - 31, -1), dft_out(1, 31)));
        h_raw.set_row(k, elem_mult(h_raw.get_row(k), conj(pss(n_id_2_est))));

        for (int t = 0; t < 62; t++)
        {
            int lt = max(0, t - 6);
            int rt = min(61, t + 6);
            h_sm(k, t) = mean(h_raw.get_row(k)(lt, rt));
        }

        pss_np(k) = sigpower(h_sm.get_row(k) - h_raw.get_row(k));

        double sss_ext_dft_location = pss_dft_location - sss_ext_offset;
        dft_in = fshift(capbuf(sss_ext_dft_location, sss_ext_dft_location + 127), -peak_freq, fs_lte / 16);
        dft_out = dft(concat(dft_in(2, -1), dft_in(0, 1)));
        sss_ext_raw.set_row(k, concat(dft_out(dft_out.length() - 31, -1), dft_out(1, 31)));

        double sss_nrm_dft_location = pss_dft_location - sss_nrm_offset;
        dft_in = fshift(capbuf(sss_nrm_dft_location, sss_nrm_dft_location + 127), -peak_freq, fs_lte / 16);
        dft_out = dft(concat(dft_in(2, -1), dft_in(0, 1)));
        sss_nrm_raw.set_row(k, concat(dft_out(dft_out.length() - 31, -1), dft_out(1, 31)));
    }

    cmat h_sm_ext_interp(h_sm);
    cmat h_sm_nrm_interp(h_sm);
    vec pss_np_ext(pss_np);
    vec pss_np_nrm(pss_np);

    vec sss_h1_np_est(62);
    vec sss_h2_np_est(62);
    // vec sss_h1_ext_np_est(62);
    // vec sss_h2_ext_np_est(62);

    cvec sss_h1_nrm_est(62);
    cvec sss_h2_nrm_est(62);
    cvec sss_h1_ext_est(62);
    cvec sss_h2_ext_est(62);
    vec pss_np_inv_h1 = 1.0 / pss_np(itpp_ext::matlab_range(0, 2, n_pss - 1));
    vec pss_np_inv_h2 = 1.0 / pss_np(itpp_ext::matlab_range(1, 2, n_pss - 1));
    for (int t = 0; t < 62; t++)
    {
        // First half (h1) and second half (h2) channel estimates.
        cvec h_sm_h1 = h_sm.get_col(t).get(itpp_ext::matlab_range(0, 2, n_pss - 1));
        cvec h_sm_h2 = h_sm.get_col(t).get(itpp_ext::matlab_range(1, 2, n_pss - 1));
        // Estimate noise power in each subcarrier
        sss_h1_np_est(t) = 1 / (1 + sum(elem_mult(sqr(h_sm_h1), pss_np_inv_h1)));
        sss_h2_np_est(t) = 1 / (1 + sum(elem_mult(sqr(h_sm_h2), pss_np_inv_h2)));
        // Estimate SSS assuming normal CP
        sss_h1_nrm_est(t) = sss_h1_np_est(t) * sum(elem_mult(conj(h_sm_h1), to_cvec(pss_np_inv_h1), sss_nrm_raw.get_col(t).get(itpp_ext::matlab_range(0, 2, n_pss - 1))));
        sss_h2_nrm_est(t) = sss_h2_np_est(t) * sum(elem_mult(conj(h_sm_h2), to_cvec(pss_np_inv_h2), sss_nrm_raw.get_col(t).get(itpp_ext::matlab_range(1, 2, n_pss - 1))));
        // Estimate SSS assuming extended CP
        sss_h1_ext_est(t) = sss_h1_np_est(t) * sum(elem_mult(conj(h_sm_h1), to_cvec(pss_np_inv_h1), sss_ext_raw.get_col(t).get(itpp_ext::matlab_range(0, 2, n_pss - 1))));
        sss_h2_ext_est(t) = sss_h2_np_est(t) * sum(elem_mult(conj(h_sm_h2), to_cvec(pss_np_inv_h2), sss_ext_raw.get_col(t).get(itpp_ext::matlab_range(1, 2, n_pss - 1))));
    }

    mat log_lik_nrm(168, 2);
    mat log_lik_ext(168, 2);
    for (int t = 0; t < 168; t++)
    {
        ivec sss_h1_try_i = sss(t, n_id_2_est, 0);
        ivec sss_h2_try_i = sss(t, n_id_2_est, 10);

        // Rotate the candiate sequence to match the received sequence.
        double ang = angle(sum(elem_mult(conj(concat(sss_h1_nrm_est, sss_h2_nrm_est)), to_cvec(concat(sss_h1_try_i, sss_h2_try_i)))));
        cvec sss_h1_try = sss_h1_try_i * exp(complex<double>(0, 1) * -ang);
        cvec sss_h2_try = sss_h2_try_i * exp(complex<double>(0, 1) * -ang);
        cvec df = concat(sss_h1_try, sss_h2_try) - concat(sss_h1_nrm_est, sss_h2_nrm_est);
        vec tmp = concat(real(df), imag(df));
        log_lik_nrm.set(t, 0, sum(elem_div(-elem_mult(tmp, tmp), repmat(concat(sss_h1_np_est, sss_h2_np_est), 2))));

        // Exchange h1 and h2 and re-do
        std::swap(sss_h1_try, sss_h2_try);
        ang = angle(sum(elem_mult(conj(concat(sss_h1_nrm_est, sss_h2_nrm_est)), concat(sss_h1_try, sss_h2_try))));
        sss_h1_try = sss_h1_try * exp(complex<double>(0, 1) * -ang);
        sss_h2_try = sss_h2_try * exp(complex<double>(0, 1) * -ang);
        df = concat(sss_h1_try, sss_h2_try) - concat(sss_h1_nrm_est, sss_h2_nrm_est);
        tmp = concat(real(df), imag(df));
        log_lik_nrm.set(t, 1, sum(elem_div(-elem_mult(tmp, tmp), repmat(concat(sss_h1_np_est, sss_h2_np_est), 2))));

        // Re-do for extended prefix
        // Rotate the candiate sequence to match the received sequence.
        std::swap(sss_h1_try, sss_h2_try);
        ang = angle(sum(elem_mult(conj(concat(sss_h1_ext_est, sss_h2_ext_est)), to_cvec(concat(sss_h1_try, sss_h2_try)))));
        sss_h1_try = sss_h1_try * exp(complex<double>(0, 1) * -ang);
        sss_h2_try = sss_h2_try * exp(complex<double>(0, 1) * -ang);
        df = concat(sss_h1_try, sss_h2_try) - concat(sss_h1_ext_est, sss_h2_ext_est);
        tmp = concat(real(df), imag(df));
        log_lik_ext.set(t, 0, sum(elem_div(-elem_mult(tmp, tmp), repmat(concat(sss_h1_np_est, sss_h2_np_est), 2))));

        std::swap(sss_h1_try, sss_h2_try);
        ang = angle(sum(elem_mult(conj(concat(sss_h1_ext_est, sss_h2_ext_est)), to_cvec(concat(sss_h1_try, sss_h2_try)))));
        sss_h1_try = sss_h1_try * exp(complex<double>(0, 1) * -ang);
        sss_h2_try = sss_h2_try * exp(complex<double>(0, 1) * -ang);
        df = concat(sss_h1_try, sss_h2_try) - concat(sss_h1_ext_est, sss_h2_ext_est);
        tmp = concat(real(df), imag(df));
        log_lik_ext.set(t, 1, sum(elem_div(-elem_mult(tmp, tmp), repmat(concat(sss_h1_np_est, sss_h2_np_est), 2))));
    }
    cp_type_t::cp_type_t cp_type_est;
    int cp_type_flag;
    mat log_lik;
    if (max(max(log_lik_nrm.get_col(0)), max(log_lik_nrm.get_col(1))) > max(max(log_lik_ext.get_col(0)), max(log_lik_ext.get_col(1))))
    {
        cp_type_est = cp_type_t::NORMAL;
        cp_type_flag = 0;
        log_lik = log_lik_nrm;
    }
    else
    {
        cp_type_est = cp_type_t::EXTENDED;
        cp_type_flag = 1;
        log_lik = log_lik_ext;
    }
    double frame_start;
    if (tdd_flag == 1)
    {
        if (cp_type_flag == 0)
        {
            frame_start = peak_loc + (-(2 * (128 + 9) + 1) - 1920 - 2) * k_factor;
        }
        else
        {
            frame_start = peak_loc + (-(2 * (128 + 32)) - 1920 - 2) * k_factor;
        }
    }
    else
    {
        frame_start = peak_loc + (128 + 9 - 960 - 2) * k_factor;
    }
    frame_start++;
    vec ll;
    if (max(log_lik.get_col(0)) > max(log_lik.get_col(1)))
    {
        ll = log_lik.get_col(0);
    }
    else
    {
        frame_start = frame_start + 9600 * k_factor;
        ll = log_lik.get_col(1);
    }
    // todo:
    frame_start = fmod(frame_start - 0.5, 2 * 9600) + 0.5;
    int n_id_1_est;
    double lik_final = max(ll, n_id_1_est);

    mat mL = concat_horizontal(log_lik_nrm, log_lik_ext);
    vec L(mL.rows() * mL.cols());
    for (int i = 0; i < mL.rows(); ++i)
    {
        for (int j = 0; j < mL.cols(); ++j)
        {
            L(i * mL.cols() + j) = mL(i, j);
        }
    }
    double L_mean = mean(L);
    double L_var = variance(L);
    if (lik_final < L_mean + sqrt(L_var) * thresh2_n_sigma)
    {
        peak.n_id_1 = -1;
        peak.cp_type = cp_type_t::UNKNOWN;
        // peak.cp_type_val = -1;
        peak.frame_start = -1;
        peak.duplex_mode = -1;
    }
    else
    {
        peak.n_id_1 = n_id_1_est;
        peak.cp_type = cp_type_est;
        // peak.cp_type_val = cp_type_flag;
        peak.frame_start = frame_start;
        peak.duplex_mode = tdd_flag;
    }
    return;
}

double refine_fo(cvec capbuf, cp_type_t::cp_type_t cp_type, int n_id_2, double freq, double fs, double frame_start, vec k_factor_vec, int &k_factor_idx)
{
	const int len_pss = 137;

    vec fo_set = "-3e3:2e3:3e3";
    fo_set = freq + fo_set;
	const int len = capbuf.length();
    vec corr_val(4);
    cmat td_pss_v(137, 3);
    pss_gen(td_pss_v);
    cvec td_pss = td_pss_v.get_col(n_id_2).transpose().get_row(0);
    for (size_t i = 0; i < 4; i++)
    {
	    const double k_factor_tmp = k_factor_vec(i);
        double pss_from_frame_start;
        if (cp_type_t::EXTENDED == cp_type)
        {
            pss_from_frame_start = k_factor_tmp * (1920 + 2 * (128 + 32));
        }
        else
        {
            pss_from_frame_start = k_factor_tmp * (1920 + 2 * (128 + 9) + 1);
        }
        double pss_sp = frame_start + pss_from_frame_start + 3;
        cvec tmp = complex<double>(0, 1) * itpp_ext::matlab_range(0, len_pss - 1);
        cvec pss_fo = conj(tmp);
        corr_val[i] = 0;
        int pss_count = 0;

        while (pss_sp + len_pss + 1 <= len)
        {
	        const double pss_idx = std::round(pss_sp);
            cvec chn_tmp = capbuf(pss_idx, pss_idx + len_pss - 1);
            double corr_tmp = pow(abs(sum(elem_mult(chn_tmp, pss_fo))), 2);
            corr_val(i) = corr_val(i) + corr_tmp;

            chn_tmp = capbuf(pss_idx + 1, pss_idx + 1 + len_pss - 1);
            corr_tmp = pow(abs(sum(elem_mult(chn_tmp, pss_fo))), 2);
            corr_val(i) = corr_val(i) + corr_tmp;

            chn_tmp = capbuf(pss_idx - 1, pss_idx - 1 + len_pss - 1);
            corr_tmp = pow(abs(sum(elem_mult(chn_tmp, pss_fo))), 2);
            corr_val(i) = corr_val(i) + corr_tmp;

            pss_count = pss_count + 1;
            pss_sp = pss_sp + k_factor_tmp * 5 * 1920;
        }

        corr_val(i) = corr_val(i) / pss_count;
    }
    max(corr_val, k_factor_idx);
    return fo_set(k_factor_idx);
}

void pss_sss_foe(Cell &peak, cvec capbuf, double fc, int sampling_carrier_twist, double tdd_flag)
{
    double k_factor;
    vec k_factor_vec(4);
    if (sampling_carrier_twist == 1)
    {
        k_factor = (fc - peak.freq) / fc;
        k_factor_vec = "3e3 : -2e3 : -3e3";
        k_factor_vec = (fc - peak.freq + k_factor_vec) / fc;
    }
    else
    {
        k_factor = peak.k_factor;
        k_factor_vec.ones();
        k_factor_vec = k_factor * k_factor_vec;
    }

    if (tdd_flag == 1)
    {
        int k_factor_idx;
        peak.freq = refine_fo(capbuf, peak.cp_type, peak.n_id_2, peak.freq, 1920000, peak.frame_start, k_factor_vec, k_factor_idx);
        k_factor = k_factor_vec(k_factor_idx);
    }
    double pss_sss_dist, first_sss_dft_location;
    if (peak.cp_type == cp_type_t::NORMAL)
    {
        if (tdd_flag == 1)
        {
            pss_sss_dist = std::round((3 * (128 + 9) + 1) * k_factor);
            first_sss_dft_location = peak.frame_start + (1920 - 128) * k_factor;
        }
        else
        {
            pss_sss_dist = std::round((128 + 9) * k_factor);
            first_sss_dft_location = peak.frame_start + (960 - 128 - 9 - 128) * k_factor;
        }
    }
    else if (peak.cp_type == cp_type_t::EXTENDED)
    {
        if (tdd_flag == 1)
        {
            pss_sss_dist = std::round(3 * (128 + 32) * k_factor);
            first_sss_dft_location = peak.frame_start + (1920 - 128) * k_factor;
        }
        else
        {
            pss_sss_dist = std::round((128 + 32) * k_factor);
            first_sss_dft_location = peak.frame_start + (960 - 128 - 32 - 128) * k_factor;
        }
    }

    first_sss_dft_location = fmod(first_sss_dft_location - 0.5, 9600 * 2) + 0.5;
    int sn;
    if (first_sss_dft_location - 9600 * k_factor > 0.5)
    {
        first_sss_dft_location = first_sss_dft_location - 9600 * k_factor;
        sn = 10;
    }
    else
    {
        sn = 0;
    }

    vec sss_dft_loc_set = itpp_ext::matlab_range(first_sss_dft_location, 9600 * k_factor, length(capbuf) - 127 - pss_sss_dist - 100);
    int n_sss = sss_dft_loc_set.length();
    cmat h_raw_fo_pss = zeros_c(n_sss, 62);
    cmat sss_raw_fo = zeros_c(n_sss, 62);
    cmat h_sm = zeros_c(n_sss, 62);
    cvec pss_np(n_sss);
    sn = (1 - sn / 10) * 10;
    complex<double> M(0, 0);
    for (int k = 0; k < n_sss; k++)
    {
        sn = (1 - sn / 10) * 10;
        double sss_dft_location = std::round(sss_dft_loc_set(k));

        double pss_dft_location = sss_dft_location + pss_sss_dist;
        cvec dft_in_pss = fshift(capbuf(pss_dft_location - 1, pss_dft_location + 127 - 1), -peak.freq, 30720000.0 / 16);
        cvec dft_in_pss_t = concat(dft_in_pss(2, -1), dft_in_pss(0, 1));
        cvec dft_out_pss = dft(dft_in_pss_t);
        cvec t = concat(dft_out_pss(dft_out_pss.length() - 31, -1), dft_out_pss(1, 31));
        h_raw_fo_pss.set_row(k, elem_mult(t, conj(pss(peak.n_id_2))));
        for (int t = 0; t < 62; t++)
        {
            int lt = max(0, t - 6);
            int rt = min(61, t + 6);
            h_sm(k, t) = mean(h_raw_fo_pss.get_row(k)(lt, rt));
        }

        pss_np(k) = sigpower(h_sm.get_row(k) - h_raw_fo_pss.get_row(k));
        cvec dft_in_sss = fshift(capbuf(sss_dft_location - 1, sss_dft_location + 127 - 1), -peak.freq, 30720000.0 / 16) * exp(complex<double>(0, 1) * pi * -peak.freq / (30720000.0 / 16 / 2) * -pss_sss_dist);

        cvec dft_in_sss_t = concat(dft_in_sss(2, -1), dft_in_sss(0, 1));
        cvec dft_out_sss = dft(dft_in_sss_t);

        t = concat(dft_out_sss(dft_out_sss.length() - 31, -1), dft_out_sss(1, 31));
        sss_raw_fo.set_row(k, elem_mult(t, conj(to_cvec(sss(peak.n_id_1, peak.n_id_2, sn)))));

        cvec tmp = elem_mult(conj(sss_raw_fo.get_row(k)), h_raw_fo_pss.get_row(k), to_cvec(absx2(h_sm.get_row(k))));


		vec h_sm_2 = 2 * absx2(h_sm.get_row(k));
		for (int i = 0; i < h_sm_2.size(); ++i) {
			complex<double> denominator = h_sm_2(i) * pss_np(k) + pow(pss_np(k), 2);
			complex<double> tmp_div_denominator = tmp[i] / denominator;
			M += tmp_div_denominator;
		}
        //M += sum(elem_div(tmp, 2 * absx2(h_sm.get_row(k)) * pss_np(k) + pow(pss_np(k), 2)));
    }

    double f_off_est = peak.freq + angle(M) / (2 * pi) / (1 / (30720000.0 / 16) * pss_sss_dist);
    peak.freq_fine = f_off_est;
}

inline cvec extract_psss(
    const cvec td_samps,
    const double foc_freq,
    const double &k_factor,
    const double &fs_programmed)
{
    // Frequency shift
    cvec dft_in = fshift(td_samps, foc_freq, fs_programmed * k_factor);
    // Remove the 2 sample time offset
    dft_in = concat(dft_in(2, -1), dft_in.left(2));
    // DFT
    const cvec dft_out = dft(dft_in);
    // Extract interesting samples.
    return concat(dft_out.right(31), dft_out.mid(1, 31));
}

// void pss_sss_foe(
//   Cell & cell_in,
//   cvec capbuf,
//   double fc,
//   int sampling_carrier_twist,
//   double tdd_flag
// ) {

//   vec k_factor_vec(4);
//   double k_factor;
//   if (sampling_carrier_twist){
// //    k_factor=(fc_requested-cell_in.freq)/fc_programmed;
//     k_factor=(fc-cell_in.freq)/fc;
//     k_factor_vec(0) = (fc-cell_in.freq+3e3)/fc;
//     k_factor_vec(1) = (fc-cell_in.freq+1e3)/fc;
//     k_factor_vec(2) = (fc-cell_in.freq-1e3)/fc;
//     k_factor_vec(3) = (fc-cell_in.freq-3e3)/fc;
//   } else {
//     k_factor = cell_in.k_factor;
//     k_factor_vec(0) = k_factor;
//     k_factor_vec(1) = k_factor;
//     k_factor_vec(2) = k_factor;
//     k_factor_vec(3) = k_factor;
//   }

//   int k_factor_idx;
//   if (tdd_flag == 1) {
//     cell_in.freq = refine_fo(capbuf, cell_in.cp_type, cell_in.n_id_2, cell_in.freq, fc, cell_in.frame_start, k_factor_vec, k_factor_idx);
//     k_factor = k_factor_vec(k_factor_idx);
//   }

//   // Determine where we can find both PSS and SSS
//   uint16 pss_sss_dist;
//   double first_sss_dft_location;
//   if (cell_in.cp_type==cp_type_t::NORMAL) {
//     if (tdd_flag==0)
//     {
//         pss_sss_dist=itpp::round_i((128+9)*16/FS_LTE*fc*k_factor); //FDD
//         first_sss_dft_location=cell_in.frame_start+(960-128-9-128)*16/FS_LTE*fc*k_factor; //FDD
//     }
//     else
//     {
//         pss_sss_dist=itpp::round_i((3*(128+9)+1)*16/FS_LTE*fc*k_factor); //TDD
//         first_sss_dft_location=cell_in.frame_start+(1920-128)*16/FS_LTE*fc*k_factor; //TDD
//     }
//   } else if (cell_in.cp_type==cp_type_t::EXTENDED) {
//     if (tdd_flag==0)
//     {
//         pss_sss_dist=round_i((128+32)*16/FS_LTE*fc*k_factor); //FDD
//         first_sss_dft_location=cell_in.frame_start+(960-128-32-128)*16/FS_LTE*fc*k_factor; //FDD
//     }
//     else
//     {
//         pss_sss_dist=round_i((3*(128+32))*16/FS_LTE*fc*k_factor); //TDD
//         first_sss_dft_location=cell_in.frame_start+(1920-128)*16/FS_LTE*fc*k_factor; //TDD
//     }
//   } else {
//     throw("Error... check code...");
//   }
//   uint8 sn;
//   first_sss_dft_location=fmod(first_sss_dft_location -0.5, 9600*2) +0.5;
//   if (first_sss_dft_location-9600*k_factor>0.5) {
//     first_sss_dft_location-=9600*k_factor;
//     sn=10;
//   } else {
//     sn=0;
//   }
//   vec sss_dft_loc_set=itpp_ext::matlab_range(first_sss_dft_location,9600*16/FS_LTE*fc*k_factor,(double)(length(capbuf)-127-pss_sss_dist-100));
//   uint16 n_sss=length(sss_dft_loc_set);

//   // Loop around for each PSS/SSS pair
//   sn=(1-(sn/10))*10;
//   complex <double> M(0,0);
//   cmat h_raw_fo_pss(n_sss,62);
//   cmat h_sm(n_sss,62);
//   cmat sss_raw_fo(n_sss,62);
//   vec pss_np(n_sss);
// #ifndef NDEBUG
//   h_raw_fo_pss=NAN;
//   h_sm=NAN;
//   sss_raw_fo=NAN;
//   pss_np=NAN;
// #endif
//   for (uint16 k=0;k<n_sss;k++) {
//     sn=(1-(sn/10))*10;
//     uint32 sss_dft_location=round_i(sss_dft_loc_set(k));

//     // Find the PSS and use it to estimate the channel.
//     uint32 pss_dft_location=sss_dft_location+pss_sss_dist;
//     h_raw_fo_pss.set_row(k,extract_psss(capbuf.mid(pss_dft_location,128),-cell_in.freq,k_factor,fc));
//     h_raw_fo_pss.set_row(k,elem_mult(h_raw_fo_pss.get_row(k),conj(ROM_TABLES.pss_fd[cell_in.n_id_2])));

//     // Smoothing... Basic...
//     for (uint8 t=0;t<62;t++) {
//       uint8 lt=max(0,t-6);
//       uint8 rt=min(61,t+6);
//       h_sm(k,t)=mean(h_raw_fo_pss.get_row(k).mid(lt,rt-lt+1));
//     }

//     // Estimate noise power.
//     pss_np(k)=sigpower(h_sm.get_row(k)-h_raw_fo_pss.get_row(k));

//     // Calculate the SSS in the frequency domain
//     sss_raw_fo.set_row(k,extract_psss(capbuf.mid(sss_dft_location,128),-cell_in.freq,k_factor,fc)*exp(complex<double>(0,1)*pi*-cell_in.freq/(FS_LTE/16/2)*-pss_sss_dist));
//     sss_raw_fo.set_row(k,elem_mult(sss_raw_fo.get_row(k),to_cvec(ROM_TABLES.sss_fd(cell_in.n_id_1,cell_in.n_id_2,sn))));

//     // Compare PSS to SSS. With no frequency offset, arg(M) is zero.
//     M=M+sum(elem_mult(
//       conj(sss_raw_fo.get_row(k)),
//       h_raw_fo_pss.get_row(k),
//       to_cvec(elem_mult(
//         sqr(h_sm.get_row(k)),
//         1.0/(2*sqr(h_sm.get_row(k))*pss_np(k)+sqr(pss_np(k)))
//       ))
//     ));
//   }

//   // Store results.
//   cell_in.freq_fine=cell_in.freq+arg(M)/(2*pi)/(1/(fc*k_factor)*pss_sss_dist);
// }

int conv_idx(double in_idx, int decimation_ratio)
{
    return in_idx * decimation_ratio;
}

void extract_tfg(Cell &peak, cvec capbuf, double fc, int sampling_carrier_twist, int nRB, cmat &tfg, vec &tfg_timestamp)
{
	const double frame_start = peak.frame_start;
	const cp_type_t::cp_type_t cp_type = peak.cp_type;
	const double freq_fine = peak.freq_fine;

    int conv_idx_ratio, decimation_ratio;
    if (nRB == 6)
    {
        conv_idx_ratio = 1;
        decimation_ratio = 16;
    }
    else
    {
        conv_idx_ratio = 16;
        decimation_ratio = 1;
    }

	const double fs = 30720000.0 / decimation_ratio;
	const double fft_size = 2048 / decimation_ratio;
	const double len_cp_extended = 512 / decimation_ratio;
	const double len_cp_normal1 = 144 / decimation_ratio;
	const double len_cp_normal2 = 160 / decimation_ratio;
	const int nSC = nRB * 12;

    double k_factor;
    if (sampling_carrier_twist == 1)
    {
        k_factor = (fc - freq_fine) / fc;
    }
    else
    {
        k_factor = peak.k_factor;
    }
    int n_symb_dl;
    double dft_location;
    if (cp_type == cp_type_t::NORMAL)
    {
        n_symb_dl = 7;
        dft_location = conv_idx(frame_start + len_cp_normal2 / conv_idx_ratio, conv_idx_ratio);
    }
    else if (cp_type == cp_type_t::EXTENDED)
    {
        n_symb_dl = 6;
        dft_location = conv_idx(frame_start + len_cp_extended / conv_idx_ratio, conv_idx_ratio);
    }
    // peak.n_symb_dl = n_symb_dl;
    if (dft_location - k_factor * fs * .01 >= 0.5)
    {
        dft_location = dft_location - k_factor * fs * .01;
    }
    capbuf = fshift(capbuf, -freq_fine, fs);
	const int n_ofdm_sym = 6 * 10 * 2 * n_symb_dl + 2 * n_symb_dl;
    tfg = cmat(n_ofdm_sym, fft_size);
    cmat tfg_t = cmat(n_ofdm_sym, fft_size);
    tfg_timestamp = vec(n_ofdm_sym);
    int sym_num = 0;
    for (int t = 0; t < n_ofdm_sym; t++)
    {
        ivec indices = to_ivec(itpp_ext::matlab_range(std::round(dft_location)-1, std::round(dft_location) + fft_size - 2));
        tfg_t.set_row(t, dft(capbuf(indices)));
        tfg_timestamp(t) = dft_location;
        if (n_symb_dl == 6)
        {
            dft_location = dft_location + k_factor * (fft_size + len_cp_extended);
        }
        else
        {
            if (sym_num == 6)
            {
                dft_location = dft_location + k_factor * (fft_size + len_cp_normal2);
            }
            else
            {
                dft_location = dft_location + k_factor * (fft_size + len_cp_normal1);
            }
            sym_num = mod(sym_num + 1, 7);
        }
    }
    tfg = concat_horizontal(tfg_t.get_cols(tfg_t.cols() -1 - (nSC / 2 - 1), tfg_t.cols() -1), tfg_t.get_cols(1, nSC / 2));
	const ivec cn = concat(itpp_ext::matlab_range(-(nSC / 2), -1), itpp_ext::matlab_range(1, nSC / 2));
    for (int t = 0; t < n_ofdm_sym; t++)
    {
	    const double ideal_offset = tfg_timestamp(t);
	    const double actual_offset = std::round(ideal_offset);
        double late = actual_offset - ideal_offset;
        tfg.set_row(t, elem_mult(tfg.get_row(t), exp(complex<double>(0, -1) * 2 * pi * cn * late / fft_size)));
    }
}

cvec rs_dl(int slot_num, int sym_num, int port_num, int n_id_cell, int n_rb_dl, cp_type_t::cp_type_t cp_type, int &shift)
{
	const int n_rb_maxdl = 110;
    int n_symb_dl, n_cp;
    double dft_location;

    cvec r(1);
    if (cp_type == cp_type_t::NORMAL)
    {
        n_symb_dl = 7;
        n_cp = 1;
    }
    else if (cp_type == cp_type_t::EXTENDED)
    {
        n_symb_dl = 6;
        n_cp = 0;
    }
    if (port_num == 0 || port_num == 1)
    {
        if (sym_num != 0 && sym_num != n_symb_dl - 3)
        {
            return r; // Exit the function or return statement in C++
        }
    }

    if (port_num == 2 || port_num == 3)
    {
        if (sym_num != 1)
        {
            return r; // Exit the function or return statement in C++
        }
    }

	const double c_init = 2 ^ 10 * (7 * (slot_num + 1) + sym_num + 1) * (2 * n_id_cell + 1) + 2 * n_id_cell + n_cp;
    vec c = to_vec(lte_pn(c_init, 4 * n_rb_maxdl));
    cvec r_l_ns(c.length() / 2);
    for (int t = 0; t < c.length() / 2; t++)
    {
        r_l_ns[t] = 1 / sqrt(2) * complex<double>(1 - 2 * c[2 * t], 1 - 2 * c[2 * t + 1]);
    }
    r = r_l_ns(1 + n_rb_maxdl - n_rb_dl, 2 * n_rb_dl + n_rb_maxdl - n_rb_dl);

    int v = 0;
    int v_shift;

    if (port_num == 0 && sym_num == 0)
    {
        v = 0;
    }
    else if (port_num == 0 && sym_num != 0)
    {
        v = 3;
    }
    else if (port_num == 1 && sym_num == 0)
    {
        v = 3;
    }
    else if (port_num == 1 && sym_num != 0)
    {
        v = 0;
    }
    else if (port_num == 2)
    {
        v = 3 * (slot_num % 2);
    }
    else if (port_num == 3)
    {
        v = 3 + 3 * (slot_num % 2);
    }

    v_shift = n_id_cell % 6;
    shift = (v + v_shift) % 6;
    return r;
}

Cell tfoec(
    // Inputs
    const Cell &cell,
    const cmat &tfg,
    const vec &tfg_timestamp,
    const double &fc,
    int sampling_carrier_twist,
    const RS_DL &rs_dl,
    // Outputs
    cmat &tfg_comp,
    vec &tfg_comp_timestamp)
{
    // cout << cell <<endl;
    // cout << "=========================" <<endl;
    // cout << tfg <<endl;
    // cout << "=========================" <<endl;
    // cout << tfg_timestamp <<endl;
    // cout << "=========================" <<endl;

    // Local shortcuts
    const int8 n_symb_dl = cell.n_symb_dl();
    uint16 n_ofdm = tfg.rows();
    uint16 n_slot = floor((double)n_ofdm / n_symb_dl);

    // Perform super-fine FOE
    complex<double> foe;
    for (uint8 sym_num = 0; sym_num <= n_symb_dl - 3; sym_num += n_symb_dl - 3)
    {
        // Extract all the RS and compensate for the known transmitted symbols.
        cmat rs_extracted(n_slot, 12);
#ifndef NDEBUG
        rs_extracted = NAN;
#endif
        for (uint16 t = 0; t < n_slot; t++)
        {
            // Extract RS
            rs_extracted.set_row(t, tfg.get_row(t * n_symb_dl + sym_num).get(itpp_ext::matlab_range((uint32)rs_dl.get_shift(mod(t, 20), sym_num, 0), (uint32)6, (uint32)71)));
            // Compensate
            rs_extracted.set_row(t, elem_mult(rs_extracted.get_row(t), conj(rs_dl.get_rs(mod(t, 20), sym_num))));
        }
        // FOE, subcarrier by subcarrier.
        for (uint16 t = 0; t < 12; t++)
        {
            cvec col = rs_extracted.get_col(t);
            foe = foe + sum(elem_mult(conj(col(0, n_slot - 2)), col(1, -1)));
        }
    }
    double k_factor;
    if (sampling_carrier_twist)
    {
        k_factor = (fc - cell.freq_fine) / fc;
    }
    else
    {
        k_factor = cell.k_factor;
    }
    double residual_f = arg(foe) / (2 * pi) / (k_factor * 0.0005);

    // Perform FOC. Does not fix ICI!
    double k_factor_residual;
    if (sampling_carrier_twist)
    {
        k_factor_residual = (fc - residual_f) / fc;
    }
    else
    {
        //    k_factor_residual = cell.k_factor;
        k_factor_residual = 1.0;
    }

    tfg_comp = cmat(n_ofdm, 72);
#ifndef NDEBUG
    tfg_comp = NAN;
#endif
    tfg_comp_timestamp = k_factor_residual * tfg_timestamp;
    ivec cn = concat(itpp_ext::matlab_range(-36, -1), itpp_ext::matlab_range(1, 36));
    for (uint16 t = 0; t < n_ofdm; t++)
    {
        tfg_comp.set_row(t, tfg.get_row(t) * exp(complex<double>(0, 1) * 2 * pi * -residual_f * tfg_comp_timestamp(t) / (FS_LTE / 16)));
        // How late were we in locating the DFT
        double late = tfg_timestamp(t) - tfg_comp_timestamp(t);
        // Compensate for the improper location of the DFT
        tfg_comp.set_row(t, elem_mult(tfg_comp.get_row(t), exp(complex<double>(0, -1) * 2 * pi * late / 128 * cn)));
    }

    // Perform TOE.
    // Implemented by comparing subcarrier k of one OFDM symbol with subcarrier
    // k+3 of another OFDM symbol. This is why FOE must be performed first.
    // Slightly less performance but faster execution time could be obtained
    // by comparing subcarrier k with subcarrier k+6 of the same OFDM symbol.
    complex<double> toe = 0;
    for (uint16 t = 0; t < 2 * n_slot - 1; t++)
    {
        // Current OFDM symbol containing RS
        uint8 current_sym_num = t & 1 ? n_symb_dl - 3 : 0;
        uint8 current_slot_num = mod(t >> 1, 20);
        uint16 current_offset = (t >> 1) * n_symb_dl + current_sym_num;
        // Since we are using port 0, the shift is the same for all slots.
        uint8 current_shift = rs_dl.get_shift(0, current_sym_num, 0);
        // Next OFDM symbol containing RS
        uint8 next_sym_num = t + 1 & 1 ? n_symb_dl - 3 : 0;
        uint8 next_slot_num = mod(t + 1 >> 1, 20);
        uint16 next_offset = (t + 1 >> 1) * n_symb_dl + next_sym_num;
        // Since we are using port 0, the shift is the same for all slots.
        uint8 next_shift = rs_dl.get_shift(0, next_sym_num, 0);

        uint16 r1_offset, r2_offset;
        uint8 r1_shift, r2_shift;
        uint8 r1_sym_num, r2_sym_num;
        uint8 r1_slot_num, r2_slot_num;
        if (current_shift < next_shift)
        {
            r1_offset = current_offset;
            r1_shift = current_shift;
            r1_sym_num = current_sym_num;
            r1_slot_num = current_slot_num;
            r2_offset = next_offset;
            r2_shift = next_shift;
            r2_sym_num = next_sym_num;
            r2_slot_num = next_slot_num;
        }
        else
        {
            r1_offset = next_offset;
            r1_shift = next_shift;
            r1_sym_num = next_sym_num;
            r1_slot_num = next_slot_num;
            r2_offset = current_offset;
            r2_shift = current_shift;
            r2_sym_num = current_sym_num;
            r2_slot_num = current_slot_num;
        }
        cvec r1v = tfg_comp.get_row(r1_offset).get(itpp_ext::matlab_range(r1_shift, 6, 71));
        r1v = elem_mult(r1v, conj(rs_dl.get_rs(r1_slot_num, r1_sym_num)));
        cvec r2v = tfg_comp.get_row(r2_offset).get(itpp_ext::matlab_range(r2_shift, 6, 71));
        r2v = elem_mult(r2v, conj(rs_dl.get_rs(r2_slot_num, r2_sym_num)));
        complex<double> toe1 = sum(elem_mult(conj(r1v), r2v));
        complex<double> toe2 = sum(elem_mult(conj(r2v(0, 10)), r1v(1, 11)));
        toe += toe1 + toe2;
    }
    double delay = -arg(toe) / 3 / (2 * pi / 128);

    // Perform TOC
    cvec comp_vector = exp(complex<double>(0, 1) * 2 * pi / 128 * delay * cn);
    for (uint16 t = 0; t < n_ofdm; t++)
    {
        tfg_comp.set_row(t, elem_mult(tfg_comp.get_row(t), comp_vector));
    }

    Cell cell_out(cell);
    cell_out.freq_superfine = cell_out.freq_fine + residual_f;
    if (sampling_carrier_twist)
    {
        //    cell_out.k_factor=(fc_requested-cell_out.freq_superfine)/fc_programmed;
        cell_out.k_factor = (fc - cell_out.freq_superfine) / fc;
    }
    return cell_out;
}

void pbch_extract_2(
    // Inputs
    const Cell &cell,
    const cmat &tfg,
    const Array<cmat> &ce,
    // Outputs
    cvec &pbch_sym,
    cmat &pbch_ce)
{
    // Shortcuts
    const int8 n_symb_dl = cell.n_symb_dl();
    const uint16 m_bit = cell.cp_type == cp_type_t::NORMAL ? 1920 : 1728;
    const uint8 v_shift_m3 = mod(cell.n_id_cell(), 3);

    pbch_sym = cvec(m_bit / 2);
    // One channel estimate from each of 4 ports for each RE.
    pbch_ce = cmat(4, m_bit / 2);
#ifndef NDEBUG
    pbch_sym = NAN;
    pbch_ce = NAN;
#endif
    uint32 idx = 0;
    for (uint8 fr = 0; fr <= 3; fr++)
    {
        for (uint8 sym = 0; sym <= 3; sym++)
        {
            for (uint8 sc = 0; sc <= 71; sc++)
            {
                // Skip if there might be an RS occupying this position.
                if (mod(sc, 3) == v_shift_m3 && (sym == 0 || sym == 1 || (sym == 3 && n_symb_dl == 6)))
                {
                    continue;
                }
                const uint16 sym_num = fr * 10 * 2 * n_symb_dl + n_symb_dl + sym;
                pbch_sym(idx) = tfg(sym_num, sc);
                pbch_ce(0, idx) = ce(0).get(sym_num, sc);
                pbch_ce(1, idx) = ce(1).get(sym_num, sc);
                pbch_ce(2, idx) = ce(2).get(sym_num, sc);
                pbch_ce(3, idx) = ce(3).get(sym_num, sc);
                idx++;
            }
        }
    }
}

void del_oob_2(
    ivec &v)
{
    int32 t = 0;
    while (t < v.length())
    {
        if (v(t) < 0 || v(t) > 11)
        {
            v.del(t);
        }
        else
        {
            t++;
        }
    }
}

void ce_interp_hex_extend_2(
    vec &row_x,
    cvec &row_val)
{
    if (row_x(0) != 0)
    {
        row_val.ins(0, row_val(0) - row_x(0) * (row_val(1) - row_val(0)) / (row_x(1) - row_x(0)));
        row_x.ins(0, 0);
    }
    if (itpp_ext::last(row_x) != 71)
    {
	    const uint16 len = length(row_val);
        row_val.ins(len, row_val(len - 1) + (71 - itpp_ext::last(row_x)) * (row_val(len - 1) - row_val(len - 2)) / (row_x(len - 1) - row_x(len - 2)));
        row_x.ins(len, 71);
    }
}

struct triangle_vertex_t
{
    uint8 x_sc;
    uint16 y_symnum;
    complex<double> val;
};

void ce_interp_hex_2(
    // Inputs
    const cmat &ce_filt,
    const ivec &shift,
    const int16 &n_ofdm,
    const int16 &n_rs_ofdm,
    const ivec &rs_set,
    // Outputs
    cmat &ce_tfg)
{
    ce_tfg = cmat(n_ofdm, 72);
#ifndef NDEBUG
    ce_tfg = NAN;
#endif
    for (uint16 t = 0; t <= n_rs_ofdm - 2; t++)
    {
        // Extract two rows and ensure that there are samples for subcarriers
        // 0 and 71.
        // In general, top_row_* is actually equal to bot_row_* from the
        // previous iteration.
        vec top_row_x = to_vec(itpp_ext::matlab_range(t & 1 ? shift(1) : shift(0), 6, 71));
        cvec top_row_val = ce_filt.get_row(t);
        ce_interp_hex_extend_2(top_row_x, top_row_val);
        vec bot_row_x = to_vec(itpp_ext::matlab_range(t & 1 ? shift(0) : shift(1), 6, 71));
        cvec bot_row_val = ce_filt.get_row(t + 1);
        ce_interp_hex_extend_2(bot_row_x, bot_row_val);

        // First row is not handled inside the main loop.
        if (t == 0)
        {
            ce_tfg.set_row(rs_set(0), interp1(top_row_x, top_row_val, itpp_ext::matlab_range(0.0, 71.0)));
        }

        // Create initial triangle
        uint8 top_row_last_used;
        uint8 bot_row_last_used;
        vector<triangle_vertex_t> triangle(3);
        if (top_row_x(1) < bot_row_x(1))
        {
            triangle[0].x_sc = top_row_x(0);
            triangle[0].y_symnum = rs_set(t);
            triangle[0].val = top_row_val(0);
            triangle[1].x_sc = bot_row_x(0);
            triangle[1].y_symnum = rs_set(t + 1);
            triangle[1].val = bot_row_val(0);
            triangle[2].x_sc = top_row_x(1);
            triangle[2].y_symnum = rs_set(t);
            triangle[2].val = top_row_val(1);
            top_row_last_used = 1;
            bot_row_last_used = 0;
        }
        else
        {
            triangle[0].x_sc = bot_row_x(0);
            triangle[0].y_symnum = rs_set(t + 1);
            triangle[0].val = bot_row_val(0);
            triangle[1].x_sc = top_row_x(0);
            triangle[1].y_symnum = rs_set(t);
            triangle[1].val = top_row_val(0);
            triangle[2].x_sc = bot_row_x(1);
            triangle[2].y_symnum = rs_set(t + 1);
            triangle[2].val = bot_row_val(1);
            top_row_last_used = 0;
            bot_row_last_used = 1;
        }

        // This loop succesively creates triangles to cover the space between
        // top_row and bot_row.
        const uint8 spacing = rs_set(t + 1) - rs_set(t);
        vec x_offset(spacing + 1);
        x_offset = 0.0;
        while (true)
        {
            // Calculate the parameters of a plane passing through all points
            // of the triangle.
            // value=a_p*x_sc+b_p*y_symnum+c_p;
            cmat M(3, 3);
            M(0, 0) = triangle[0].x_sc;
            M(1, 0) = triangle[1].x_sc;
            M(2, 0) = triangle[2].x_sc;
            M(0, 1) = triangle[0].y_symnum;
            M(1, 1) = triangle[1].y_symnum;
            M(2, 1) = triangle[2].y_symnum;
            M(0, 2) = 1;
            M(1, 2) = 1;
            M(2, 2) = 1;
            cvec V(3);
            V(0) = triangle[0].val;
            V(1) = triangle[1].val;
            V(2) = triangle[2].val;
            // In the future, inv(M) can be calculated directly for speed
            // since the last column of M is all ones.
            cvec abc = inv(M) * V;
            complex<double> a_p = abc(0);
            complex<double> b_p = abc(1);
            complex<double> c_p = abc(2);

            // Calculate the parameters of the line defining the rightmost
            // edge of the triangle.
            // x_sc=a_l*y_symnum+b_l;
            const double x1 = triangle[1].x_sc;
            const double x2 = triangle[2].x_sc;
            const double y1 = triangle[1].y_symnum;
            const double y2 = triangle[2].y_symnum;
            const double a_l = (x1 - x2) / (y1 - y2);
            const double b_l = (y1 * x2 - y2 * x1) / (y1 - y2);

            for (uint8 r = 1; r <= spacing; r++)
            {
                while (x_offset(r) <= a_l * (rs_set(t) + r) + b_l)
                {
                    ce_tfg(rs_set(t) + r, x_offset(r)) = a_p * x_offset(r) + b_p * (rs_set(t) + r) + c_p;
                    x_offset(r)++;
                }
            }

            if (x_offset(1) == 72 && itpp_ext::last(x_offset) == 72)
            {
                break;
            }

            // We are not done yet. Choose the points for the next triangle.
            if (triangle[2].y_symnum == rs_set(t))
            {
                triangle[0] = triangle[1];
                triangle[1] = triangle[2];
                bot_row_last_used++;
                triangle[2].x_sc = bot_row_x(bot_row_last_used);
                triangle[2].y_symnum = rs_set(t + 1);
                triangle[2].val = bot_row_val(bot_row_last_used);
            }
            else
            {
                triangle[0] = triangle[1];
                triangle[1] = triangle[2];
                top_row_last_used++;
                triangle[2].x_sc = top_row_x(top_row_last_used);
                triangle[2].y_symnum = rs_set(t);
                triangle[2].val = top_row_val(top_row_last_used);
            }
        }
    }

    // Rows before the first and after the last OFDM symbol with RS are
    // created simply by copying the nearest OFDM symbol with RS.
    for (uint8 t = 0; t < rs_set(0); t++)
    {
        ce_tfg.set_row(t, ce_tfg.get_row(rs_set(0)));
    }
    for (uint16 t = itpp_ext::last(rs_set) + 1; t < n_ofdm; t++)
    {
        ce_tfg.set_row(t, ce_tfg.get_row(itpp_ext::last(rs_set)));
    }
}

void chan_est_2(
    // Inputs
    const Cell &cell,
    const RS_DL &rs_dl,
    const cmat &tfg,
    const uint8 &port,
    // Outputs
    cmat &ce_tfg,
    double &np)
{
    const int8 n_symb_dl = cell.n_symb_dl();
    const uint16 n_ofdm = tfg.rows();

    // Set of OFDM symbols containing reference symbols.
    ivec rs_set;
    if (port <= 1)
    {
        // There are better ways to implement this...
        Sort<int> sort;
        rs_set = concat(itpp_ext::matlab_range(0, n_symb_dl, n_ofdm - 1), itpp_ext::matlab_range(n_symb_dl - 3, n_symb_dl, n_ofdm - 1));
        sort.sort(0, rs_set.length() - 1, rs_set);
        // rs_set=reverse(rs_set);
    }
    else
    {
        rs_set = itpp_ext::matlab_range(1, n_symb_dl, n_ofdm - 1);
    }
    const uint16 n_rs_ofdm = length(rs_set);

    // Extract the raw channel estimates. 12 raw channel estimates per OFDM
    // symbol containing RS.
    cmat ce_raw(n_rs_ofdm, 12);
#ifndef NDEBUG
    ce_raw = NAN;
#endif
    uint8 slot_num = 0;
    ivec shift(2);
    shift = -1000;
    for (uint16 t = 0; t < n_rs_ofdm; t++)
    {
        uint8 sym_num = mod(rs_set(t), n_symb_dl);
        if (t <= 1)
        {
            shift(t) = rs_dl.get_shift(mod(slot_num, 20), sym_num, port);
        }

        cvec rs = rs_dl.get_rs(slot_num, sym_num);
        // Extract
        cvec raw_row = tfg.get_row(rs_set(t)).get(itpp_ext::matlab_range((int32)rs_dl.get_shift(mod(slot_num, 20), sym_num, port), 6, 71));
        ce_raw.set_row(t, raw_row);
        // Compensate for known RS
        ce_raw.set_row(t, elem_mult(ce_raw.get_row(t), conj(rs)));
        if ((t & 1) == 1 || port >= 2)
        {
            slot_num = mod(slot_num + 1, 20);
        }
    }

    // Simple filtering.
    //
    // 1   2   3   4   5   6
    //   7   8   9   A   B
    // C   D   E   F   G   H
    //
    // If the above is a representation of the location of the raw channel
    // estimates in the TFG, the filtered channel estimate for position 8
    // is the mean of the raw channel estimates for positions 2,3,7,8,9,D,
    // and E.
    cmat ce_filt(n_rs_ofdm, 12);
    bool current_row_leftmost = shift(0) < shift(1);
    for (uint16 t = 0; t < n_rs_ofdm; t++)
    {
        complex<double> total;
        uint8 n_total;
        for (uint16 k = 0; k < 12; k++)
        {
            // Current time offset
            ivec ind;
            ind = itpp_ext::matlab_range(k - 1, k + 1);
            del_oob_2(ind);
            total = sum(ce_raw.get_row(t).get(ind));
            n_total = length(ind);

            // Add in the previous and next time offset (if they exist)
            if (shift(0) == shift(1))
            {
                ind = itpp_ext::matlab_range(k - 1, k + 1);
            }
            else
            {
                if (current_row_leftmost)
                {
                    ind = itpp_ext::matlab_range(k - 1, k);
                }
                else
                {
                    ind = itpp_ext::matlab_range(k, k + 1);
                }
            }
            del_oob_2(ind);
            // Previous time offset
            if (t != 0)
            {
                total += sum(ce_raw.get_row(t - 1).get(ind));
                n_total += length(ind);
            }
            if (t != n_rs_ofdm - 1)
            {
                total += sum(ce_raw.get_row(t + 1).get(ind));
                n_total += length(ind);
            }
            ce_filt(t, k) = total / n_total;
        }
        current_row_leftmost = !current_row_leftmost;
    }

    // Estimate noise power
    np = sigpower(cvectorize(ce_filt) - cvectorize(ce_raw));

    // There is no appreciable difference in performance between these
    // algorithms for high SNR values.
    // ce_interp_2stage(ce_filt,shift,n_ofdm,n_rs_ofdm,rs_set,ce_tfg);
    // ce_interp_freq_time(ce_filt,shift,n_ofdm,n_rs_ofdm,rs_set,ce_tfg);
    ce_interp_hex_2(ce_filt, shift, n_ofdm, n_rs_ofdm, rs_set, ce_tfg);
}

Cell decode_mib_2(
    const Cell &cell,
    const cmat &tfg,
    const RS_DL &rs_dl)
{
    // Local shortcuts
    const int8 n_symb_dl = cell.n_symb_dl();

    Cell cell_out = cell;

    // Channel estimation. This is automatically performed for four antennas
    // and for every RE, not only the RE's that contain an MIB!!!
    Array<cmat> ce_tfg(4);
    vec np_v(4);
    chan_est_2(cell, rs_dl, tfg, 0, ce_tfg(0), np_v(0));
    chan_est_2(cell, rs_dl, tfg, 1, ce_tfg(1), np_v(1));
    chan_est_2(cell, rs_dl, tfg, 2, ce_tfg(2), np_v(2));
    chan_est_2(cell, rs_dl, tfg, 3, ce_tfg(3), np_v(3));

    // Try various frame offsets and number of TX antennas.
    bvec c_est;
    for (uint8 frame_timing_guess = 0; frame_timing_guess <= 3; frame_timing_guess++)
    {
        const uint16 ofdm_sym_set_start = frame_timing_guess * 10 * 2 * n_symb_dl;
        ivec ofdm_sym_set = itpp_ext::matlab_range(ofdm_sym_set_start, ofdm_sym_set_start + 3 * 10 * 2 * n_symb_dl + 2 * n_symb_dl - 1);

        // Extract only the portion of the TFG containing the four frames
        // we are interested in.
        cmat tfg_try = tfg.get_rows(ofdm_sym_set);
        Array<cmat> ce_try(4);
        for (uint8 t = 0; t < 4; t++)
        {
            ce_try(t) = ce_tfg(t).get_rows(ofdm_sym_set);
        }

        // Extract symbols and channel estimates for the PBCH
        cvec pbch_sym;
        cmat pbch_ce;
        pbch_extract_2(cell, tfg_try, ce_try, pbch_sym, pbch_ce);
        // Try 1, 2, and 4 ports.
        vec np;
        cvec syms;
        for (uint8 n_ports_pre = 1; n_ports_pre <= 3; n_ports_pre++)
        {
            const uint8 n_ports = n_ports_pre == 3 ? 4 : n_ports_pre;
            // Perform channel compensation and also estimate noise power in each
            // symbol.
            if (n_ports == 1)
            {
                cvec gain = conj(elem_div(pbch_ce.get_row(0), to_cvec(sqr(pbch_ce.get_row(0)))));
                syms = elem_mult(pbch_sym, gain);
                np = np_v(0) * sqr(gain);
            }
            else
            {
                syms.set_size(length(pbch_sym));
                np.set_size(length(pbch_sym));
#ifndef NDEBUG
                syms = NAN;
                np = NAN;
#endif
                for (int32 t = 0; t < length(syms); t += 2)
                {
                    // Simple zero-forcing
                    // http://en.wikipedia.org/wiki/Space-time_block_coding_based_transmit_diversity
                    complex<double> h1, h2;
                    double np_temp;
                    if (n_ports == 2)
                    {
                        h1 = (pbch_ce(0, t) + pbch_ce(0, t + 1)) / 2;
                        h2 = (pbch_ce(1, t) + pbch_ce(1, t + 1)) / 2;
                        np_temp = mean(np_v(0, 1));
                    }
                    else
                    {
                        if (mod(t, 4) == 0)
                        {
                            h1 = (pbch_ce(0, t) + pbch_ce(0, t + 1)) / 2;
                            h2 = (pbch_ce(2, t) + pbch_ce(2, t + 1)) / 2;
                            np_temp = (np_v(0) + np_v(2)) / 2;
                        }
                        else
                        {
                            h1 = (pbch_ce(1, t) + pbch_ce(1, t + 1)) / 2;
                            h2 = (pbch_ce(3, t) + pbch_ce(3, t + 1)) / 2;
                            np_temp = (np_v(1) + np_v(3)) / 2;
                        }
                    }
                    complex<double> x1 = pbch_sym(t);
                    complex<double> x2 = pbch_sym(t + 1);
                    double scale = pow(h1.real(), 2) + pow(h1.imag(), 2) + pow(h2.real(), 2) + pow(h2.imag(), 2);
                    syms(t) = (conj(h1) * x1 + h2 * conj(x2)) / scale;
                    syms(t + 1) = conj((-conj(h2) * x1 + h1 * conj(x2)) / scale);
                    np(t) = (pow(abs(h1) / scale, 2) + pow(abs(h2) / scale, 2)) * np_temp;
                    np(t + 1) = np(t);
                }
                // 3dB factor comes from precoding for transmit diversity
                syms = syms * pow(2, 0.5);
            }

            // Extract the bits from the complex modulated symbols.
            vec e_est = lte_demodulate(syms, np, modulation_t::QAM);
            // Unscramble
            bvec scr = lte_pn(cell.n_id_cell(), length(e_est));
            for (int32 t = 0; t < length(e_est); t++)
            {
                if (scr(t))
                    e_est(t) = -e_est(t);
            }
            // Undo ratematching
            mat d_est = lte_conv_deratematch(e_est, 40);
            // Decode
            c_est = lte_conv_decode(d_est);
            // Calculate received CRC
            bvec crc_est = lte_calc_crc(c_est(0, 23), CRC16);
            // Apply CRC mask
            if (n_ports == 2)
            {
                for (uint8 t = 0; t < 16; t++)
                {
                    crc_est(t) = 1 - (int)crc_est(t);
                }
            }
            else if (n_ports == 4)
            {
                for (uint8 t = 1; t < length(crc_est); t += 2)
                {
                    crc_est(t) = 1 - (int)crc_est(t);
                }
            }
            // Did we find it?
            if (crc_est == c_est(24, -1))
            {
                // YES!
                cell_out.n_ports = n_ports;
                // Unpack the MIB
                ivec c_est_ivec = to_ivec(c_est);
                // DL bandwidth
                const uint8 bw_packed = c_est_ivec(0) * 4 + c_est_ivec(1) * 2 + c_est_ivec(2);
                switch (bw_packed)
                {
                case 0:
                    cell_out.n_rb_dl = 6;
                    break;
                case 1:
                    cell_out.n_rb_dl = 15;
                    break;
                case 2:
                    cell_out.n_rb_dl = 25;
                    break;
                case 3:
                    cell_out.n_rb_dl = 50;
                    break;
                case 4:
                    cell_out.n_rb_dl = 75;
                    break;
                case 5:
                    cell_out.n_rb_dl = 100;
                    break;
                }
                // PHICH duration
                cell_out.phich_duration = c_est_ivec(3) ? phich_duration_t::EXTENDED : phich_duration_t::NORMAL;
                // PHICH resources
                uint8 phich_res = c_est_ivec(4) * 2 + c_est_ivec(5);
                switch (phich_res)
                {
                case 0:
                    cell_out.phich_resource = phich_resource_t::oneSixth;
                    break;
                case 1:
                    cell_out.phich_resource = phich_resource_t::half;
                    break;
                case 2:
                    cell_out.phich_resource = phich_resource_t::one;
                    break;
                case 3:
                    cell_out.phich_resource = phich_resource_t::two;
                    break;
                }
                // Calculate SFN
                int8 sfn_temp = 128 * c_est_ivec(6) + 64 * c_est_ivec(7) + 32 * c_est_ivec(8) + 16 * c_est_ivec(9) + 8 * c_est_ivec(10) + 4 * c_est_ivec(11) + 2 * c_est_ivec(12) + c_est_ivec(13);
                cell_out.sfn = itpp_ext::matlab_mod(sfn_temp * 4 - frame_timing_guess, 1024);
                return cell_out;
            }
        }
    }

    return cell_out;
}

// void tfoec(Cell &peak, cmat &tfg, vec &tfg_timestamp, double fc, int sampling_carrier_twist, int nRB, vec &tfg_comp_timestamp)
// {
//     int n_id_1 = peak.n_id_1;
//     int n_id_2 = peak.n_id_2;
//     cp_type_t::cp_type_t cp_type = peak.cp_type;

//     int decimation_ratio;
//     if (nRB == 6)
//     {
//         decimation_ratio = 16;
//     }
//     else if (nRB == 100)
//     {
//         decimation_ratio = 1;
//     }
//     double fs = 30720000.0 / decimation_ratio;
//     double fft_size = 2048 / decimation_ratio;
//     int nSC = nRB * 12;
//     int nRS = nRB * 2;
//     double k_factor;
//     if (sampling_carrier_twist == 1)
//     {
//         k_factor = (fc - peak.freq_fine) / fc;
//     }
//     else
//     {
//         k_factor = peak.k_factor;
//     }

//     int n_ofdm = tfg.rows();
//     int n_symb_dl;
//     double dft_location;
//     if (cp_type == cp_type_t::NORMAL)
//     {
//         n_symb_dl = 7;
//     }
//     else if (cp_type == cp_type_t::EXTENDED)
//     {
//         n_symb_dl = 6;
//     }
//     int n_id_cell = n_id_2 + 3 * n_id_1;
//     cmat rs0_start(20, nRS);
//     cmat rs0_mid(20, nRS);
//     int shift_start, shift_mid;
//     for (int slot_num = 0; slot_num < 20; slot_num++)
//     {
//         cvec r = rs_dl(slot_num, 0, 0, n_id_cell, nRB, cp_type, shift_start);
//         rs0_start.set_row(slot_num, r);
//         cvec r2 = rs_dl(slot_num, n_symb_dl - 3, 0, n_id_cell, nRB, cp_type, shift_mid);
//         rs0_mid.set_row(slot_num, r);
//     }

//     Sort<int> sort;
//     ivec rs_set = concat(itpp_ext::matlab_range(0, n_symb_dl, n_ofdm - 1), itpp_ext::matlab_range(n_symb_dl - 3, n_symb_dl, n_ofdm - 1));
//     sort.sort(0, rs_set.length() - 1, rs_set);
//     int slot_num = 0;
//     for (int t = 0; t < rs_set.length(); t++)
//     {
//         int idx = rs_set(t);
//         cvec tmp = tfg.get_row(idx);
//         if (t % 2 == 1)
//         {
//             cvec tt = conj(rs0_start.get_row(slot_num));
//             for (int j = shift_start; j < tmp.length(); t += 6)
//             {
//                 tmp[j] = tmp[j] * tt(j - shift_start);
//             }
//         }
//         else
//         {
//             cvec tt = conj(rs0_mid.get_row(slot_num));
//             for (int j = shift_mid; j < tmp.length(); t += 6)
//             {
//                 tmp[j] = tmp[j] * tt(j - shift_start);
//             }
//         }

//         if (t % 2 == 0)
//         {
//             slot_num = mod(slot_num, 20);
//         }
//     }

//     complex<double> foe(0, 0);
//     for (int t = 0; t < nRS; t++)
//     {
//         int idx = (t - 1) * 6 + shift_start + 1;
//         cvec rs_extracted = (tfg.get_col(idx)(itpp_ext::matlab_range(0, n_symb_dl, tfg.rows()))).transpose().get_row(0);
//         cvec rs_comp(rs_extracted);
//         foe += sum(elem_mult(conj(rs_comp(0, -2)), rs_comp(1, -1)));
//     }

//     for (int t = 0; t < nRS; t++)
//     {
//         int idx = (t - 1) * 6 + shift_mid + 1;
//         cvec rs_extracted = (tfg.get_col(idx)(itpp_ext::matlab_range(n_symb_dl - 3, n_symb_dl, tfg.rows()))).transpose().get_row(0);
//         cvec rs_comp(rs_extracted);
//         foe += sum(elem_mult(conj(rs_comp(0, -2)), rs_comp(1, -1)));
//     }

//     double residual_f = angle(foe) / (2 * pi) / (k_factor * .0005);
//     peak.freq_superfine = peak.freq_fine + residual_f;
//     if (sampling_carrier_twist)
//     {
//         peak.k_factor = (fc - peak.freq_superfine) / fc;
//     }
//     double k_factor_residual;
//     if (sampling_carrier_twist == 1)
//     {
//         k_factor_residual = (fc - residual_f) / fc;
//     }
//     else
//     {
//         k_factor_residual = 1;
//     }

//     cmat tfg_comp(tfg.rows(), tfg.cols());
//     tfg_comp_timestamp = (tfg_timestamp - 1) * k_factor_residual + 1;

// }