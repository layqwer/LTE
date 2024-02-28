#include <itpp/itbase.h>
#include <itpp/comm/convcode.h>
#include <itpp/comm/crc.h>
#include <itpp/signal/transforms.h>
#include <itpp/stat/misc_stat.h>
#include <itpp/signal/freq_filt.h>
#include <itpp/base/vec.h>
#include <itpp/signal/filter.h>
#include <list>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <tuple>
#include <string>
#include <filesystem>
#include <numeric>

#include "common.h"
#include "lte_lib.h"
#include "constants.h"
#include "macros.h"
#include "itpp_ext.h"
#include "dsp.h"
#include "utils.h"
#include "itpp1.h"


#include "stdafx.h"
//#include "itpp/itcomm.h"
//
//using std::cin;
//using std::cout;
//using std::endl;
//using namespace itpp;
//
// 当使用预编译的头时，需要使用此源文件，编译才能成功。
using namespace itpp;
using namespace std;

struct Peak
{
	int duplex_mode;
	int n_id_cell;
	int n_ports;
	double freq_superfine;
	double pow;
	std::string cp_type;
	int n_rb_dl;
	std::string phich_dur;
	int phich_res;
};

std::tuple<double, std::string> get_freq_hardware_from_filename(const std::string &filename)
{
	double fc = 0.0;
	std::string hardware;

	size_t sp = filename.find_last_of('\\');
	if (sp != std::string::npos)
	{
		std::string shortFilename = filename.substr(sp + 1);
		sp = shortFilename.find('f');
		if (sp == std::string::npos)
		{
			std::cout << "Can not find 'f'!" << std::endl;
			return{ fc, hardware };
		}

		size_t ep = shortFilename.find('_');
		if (ep == std::string::npos)
		{
			std::cout << "Can not find '_'!" << std::endl;
			return{ fc, hardware };
		}

		sp += 1;
		ep -= 1;
		fc = std::stod(shortFilename.substr(sp, ep - sp + 1)) * 1e6;

		if (shortFilename.find("rtlsdr") != std::string::npos)
		{
			hardware = "rtlsdr";
		}
		else if (shortFilename.find("hackrf") != std::string::npos)
		{
			hardware = "hackrf";
		}
		else if (shortFilename.find("bladerf") != std::string::npos)
		{
			hardware = "bladerf";
		}
		else if (shortFilename.find("usrp") != std::string::npos)
		{
			hardware = "usrp";
		}
		else if (shortFilename.find("rec") != std::string::npos)
		{
			hardware = "rec";
		}
		else
		{
			std::cout << "The filename does not have hardware information (rtlsdr/hackrf/bladerf/usrp)!" << std::endl;
		}
	}

	return{ fc, hardware };
}

cvec get_signal_from_bin(const std::string &filename, size_t num_sample_read, const std::string &dev)
{
	cvec s;
	std::ifstream file(filename, std::ios::binary);
	if (!file.is_open())
	{
		std::cout << "get_signal_from_bin: Can not open file!" << std::endl;
		return s;
	}

	if (dev == "hackrf")
	{
		// 获取文件大小
		file.seekg(0, std::ios::end);
		std::streampos fileSize = file.tellg();
		file.seekg(0, std::ios::beg);

		s = zeros_c(fileSize / 2);
		// 创建足够大的缓冲区
		vector<int8_t> buffer(fileSize);

		// 一次性读取整个文件
		file.read(reinterpret_cast<char *>(buffer.data()), fileSize);

		// 将数据转换为 complex<double>
		for (size_t i = 0; i < fileSize / 2; ++i)
		{
			s[i] = complex<double>(buffer[2 * i] / 128.0, buffer[2 * i + 1] / 128.0);
		}
	}
	// else if (dev == "rtlsdr") {
	// 	vector<size_t_t> buffer(num_sample_read * 2);
	// 	file.read(reinterpret_cast<char*>(buffer.data()), buffer.size());

	// 	for (size_t i = 0; i < num_sample_read; ++i)  {
	// 		s.push_back(raw2iq(buffer[2 * i], buffer[2     * i + 1]));
	// 	}
	// }
	// else if (dev == "bladerf" || dev == "usrp")
	// {
	// 	vector<size_t> buffer(num_sample_read * 2);
	// 	file.read(reinterpret_cast<char *>(buffer.data()), buffer.size());

	// 	for (size_t i = 0; i < num_sample_read; ++i)
	// 	{
	// 		s.ins(i, complex<double>(buffer[2 * i] / 32768.0, buffer[2 * i + 1] / 32768.0));
	// 	}
	// }
	// else if (dev == "rec")
	// {
	// 	size_t a;
	// 	size_t i = 0;
	// 	while (file.read(reinterpret_cast<char *>(&a), sizeof(a)))
	// 	{
	// 		s.ins(i++, complex<double>(a, 0.0));
	// 	}
	// }

	file.close();

	// if (num_sample_read != std::numeric_limits<size_t>::max() && s.size() != num_sample_read)
	// {
	// 	s.clear();
	// 	std::cout << "get_signal_from_bin: Not enough samples in the file!" << std::endl;
	// }

	return s;
}

cvec filter_wo_tail(cvec s, vec &coef, double convert_ratio)
{
	cvec r;
	if (s.size() == 1)
	{
		std::cerr << "filter_wo_tail: input signal must be column vector" << std::endl;
		return r;
	}
	size_t len_coef = coef.size();
	if (len_coef % 2 == 0)
	{
		std::cerr << "filter_wo_tail: length of coef must be odd!" << std::endl;
		return r;
	}

	if (convert_ratio <= 1)
	{
		int downsampling_ratio = static_cast<int>(std::round(1 / convert_ratio));
		cvec tmp(s);
		r = filter(coef, 1, concat(tmp, zeros_c(len_coef - 1)));
		r = r(((len_coef - 1) / 2), (r.size() - (len_coef - 1) / 2 - 1));
		// cout << "r " << r(4003) << " " << r(9994) << " " << r(30000) << " " << r(150000) << " " << r(153599) << endl;
		r = r(itpp_ext::matlab_range(0, downsampling_ratio, r.size() - 1));
		// cout << "r1 " << r(4003) << " " << r(9994) << " " << r(30000) << " " << r(150000) << " " << r(153599) << endl;
	}
	else
	{
		cvec upsampled_s = zeros_c(s.size() * convert_ratio);
		for (size_t i = 0; i < s.size(); ++i)
		{
			upsampled_s[i * static_cast<int>(convert_ratio)] = s[i];
		}

		cvec tmp(upsampled_s);
		r = filter(coef, 1, concat(tmp, zeros_c(len_coef - 1)));
		r = r(((len_coef - 1) / 2), (r.size() - (len_coef - 1) / 2 - 1));
	}
	return r;
}

// 等同[period_ppm, dynamic_f_search_set, xc, ~, ~, ~, ~, extra_info] =
// sampling_ppm_f_search_set_by_pss(capbuf_pbch.', f_search_set, td_pss, pss_fo_set, sampling_carrier_twist, pss_peak_max_reserve, num_pss_period_try, combined_pss_peak_range, par_th, num_peak_th);
// 返回4个值
void sampling_ppm_f_search_set_by_pss(
	// input
	cvec s, const vec &fo_search_set, cmat &td_pss, cmat pss_fo_set, int sampling_carrier_twist, complex<double> pss_peak_max_reserve,
	int num_pss_period_try, int combined_pss_peak_range, complex<double> par_th, int num_peak_th,
	// output
	size_t &ppm, vec &f_set, vector<mat> &xc)
{
	size_t len_pss = pss_fo_set.rows();
	size_t num_fo_pss = pss_fo_set.cols();

	size_t len = s.size();
	size_t len_short = len - (len_pss - 1);

	mat corr_store = fft_corr(s, td_pss, fo_search_set);

	if (sampling_carrier_twist == 1)
	{
		ppm = std::numeric_limits<size_t>::max();
		f_set = fo_search_set;
		ivec fo_idx_set = itpp_ext::matlab_range(0, 1, f_set.size() - 1);
		size_t n_f = f_set.size();
		xc[0] = mat(len_short, n_f);
		xc[1] = mat(len_short, n_f);
		xc[2] = mat(len_short, n_f);
		for (size_t foi = 0; foi < n_f; foi++)
		{
			for (size_t t = 0; t < 3; t++)
			{
				int col_idx = t * fo_search_set.size() + fo_idx_set[foi];
				for (int i = 0; i < len_short; ++i)
				{
					double tttt = corr_store(i, col_idx);
					xc[t](i, foi) = corr_store(i, col_idx);
				}
			}
		}
		return;
	}
}

void xcorr_pss(
	// input
	const cvec &capbuf, const vec &f_search_set, int DS_COMB_ARM, double &fc, int sampling_carrier_twist, double k_factor, vector<mat> &xc,
	// output
	mat &xc_incoherent_collapsed_pow, vector<ivec> &xc_incoherent_collapsed_frq, double &n_comb_xc, vec &sp_incoherent, vec &sp)
{
	int n_cap = capbuf.size();
	size_t n_f = f_search_set.size();
	sp = zeros(n_cap - 136 - 137);
	vec capbufx2 = absx2(capbuf);
	sp[0] = sum(capbufx2(0, 137 * 2 - 1));
	for (size_t k = 1; k < n_cap - 136 - 137 - 1; k++)
	{
		sp[k] = sp[k - 1] - capbufx2[k - 1] + capbufx2(k + 137 * 2 - 1);
	}
	sp /= (137.0 * 2.0);
	// cout << sp(0) << " " << sp(4999) << " " << sp(9500) << endl;
	vector<mat> xc_incoherent(3);
	xc_incoherent[0] = zeros(9600, n_f);
	xc_incoherent[1] = zeros(9600, n_f);
	xc_incoherent[2] = zeros(9600, n_f);

	vector<mat> xc_incoherent_single(3);
	xc_incoherent_single[0] = zeros(9600, n_f);
	xc_incoherent_single[1] = zeros(9600, n_f);
	xc_incoherent_single[2] = zeros(9600, n_f);

	n_comb_xc = floor((xc[0].rows() - 100) / 9600);

	for (size_t foi = 0; foi < n_f; foi++)
	{
		if (sampling_carrier_twist == 1)
		{
			double f_off = f_search_set[foi];
			k_factor = (fc - f_off) / fc;
		}
		for (size_t t = 0; t < 3; t++)
		{
			// xc_incoherent_single[t].set_col(foi, zeros(9600));
			for (size_t m = 0; m < n_comb_xc; m++)
			{

				// Because of the large supported frequency offsets and the large
				// amount of time represented by the capture buffer, the length
				// in samples, of a frame varies by the frequency offset.
				double actual_time_offset = m * 0.005 * k_factor;
				double actual_start_index = itpp::round_i(actual_time_offset * FS_LTE / 16);

				// double actual_start_index = itpp::round_i(m * .005 * k_factor * fs_programmed);
				for (size_t idx = 0; idx < 9600; idx++)
				{
					xc_incoherent_single[t](idx, foi) += xc[t](idx + actual_start_index, foi);
				}
				// cout << "xc_incoherent_single=" << xc_incoherent_single[t] << endl;
			}
			xc_incoherent_single[t].set_col(foi, xc_incoherent_single[t].get_col(foi) / n_comb_xc);
			xc_incoherent[t].set_col(foi, xc_incoherent_single[t].get_col(foi));
		}
		// cout << xc_incoherent_single[0](0, 0) << " " << xc_incoherent_single[0](4999, 0) << " " << xc_incoherent_single[0](9599, 0) << endl;
		// cout << xc_incoherent_single[1](0, 0) << " " << xc_incoherent_single[1](4999, 0) << " " << xc_incoherent_single[1](9599, 0) << endl;
		// cout << xc_incoherent_single[2](0, 0) << " " << xc_incoherent_single[2](4999, 0) << " " << xc_incoherent_single[2](9599, 0) << endl;
		for (int t = 0; t < DS_COMB_ARM; t++)
		{
			for (size_t k = 0; k < 3; k++)
			{
				vec t1(xc_incoherent_single[k].get_col(foi));
				vec t2(xc_incoherent_single[k].get_col(foi));
				t1 = tshiftmy(t1, -(t + 1));
				// cout << t1(0) << " " << t1(4999) << " " << t1(9500) << endl;
				t2 = tshiftmy(t2, t + 1);
				// cout << t2(0) << " " << t2(4999) << " " << t2(9500) << endl;
				for (size_t idx = 0; idx < 9600; idx++)
				{
					xc_incoherent[k](idx, foi) += (t1[idx] + t2[idx]);
				}
				// cout << "xc_incoherent=" << xc_incoherent[k](9599, foi) << " k=" << k << " idx=" << 9599 << " foi=" << foi << endl;
			}
		}
		for (size_t i = 0; i < 3; i++)
		{
			xc_incoherent[i].set_col(foi, xc_incoherent[i].get_col(foi) / (2 * DS_COMB_ARM + 1));
		}
	}
	// cout << "xc_incoherent " << xc_incoherent[0](0, 0) << " " << xc_incoherent[0](4999, 0) << " " << xc_incoherent[0](9599, 0) << endl;
	// cout << "xc_incoherent " << xc_incoherent[1](0, 0) << " " << xc_incoherent[1](4999, 0) << " " << xc_incoherent[1](9599, 0) << endl;
	// cout << "xc_incoherent " << xc_incoherent[2](0, 0) << " " << xc_incoherent[2](4999, 0) << " " << xc_incoherent[2](9599, 0) << endl;
	int n_comb_sp = floor(length(sp) / 9600);
	sp_incoherent = sum_by_columns(reshape(sp, 9600, n_comb_sp).transpose()) / n_comb_sp; // todo:
	sp_incoherent = tshiftmy(sp_incoherent, 137);
	// cout << sp_incoherent(0) << sp_incoherent(1000) << sp_incoherent(6000) << sp_incoherent(9599) << endl;
	for (size_t t = 0; t < 3; t++)
	{
		xc_incoherent_collapsed_pow.set_row(t, itpp::max(transpose(xc_incoherent[t]), xc_incoherent_collapsed_frq[t]));
		vec pow = xc_incoherent_collapsed_pow.get_row(t);
		ivec frq = xc_incoherent_collapsed_frq[t];
		// cout << pow(0) << " " << pow(4999) << " " << pow(9500) << endl;
		// cout << frq(0) << " " << frq(4999) << " " << frq(9500) << endl;
	}
}

mat myudb10(mat &db)
{
	mat res = zeros(db.rows(), db.cols());
	for (size_t row = 0; row < db.rows(); row++)
	{
		for (size_t col = 0; col < db.cols(); col++)
		{
			res.set(row, col, std::pow(10, db.get(row, col) / 10.0));
		}
	}
	return res;
}

double myudb10(double s)
{
	return std::pow(10, s / 10.0);
}

void peak_search(
	// Inputs
	mat &xc_incoherent_collapsed_pow,
	vector<ivec> &xc_incoherent_collapsed_frq,
	const vec &Z_th1,
	const vec &f_search_set,
	double fc,
	int sampling_carrier_twist,

	// Outputs
	list<Cell> &cells)
{
	mat xc_incoherent_working(xc_incoherent_collapsed_pow);

	for (;;)
	{
		// Search for the largest peak. (Not the largest peak relative to
		// the detection threshold Z_th1.)
		ivec peak_ind_v;
		vec peak_pow_v = max(transpose(xc_incoherent_working), peak_ind_v);
		int32 peak_n_id_2;
		double peak_pow = max(peak_pow_v, peak_n_id_2);
		int32 peak_ind = peak_ind_v(peak_n_id_2);
		if (peak_pow < Z_th1(peak_ind))
		{
			// This peak has too low of a received power. There are no more
			// interesting peaks. Break!
			break;
		}
		//    cout << "peak_search " << peak_pow << " " << peak_ind << " " << Z_th1(peak_ind) << "\n";

		// A peak was found at location peak_ind and has frequency index
		// xc_incoherent_collapsed_frq(peak_n_id_2,peak_ind). This peak
		// is the sum of the energy within ds_comb_arm samples around this
		// peak location. From the samples within ds_comb_arm samples
		// around peak_ind, find the index with the highest power.
		// double best_pow=-INFINITY;
		// int16 best_ind=-1;
		// for (uint16 t=peak_ind-DS_COMB_ARM;t<=peak_ind+DS_COMB_ARM;t++) {
		//   uint16 t_wrap=mod(t,9600);
		//   if (xc_incoherent_single[peak_n_id_2](xc_incoherent_collapsed_frq(peak_n_id_2,peak_ind),t_wrap)>best_pow) {
		//     best_pow=xc_incoherent_single[peak_n_id_2](xc_incoherent_collapsed_frq(peak_n_id_2,peak_ind),t_wrap);
		//     best_ind=t_wrap;
		//   }
		// }

		// Record this peak for further processing
		Cell cell;
		// cell.fc_requested=fc_requested;
		// cell.fc_programmed=fc_programmed;
		cell.pss_pow = peak_pow;
		cell.ind = peak_ind;
		// cell.ind=best_ind;
		cell.freq = f_search_set(xc_incoherent_collapsed_frq[peak_n_id_2].get(peak_ind));
		cell.n_id_2 = peak_n_id_2;

		if (sampling_carrier_twist)
		{
			//      cell.k_factor = (fc_requested-cell.freq)/fc_programmed;
			cell.k_factor = (fc - cell.freq) / fc;
		}
		else
		{
			//   cell.k_factor = k_factor;
		}

		cells.push_back(cell); // for tdd test
		cells.push_back(cell); // for fdd test

		ivec cancel = wrap(itpp_ext::matlab_range(peak_ind - 137 * 2, peak_ind + 137 * 2), 0, 9600);
		for (size_t i = 0; i < cancel.size(); i++)
		{
			xc_incoherent_working.set(peak_n_id_2, cancel[i], 0);
		}
		mat t("0 -8 -8;-8 0 -8;-8 -8 0");
		mat pss_xc_matrix = myudb10(t);
		for (int t = 0; t <= 2; ++t)
		{
			if (t != peak_n_id_2)
			{
				double tmp = peak_pow * pss_xc_matrix(peak_n_id_2, t);
				for (size_t i = 0; i < cancel.size(); i++)
				{
					if (xc_incoherent_working.get(t, cancel[i]) < tmp)
					{
						xc_incoherent_working.set(t, cancel[i], 0);
					}
				}
			}
		}
		double tmp = peak_pow * myudb10(-12);
		for (int r = 0; r < xc_incoherent_working.rows(); r++)
		{
			for (int c = 0; c < xc_incoherent_working.rows(); c++)
			{
				if (xc_incoherent_working.get(r, c) < tmp)
				{
					xc_incoherent_working.set(r, c, 0);
				}
			}
		}
	}
}

	string CellSearch(const cvec &r_pbch, const cvec &r_20M,
		const vec &f_search_set, double fc,
		int sampling_carrier_twist, complex<double> pss_peak_max_reserve,
		int num_pss_period_try, int combined_pss_peak_range,
		complex<double> par_th, int num_peak_th)
	{
		cmat td_pss(128 + 9, 3);
		pss_gen(td_pss);
		cmat pss_fo_set = pss_fo_set_gen(td_pss, f_search_set);

		int num_radioframe = 8;
		int pbch_sampling_ratio = 16;
		double sampling_rate = 30.72e6;
		double sampling_rate_pbch = sampling_rate / pbch_sampling_ratio;
		int num_subframe_per_radioframe = 10;
		double len_time_subframe = 1e-3;
		double num_sample_per_radioframe = num_subframe_per_radioframe * len_time_subframe * sampling_rate_pbch;
		double num_sample_pbch = num_radioframe * num_sample_per_radioframe;

		int DS_COMB_ARM = 2;
		int thresh1_n_nines = 12;
		double rx_cutoff = (6 * 12 * 15e3 / 2 + 4 * 15e3) / (FS_LTE / 16 / 2);
		int THRESH2_N_SIGMA = 3;
		int ex_gain = 2;

		size_t sp = 0;
		size_t ep = sp + num_sample_pbch - 1;

		cvec capbuf_pbch = r_pbch(sp, ep);
		// cout << "capbuf_pbch1 " << capbuf_pbch(4003) << " " << capbuf_pbch(9994) << " " << capbuf_pbch(30000) << " " << capbuf_pbch(150000) << " " << capbuf_pbch(153599) << endl;

		// size_t sp_20M = (sp - 1) * pbch_sampling_ratio;
		// size_t ep_20M = ep * pbch_sampling_ratio - 1;
		capbuf_pbch = capbuf_pbch - mean(capbuf_pbch);
		// cout << "capbuf_pbch2 " << capbuf_pbch(4003) << " " << capbuf_pbch(9994) << " " << capbuf_pbch(30000) << " " << capbuf_pbch(150000) << " " << capbuf_pbch(153599) << endl;
		// cvec r_pbch_sub = capbuf_pbch;

		// 只是返回，可以没有
		// cvec r_20M_sub;
		// if (!r_20M.length() == 0)
		// {
		// 	r_20M_sub = r_20M(sp_20M, ep_20M);
		// }

		// period_ppm, dynamic_f_search_set, xc, ~, ~, ~, ~, extra_info
		size_t period_ppm;
		vec dynamic_f_search_set;
		vector<mat> xc(3);
		// Peak extra_info[];
		// std::cout << "pss_fo_set" << pss_fo_set(0, 0) << pss_fo_set(pss_fo_set.rows() - 1, 0) << pss_fo_set(0, pss_fo_set.cols() - 1) << pss_fo_set(pss_fo_set.rows() - 1, pss_fo_set.cols() - 1);
		sampling_ppm_f_search_set_by_pss(capbuf_pbch, f_search_set, td_pss, pss_fo_set, sampling_carrier_twist, pss_peak_max_reserve, num_pss_period_try, combined_pss_peak_range, par_th, num_peak_th, period_ppm, dynamic_f_search_set, xc);
		double test = dynamic_f_search_set[dynamic_f_search_set.length() - 1];

		mat xc_incoherent_collapsed_pow(3, 9600);
		vector<ivec> xc_incoherent_collapsed_frq(3);
		xc_incoherent_collapsed_frq[0] = zeros_i(9600);
		xc_incoherent_collapsed_frq[1] = zeros_i(9600);
		xc_incoherent_collapsed_frq[2] = zeros_i(9600);
		double n_comb_xc;
		vec sp_incoherent(9600);
		vec sp_vec(9600);
		xcorr_pss(capbuf_pbch, f_search_set, DS_COMB_ARM, fc, sampling_carrier_twist, 0.0, xc, xc_incoherent_collapsed_pow, xc_incoherent_collapsed_frq, n_comb_xc, sp_incoherent, sp_vec);

		double R_th1 = chi2cdf_inv(1 - pow(10.0, -thresh1_n_nines), 2 * n_comb_xc * (2 * DS_COMB_ARM + 1));
		vec Z_th1 = R_th1 * sp_incoherent / rx_cutoff / 137 / n_comb_xc / (2 * DS_COMB_ARM + 1); // remove /2 to avoid many false alarm

		list<Cell> peaks;
		peak_search(xc_incoherent_collapsed_pow, xc_incoherent_collapsed_frq, Z_th1, dynamic_f_search_set, fc, sampling_carrier_twist, peaks);

		mat t1(1, peaks.size() / 2);
		t1.ones();
		mat t2("0 1");
		vec tdd_flags = kron(t1, t2).get_row(0);
		vec detect_flag(peaks.size());
		detect_flag.zeros();
		size_t i = 0;
		for (auto &peak : peaks)
		{
			double tdd_flag = tdd_flags.get(i);
			sss_detect(peak, capbuf_pbch, THRESH2_N_SIGMA, fc, sampling_carrier_twist, tdd_flag);
			if (peak.n_id_1 != 0)
			{
				pss_sss_foe(peak, capbuf_pbch, fc, sampling_carrier_twist, tdd_flag);
				cmat tfg;
				vec tfg_timestamp;
				extract_tfg(peak, capbuf_pbch, fc, sampling_carrier_twist, 6, tfg, tfg_timestamp);
				// cout << tfg_timestamp << endl;
				cmat tfg_comp;
				vec tfg_comp_timestamp;
				RS_DL rs_dl(peak.n_id_cell(), 6, peak.cp_type);
				peak = tfoec(peak, tfg, tfg_timestamp, fc, sampling_carrier_twist, rs_dl, tfg_comp, tfg_comp_timestamp);
				// cout << tfg_comp_timestamp << endl;
				peak = decode_mib_2(peak, tfg, rs_dl);
				if (peak.n_rb_dl != 0)
				{
					continue;
				}
				else if (tdd_flag == 1)
				{
					cout << "  Detected a TDD cell!" << endl;
				}
				else
				{
					cout << "    cell ID: " << peak.n_id_cell() << endl;
				}
				cout << "    PSS  ID: " << (peak.n_id_2 + 1) << endl;
				cout << "    RX power level: " << (10 * log10(peak.pss_pow)) << endl;
				cout << "    residual frequency offset: " << (peak.freq_superfine) << endl;
				detect_flag(i) = 1;
			}
		}

		if (sum(detect_flag) == 0)
		{
			cout << "No LTE cells were found..." << endl;
		}
		else
		{
			return "sum(detect_flag) != 0";
		}

		std::stringstream ss;
		ss << "-------------------------------Cells information summary-------------" << peaks.size() << std::endl;
		for (const auto &peak : peaks)
		{
			ss << "Cell " << i++ << " information:--------------------------------------------------------" << std::endl;
			ss << "            Cell mode: " << (peak.duplex_mode == 1 ? "TDD" : "FDD") << std::endl;
			ss << "              Cell ID: " << peak.n_id_cell() << endl;
			ss << "   Num. eNB Ant ports: " << peak.n_ports << endl;
			ss << "    Carrier frequency: " << fc / 1e6 << "MHz" << endl;
			ss << "Residual freq. offset: " << peak.freq_superfine / 1e3 << "kHz" << endl;
			ss << "       RX power level: " << 10 * log10(peak.pss_pow) << std::endl;
			ss << "              CP type: " << peak.cp_type << endl;
			ss << "              Num. RB: " << peak.n_rb_dl << endl;
			ss << "       PHICH duration: " << peak.phich_duration << endl;
			ss << "  PHICH resource type: " << peak.phich_resource << endl;
		}

		return ss.str();
	}

// __declspec(dllexport) const char* search(char* s)
// {
// 	std::stringstream ss;
// 	ss << "-------------------------------Cells information summary-------------" << std::endl;
// 	const char* return_char;
// 	return_char = strdup(ss.str().c_str());
// 	return return_char;
// 	std::cout << "123123" << '\n';
// }

__declspec(dllexport) const char* search(char* s)
{
	 // std::string filename = "C:\\Users\\yue\\Documents\\Visual Studio 2015\\Projects\\itpp1\\Win32\\Debug\\f1815.3_s19.2_bw20_0.08s_hackrf-1.bin";
	 // double gain1 = -1;
	 // double gain2 = -1;
	std::string filename(s);
	 int sampling_carrier_twist = 1;
	 double num_radioframe = 8;
	 double num_second = num_radioframe * 10e-3;
	 double raw_sampling_rate = 19.2e6;
	 int nRB = 100;
	 double sampling_rate = 30.72e6;
	 double sampling_rate_pbch = sampling_rate / 16; // LTE spec. 30.72MHz/16.
	 double bandwidth = 20e6;
	
	 int pss_peak_max_reserve = 2;
	 int num_pss_period_try = 1;
	 int combined_pss_peak_range = -1;
	 double par_th = 8.5;
	 double num_peak_th = 1 / 2;
	
	 double fc;
	 std::string sdr_board = "hackrf";
	
	 // 获取频率和硬件信息
	 std::tie(fc, sdr_board) = get_freq_hardware_from_filename(filename);
	
	 if (fc == 0.0 || sdr_board.empty())
	 {
	 	std::cout << filename << " does not include valid frequency or hardware info!" << std::endl;
	 	return "Does not include valid frequency or hardware info!"; // 返回非零表示错误
	 }
	
	 // 获取信号数据
	 cvec r_raw = get_signal_from_bin(filename, std::numeric_limits<size_t>::max(), sdr_board);
	
	 // 处理 rtlsdr 特定参数
	 if (sdr_board == "rtlsdr")
	 {
	 	raw_sampling_rate = 1.92e6;
	 	nRB = 6; // PBCH only, for rtlsdr
	 			 // 进一步处理 rtlsdr 相关逻辑
	 }
	
	 // std::cout << "fc " << fc << "; IQ from " << sdr_board << " " << filename << std::endl;
	 vec coef_pbch = fir1(254, (0.18e6 * 6 + 150e3) / raw_sampling_rate);
	 vec coef_8x_up = fir1(254, 20e6 / (raw_sampling_rate * 8));
	
	 vec f_search_set = "-140e3:5e3:135e3";
	
	 r_raw = r_raw - mean(r_raw);
	 // cout << "r_raw " << r_raw(4003) << " " << r_raw(9994) << " " << r_raw(30000) << " " << r_raw(150000) << " " << r_raw(153599) << endl;
	 cvec r_pbch, r_20M;
	 if (sdr_board == "rtlsdr")
	 {
	 	r_pbch = r_raw;
	 	r_20M.clear();
	 }
	 else
	 {
	 	vec t1 = coef_pbch * 5;
	 	vec t2 = coef_8x_up * 8;
	 	r_pbch = filter_wo_tail(r_raw(0, 80e-3 * raw_sampling_rate - 1), t1, sampling_rate_pbch / raw_sampling_rate);
	 	r_20M = filter_wo_tail(r_raw(0, 80e-3 * raw_sampling_rate - 1), t2, 8);
	 	ivec indices = itpp_ext::matlab_range(0, 5, r_20M.size());
	 	// ivec indices = to_ivec(linspace(0, r_20M.size(), r_20M.size() / 5));
	 	r_20M = r_20M((indices));
	 }
	 // cout << "r_pbch " << r_pbch(4003) << " " << r_pbch(9994) << " " << r_pbch(30000) << " " << r_pbch(150000) << " " << r_pbch(153599) << endl;
	 //
	 // 可以不返回东西, 对应[cell_info, r_pbch, r_20M] = CellSearch(r_pbch, r_20M, f_search_set, fc, sampling_carrier_twist, pss_peak_max_reserve, num_pss_period_try, combined_pss_peak_range, par_th, num_peak_th);
	 //  std::string cell_info = CellSearch(r_pbch, r_20M, f_search_set, fc, 1, 2, 1, -1, 8.5, 0.5);

	 string res = CellSearch(r_pbch, r_20M, f_search_set, fc, sampling_carrier_twist, pss_peak_max_reserve,
	                         num_pss_period_try, combined_pss_peak_range, par_th, num_peak_th);
	 const char* return_char;
	 return_char = strdup(res.c_str());
	 return return_char;
}
