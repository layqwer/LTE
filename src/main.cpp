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

// 当使用预编译的头时，需要使用此源文件，编译才能成功。
using namespace itpp;
using namespace std;
namespace mynamespace
{
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

		size_t sp = filename.find('/');
		if (sp != std::string::npos)
		{
			std::string shortFilename = filename.substr(sp + 1);
			sp = shortFilename.find('f');
			if (sp == std::string::npos)
			{
				std::cout << "Can not find 'f'!" << std::endl;
				return {fc, hardware};
			}

			size_t ep = shortFilename.find('_');
			if (ep == std::string::npos)
			{
				std::cout << "Can not find '_'!" << std::endl;
				return {fc, hardware};
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

		return {fc, hardware};
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

			// 创建足够大的缓冲区
			vector<int8_t> buffer(fileSize);

			// 一次性读取整个文件
			file.read(reinterpret_cast<char *>(buffer.data()), fileSize);
			file.close();

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
		else if (dev == "bladerf" || dev == "usrp")
		{
			vector<size_t> buffer(num_sample_read * 2);
			file.read(reinterpret_cast<char *>(buffer.data()), buffer.size());

			for (size_t i = 0; i < num_sample_read; ++i)
			{
				s.ins(i, complex<double>(buffer[2 * i] / 32768.0, buffer[2 * i + 1] / 32768.0));
			}
		}
		else if (dev == "rec")
		{
			size_t a;
			size_t i = 0;
			while (file.read(reinterpret_cast<char *>(&a), sizeof(a)))
			{
				s.ins(i++, complex<double>(a, 0.0));
			}
		}

		file.close();

		// if (num_sample_read != std::numeric_limits<size_t>::max() && s.size() != num_sample_read)
		// {
		// 	s.clear();
		// 	std::cout << "get_signal_from_bin: Not enough samples in the file!" << std::endl;
		// }

		return s;
	}

	cvec filter_wo_tail(cvec s, vec &&coef, double convert_ratio)
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
			r = r((((len_coef - 1) / 2) + 1), (r.size() - ((len_coef - 1) / 2)));
			r = r(to_ivec(linspace(1, r.size(), downsampling_ratio)));
		}
		else
		{
			// Upsampling
			cvec upsampled_s;
			upsampled_s.set_size(s.size() * convert_ratio);
			upsampled_s.zeros();
			for (size_t i = 0; i < s.size(); ++i)
			{
				upsampled_s[i * static_cast<int>(convert_ratio)] = s[i];
			}

			cvec tmp(s);
			r = filter(coef, 1, concat(tmp, zeros_c(len_coef - 1)));
			r = r((((len_coef - 1) / 2) + 1), (r.size() - ((len_coef - 1) / 2)));
		}

		return r;
	}

	// 等同[period_ppm, dynamic_f_search_set, xc, ~, ~, ~, ~, extra_info] =
	// sampling_ppm_f_search_set_by_pss(capbuf_pbch.', f_search_set, td_pss, pss_fo_set, sampling_carrier_twist, pss_peak_max_reserve, num_pss_period_try, combined_pss_peak_range, par_th, num_peak_th);
	// 返回4个值
	void sampling_ppm_f_search_set_by_pss(
		// input
		cvec s, const vec &fo_search_set, cmat & td_pss, cmat pss_fo_set, int sampling_carrier_twist, complex<double> pss_peak_max_reserve,
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
			vec fo_idx_set = itpp::linspace(0, f_set.size() - 1, 1);
			size_t n_f = f_set.size();
			xc[0] = mat(n_f, len_short);
			xc[1] = mat(n_f, len_short);
			xc[2] = mat(n_f, len_short);
			for (size_t foi = 0; foi < n_f; foi++)
			{
				for (size_t t = 0; t < 3; t++)
				{
					int col_idx = (t - 1) * fo_search_set.size() + fo_idx_set[foi];
					for (int i = 0; i < len_short; ++i)
					{
						xc[t](i, foi) = corr_store(i, col_idx);
					}
				}
			}
			// 不需要
			// pss_idx_set = std::numeric_limits<size_t>::max();
			// fo_pss_idx_set = std::numeric_limits<size_t>::max();
			// fo_with_all_pss_idx = std::numeric_limits<size_t>::max();
			return;
		}
	}

	// todo: sp_incoherent有没有复数
	void xcorr_pss(
		// input
		const cvec &capbuf, const vec &f_search_set, int DS_COMB_ARM, double &fc, int sampling_carrier_twist, double k_factor, vector<mat> &xc,
		// output
		mat &xc_incoherent_collapsed_pow, vector<ivec> &xc_incoherent_collapsed_frq, double &n_comb_xc, vec &sp_incoherent, vec &sp)
	{
		int n_cap = capbuf.size();
		size_t n_f = f_search_set.size();
		sp = zeros(n_cap - 136 - 137 - 1);
		vec capbufx2 = absx2(capbuf);
		sp[0] = sum(capbufx2(0, 137 * 2 - 1));
		for (size_t k = 1; k < n_cap - 136 - 137 - 1; k++)
		{
			sp[k] = sp[k - 1] - capbufx2[k - 1] + capbufx2(k + 137 * 2 - 2);
		}
		sp /= (137.0 * 2.0);

		vector<mat> xc_incoherent;
		xc_incoherent[0].set_size(n_f, 9600);
		xc_incoherent[1].set_size(n_f, 9600);
		xc_incoherent[2].set_size(n_f, 9600);

		vector<mat> xc_incoherent_single;
		xc_incoherent_single[0].set_size(n_f, 9600);
		xc_incoherent_single[1].set_size(n_f, 9600);
		xc_incoherent_single[2].set_size(n_f, 9600);

		n_comb_xc = floor((xc[0].size() - 100) / 9600);
		for (size_t foi = 0; foi < n_f; foi++)
		{
			if (sampling_carrier_twist == 1)
			{
				double f_off = f_search_set[foi];
				k_factor = (fc - f_off) / fc;
			}
			for (size_t t = 0; t < 3; t++)
			{
				xc_incoherent_single[t].set_row(foi, zeros(9600));
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
						xc_incoherent_single[t](foi, idx) += xc[t](foi, idx + actual_start_index);
					}
				}
				xc_incoherent_single[t].set_row(foi, xc_incoherent_single[t].get_row(foi) / n_comb_xc);
				xc_incoherent[t].set_row(foi, xc_incoherent_single[t].get_row(foi));
			}
			for (int t = 0; t < DS_COMB_ARM; t++)
			{
				for (size_t k = 0; k < 3; k++)
				{
					vec t1(xc_incoherent_single[k].get_row(foi));
					vec t2(xc_incoherent_single[k].get_row(foi));
					tshift(t1, -t);
					tshift(t2, t);
					for (size_t idx = 0; idx < 9600; idx++)
					{
						xc_incoherent[k](foi, idx) += (t1[idx] + t2[idx]);
					}
				}
				xc_incoherent[t].set_row(foi, xc_incoherent[t].get_row(foi) / (2 * DS_COMB_ARM + 1));
			}
		}
		// sp_incoherent = zeros_c(9600);
		int n_comb_sp = floor(length(sp) / 9600);
		sp_incoherent = sum_by_columns(reshape(sp, 9600, n_comb_sp).transpose()) / n_comb_sp; // todo:
		tshift(sp_incoherent, 137);

		for (size_t t = 0; t < 3; t++)
		{
			xc_incoherent_collapsed_pow.set_row(t, itpp::max(xc_incoherent[t], xc_incoherent_collapsed_frq[t]));
		}
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
		// pss_xc_matrix=udb10([NaN -8 -8; -8 NaN -8; -8 -8 NaN]);//todo:
		// Create local copy we can write to and destroy.
		mat xc_incoherent_working(xc_incoherent_collapsed_pow);

		for (;;)
		{
			// Search for the largest peak. (Not the largest peak relative to
			// the detection threshold Z_th1.)
			ivec peak_ind_v;
			vec peak_pow_v = max(transpose(xc_incoherent_working), peak_ind_v);
			int peak_n_id_2;
			double peak_pow = max(peak_pow_v, peak_n_id_2);
			int peak_ind = peak_ind_v(peak_n_id_2);
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
			// double best_pow = -INFINITY;
			// size_t best_ind = -1;
			// for (size_t t = peak_ind - 2; t <= peak_ind + 2; t++)
			// {
			// 	size_t t_wrap = mod(t, 9600);
			// 	if (xc_incoherent_single[peak_n_id_2](xc_incoherent_collapsed_frq(peak_n_id_2, peak_ind), t_wrap) > best_pow)
			// 	{
			// 		best_pow = xc_incoherent_single[peak_n_id_2](xc_incoherent_collapsed_frq(peak_n_id_2, peak_ind), t_wrap);
			// 		best_ind = t_wrap;
			// 	}
			// }

			// Record this peak for further processing
			Cell cell;
			// cell.fc_requested = fc_requested;
			// cell.fc_programmed = fc_programmed;
			cell.pss_pow = peak_pow;
			cell.ind = peak_ind;
			// cell.ind = best_ind;
			// cell.freq = f_search_set(xc_incoherent_collapsed_frq(peak_n_id_2, peak_ind));
			cell.n_id_2 = peak_n_id_2;

			// if (sampling_carrier_twist)
			// {
			// 	//      cell.k_factor = (fc_requested-cell.freq)/fc_programmed;
			// 	cell.k_factor = (fc_programmed - cell.freq) / fc_programmed;
			// }
			// else
			// {
			// 	cell.k_factor = k_factor;
			// }

			cells.push_back(cell); // for tdd test
			cells.push_back(cell); // for fdd test

			// Cancel out the false peaks around this one.
			// No peaks with the same pss sequence are allowed within 274 samples of
			// this one.
			for (int t = -274; t <= 274; t++)
			{
				// cout <<mod(peak_ind+t,9600)<<endl;
				xc_incoherent_working(peak_n_id_2, itpp_ext::matlab_mod(peak_ind + t, 9600)) = 0;
			}
			// Cancel out other PSS sequences whose power is within 8dB of the current
			// sequence.
			double thresh = peak_pow * udb10(-8.0);
			for (size_t n = 0; n <= 3; n++)
			{
				if (n == peak_n_id_2)
				{
					continue;
				}
				for (int t = -274; t <= 274; t++)
				{
					if (xc_incoherent_working(peak_n_id_2, itpp_ext::matlab_mod(peak_ind + t, 9600)) < thresh)
					{
						xc_incoherent_working(peak_n_id_2, itpp_ext::matlab_mod(peak_ind + t, 9600)) = 0;
					}
				}
			}
			// Because of the repetitive nature of the CRS, a PSS at offset I with power
			// P will produce correlation peaks during all CRS OFDM symbols with power
			// approximately P-14dB. Cancel them out.
			thresh = peak_pow * udb10(-12.0);
			for (size_t r = 0; r < 3; r++)
			{
				for (size_t c = 0; c < 9600; c++)
				{
					if (xc_incoherent_working(r, c) < thresh)
					{
						xc_incoherent_working(r, c) = 0;
					}
				}
			}
		}
	}

	void CellSearch(const cvec &r_pbch, const cvec &r_20M,
											   const vec &f_search_set, double fc,
											   int sampling_carrier_twist, complex<double> pss_peak_max_reserve,
											   int num_pss_period_try, int combined_pss_peak_range,
											   complex<double> par_th, int num_peak_th)
	{
		cmat td_pss(128 + 9, 3);
		pss_gen(td_pss);
		cmat pss_fo_set = pss_fo_set_gen(td_pss, f_search_set);
		// pss_fo_set_gen(f_search_set, pss_fo_set);

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

		size_t sp = 1;
		size_t ep = sp + num_sample_pbch - 1;

		cvec capbuf_pbch = r_pbch(sp, ep).transpose().get_row(0);

		size_t sp_20M = (sp - 1) * pbch_sampling_ratio;
		size_t ep_20M = ep * pbch_sampling_ratio - 1;
		capbuf_pbch = capbuf_pbch - mean(capbuf_pbch);
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
		sampling_ppm_f_search_set_by_pss(capbuf_pbch.transpose().get_row(0), f_search_set, td_pss, pss_fo_set, sampling_carrier_twist, pss_peak_max_reserve, num_pss_period_try, combined_pss_peak_range, par_th, num_peak_th, period_ppm, dynamic_f_search_set, xc);

		// [xc_incoherent_collapsed_pow, xc_incoherent_collapsed_frq, n_comb_xc, ~, ~, ~, sp_incoherent, sp]= ...
		// xcorr_pss(capbuf_pbch,dynamic_f_search_set,DS_COMB_ARM,fc,sampling_carrier_twist,NaN, xc);
		mat xc_incoherent_collapsed_pow;
		vector<ivec> xc_incoherent_collapsed_frq(3);
		double n_comb_xc;
		vec sp_incoherent(9600);
		vec sp_vec(9600);
		xcorr_pss(capbuf_pbch, f_search_set, DS_COMB_ARM, fc, sampling_carrier_twist, 0.0, xc, xc_incoherent_collapsed_pow, xc_incoherent_collapsed_frq, n_comb_xc, sp_incoherent, sp_vec);

		double R_th1 = chi2cdf_inv(1 - pow(10.0, -thresh1_n_nines), 2 * n_comb_xc * (2 * DS_COMB_ARM + 1));
		vec Z_th1 = R_th1 * sp_incoherent / rx_cutoff / 137 / n_comb_xc / (2 * DS_COMB_ARM + 1); // remove /2 to avoid many false alarm

		list<Cell> cells;
		peak_search(xc_incoherent_collapsed_pow, xc_incoherent_collapsed_frq, Z_th1, dynamic_f_search_set, fc, sampling_carrier_twist, cells);

		std::cout << "-------------------------------Cells information summary-------------------------------" << std::endl;
		int i = 1;
		for (const auto &value : cells)
		{
			std::cout << "Cell " << i++ << " information:--------------------------------------------------------" << std::endl;
			std::cout << "Cell mode: " << (value.duplex_mode == 1 ? "TDD" : "FDD") << std::endl;
			std::cout << "RX power level: " << 10 * log10(value.pss_pow) << std::endl;
		}
	}

	// std::string lteDLReceiver() {
	int main(int arg, char* args[])
	{
		std::string filename = "../regression_test_signal_file/f1815.3_s19.2_bw20_0.08s_hackrf-1.bin";
		double gain1 = -1;
		double gain2 = -1;
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
		std::string sdr_board;

		// 获取频率和硬件信息
		std::tie(fc, sdr_board) = get_freq_hardware_from_filename(filename);

		if (fc == 0.0 || sdr_board.empty())
		{
			std::cout << filename << " does not include valid frequency or hardware info!" << std::endl;
			return 1; // 返回非零表示错误
		}

		std::cout << filename << std::endl;

		// 获取信号数据
		cvec r_raw = get_signal_from_bin(filename, std::numeric_limits<size_t>::max(), sdr_board);

		// 处理 rtlsdr 特定参数
		if (sdr_board == "rtlsdr")
		{
			raw_sampling_rate = 1.92e6;
			nRB = 6; // PBCH only, for rtlsdr
					 // 进一步处理 rtlsdr 相关逻辑
		}

		std::cout << "fc " << fc << "; IQ from " << sdr_board << " " << filename << std::endl;
		vec coef_pbch = fir1(254, (0.18e6 * 6 + 150e3) / raw_sampling_rate);
		vec coef_8x_up = fir1(254, 20e6 / (raw_sampling_rate * 8));

		vec f_search_set = "-140e3:5e3:135e3";

		r_raw = r_raw - mean(r_raw);
		cvec r_pbch, r_20M;
		if (sdr_board == "rtlsdr")
		{
			r_pbch = r_raw;
			r_20M.clear();
		}
		else
		{
			r_pbch = filter_wo_tail(r_raw(0, 80e-3 * raw_sampling_rate - 1), coef_pbch * 5, sampling_rate_pbch / raw_sampling_rate);
			r_20M = filter_wo_tail(r_raw(0, 80e-3 * raw_sampling_rate - 1), coef_8x_up * 8, 8);
			vec indices = linspace(0, r_20M.size() - 1, 5);
			r_20M = r_20M(to_ivec(indices));
		}

		// 可以不返回东西, 对应[cell_info, r_pbch, r_20M] = CellSearch(r_pbch, r_20M, f_search_set, fc, sampling_carrier_twist, pss_peak_max_reserve, num_pss_period_try, combined_pss_peak_range, par_th, num_peak_th);
		//  std::string cell_info = CellSearch(r_pbch, r_20M, f_search_set, fc, 1, 2, 1, -1, 8.5, 0.5);
		CellSearch(r_pbch, r_20M, f_search_set, fc, sampling_carrier_twist, pss_peak_max_reserve, num_pss_period_try, combined_pss_peak_range, par_th, num_peak_th);

		return 0;
	}

}
